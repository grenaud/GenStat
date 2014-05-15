/*
 * ComputeDivergence
 * Date: Aug-15-2012 
 * Author : Gabriel Renaud gabriel.reno@gmail.com
 *
 */

#include "ComputeDivergence.h"

// #define PRINTPROGRESS
// #define DEBUG
// #define DEBUG2
// #define DEBUG3
//#define DEBUGCOORD

static unsigned int limitToReOpenFP = 200; //if the coordinate is this far away, we will seek the file pointer




int computeDivergence(string refereVCF,string refereVCFidx,string sampleVCF,string sampleVCFidx, bool   useBedFileRegions,vector<GenomicRange> *  bedRegionsToFilter,string epoFile,string epoFileidx,int minimumNumberOfGoodSites,GenomicRange grc,SetVCFFilters * filtersVCFREF,SetVCFFilters * filtersVCFSMP,int bpForIndels,bool cpgInData,bool useHumanGenomeAsReference,int minPLdiffind,bool maximizeDiv){
    return computeDivergence( refereVCF, refereVCFidx, sampleVCF, sampleVCFidx, useBedFileRegions, bedRegionsToFilter, epoFile, epoFileidx, minimumNumberOfGoodSites, grc.getChrName(), grc.getStartCoord(), grc.getEndCoord(),  filtersVCFREF,filtersVCFSMP,bpForIndels,cpgInData,useHumanGenomeAsReference,minPLdiffind,maximizeDiv);
}


int computeDivergence( string refereVCF,string refereVCFidx,string sampleVCF,string sampleVCFidx,bool   useBedFileRegions,vector<GenomicRange> *  bedRegionsToFilter,string epoFile,string epoFileidx,int minimumNumberOfGoodSites,string chrName, int startChrCoord,  int endChrCoord,SetVCFFilters * filtersVCFREF,SetVCFFilters * filtersVCFSMP,int bpForIndels,bool cpgInData,bool useHumanGenomeAsReference,int minPLdiffind,bool maximizeDiv){
	
    //int minGQcutoff,int minCovcutoffREF,int maxCovcutoffREF,int minCovcutoffSMP,int maxCovcutoffSMP,int minMQcutoff,double minMapabilitycutoff,
    //int bpForIndels,bool cpgInData,bool useHumanGenomeAsReference,bool specifiedPL,int minPLdiffind){

    int numberOfGoodSites=0;

#ifdef PRINTPROGRESS
    vector<int>  milestonesCoord  (10,0);
    //vector<bool> milestonesPrint  (10,false);
    int difference= endChrCoord-startChrCoord;
    int currentMilestone=0;
    for(int i=0;i<10;i++){
	// cout<<i<<endl;
	milestonesCoord[i] = startChrCoord + (i+1)*int(double(difference)/double(10));
	// cout<<milestonesCoord[i]<<endl;
    }
    //    exit(1);
#endif		    
    unsigned int bedRegionsToFilterIndex=0;


    if(useBedFileRegions)
	if(bedRegionsToFilter->empty())     //if the vector is empty, might as well set the flag to false
	    useBedFileRegions=false;
    
    bool refIfVCF=false;

    //VCFreader * vcfREF;
    AlleleInfoReader * vcfREF;

    if(!useHumanGenomeAsReference){
	//READING VCF DATA
	if( strEndsWith(refereVCF,"vcf.gz")){
	    // cerr<<"ERROR: the reference must be a tabix indexed vcf"<<endl;
	    // exit(1);	
	    vcfREF=new VCFreader ( refereVCF.c_str(), refereVCFidx.c_str(), chrName, startChrCoord,endChrCoord,bpForIndels );    
	    refIfVCF=true;
	}else{
	    if( strEndsWith(refereVCF,"bed.gz")){
		// if(specifiedPL){
		cerr<<"WARNING: pl filtering is NOT used when using BAMTABLE"<<endl;
		//     //exit(1);
		// }
		vcfREF=new BAMTABLEreader ( refereVCF.c_str(), refereVCFidx.c_str(), chrName, startChrCoord,endChrCoord );    
	    }else{
		vcfREF=0;
		cerr<<"ERROR: the reference must be a tabix indexed vcf or bamtable file"<<endl;
		exit(1);
	    }
	}
    }else{
	vcfREF=0;
    }

    AlleleInfoReader * vcfSMP;
    bool smpIfVCF=false;

    if( strEndsWith(sampleVCF,"vcf.gz")){
	vcfSMP = new VCFreader ( sampleVCF.c_str(), sampleVCFidx.c_str(), chrName, startChrCoord,endChrCoord,bpForIndels );
	smpIfVCF=true;
    }else{
	if( strEndsWith(sampleVCF,"bed.gz")){
	    //if(specifiedPL){
	    //cerr<<"ERROR: cannot specify pl filtering when using BAMTABLE"<<endl;
	    cerr<<"WARNING: pl filtering is NOT used when using BAMTABLE"<<endl;
		//	exit(1);
		//}
	    vcfSMP = new BAMTABLEreader (sampleVCF.c_str(), sampleVCFidx.c_str(), chrName, startChrCoord,endChrCoord);
	}else{
	    cerr<<"ERROR: the sample must be a tabix indexed vcf or bamtable file"<<endl;
	    exit(1);
	}
    }
    // VCFreader vcfSMP ( sampleVCF.c_str(), sampleVCFidx.c_str(), chrName, startChrCoord,endChrCoord,bpForIndels );

    //SimpleVCF * smvcfREF;
    // SimpleVCF smvcfSMP;
    //AlleleInfo * smvcfREF;
    AlleleInfo * smvcfREF;
    AlleleInfo * smvcfSMP;

    bool lineLeftREF;
    if(useHumanGenomeAsReference){
	lineLeftREF=true;
    }else{
	lineLeftREF=vcfREF->hasData();
    }
    bool lineLeftSMP=vcfSMP->hasData();
    bool lineLeftEPO;

    if(!useHumanGenomeAsReference){
	if(lineLeftREF){ 
	    smvcfREF=vcfREF->getData();  
	}else{  
	    smvcfREF=0; 
	}
    }else{
	smvcfREF=0; 
    }

    if(lineLeftSMP){     
	smvcfSMP=vcfSMP->getData();  
    }else{  
	smvcfSMP=0; 
    }


    //READING EPO DATA CHECK IF FLAG FAIL

    ReadTabix rtEPO ( epoFile.c_str()  , epoFileidx.c_str()  , chrName, startChrCoord, endChrCoord ); //the destructor should be called automatically
    string lineFromEPO;
    lineLeftEPO=(rtEPO.readLine( lineFromEPO ));


	 

	    
    //Quit if there is no data to start with
    if(!( lineLeftREF &&  //reference
	  lineLeftSMP && //sample
	  lineLeftEPO )){ //EPO
	if(!useHumanGenomeAsReference)
	    delete vcfREF;
	delete vcfSMP;		

	return -1;
    }

    unsigned int coordCurrent=startChrCoord;
    //unsigned int previousCoord=0;

    DivergenceResult divr;
    // //This is everything
    // AlleleCounter all;
    // //This is without the ones marked as CpG
    // AlleleCounter noCpg;
    // //This is only with the ones marked as CpG
    // AlleleCounter onlyCpg;
    // //This excludes the following cases:
    // // S = C, R or A = T
    // // S = T, R or A = C
    // // S = A, R or A = G
    // // S = G, R or A = A    
    // AlleleCounter noTransitions;
    // //This excludes the following cases:
    // // S = T, R or A = C
    // // S = A, R or A = G    
    // AlleleCounter noDamage;

    //DEBUG 
 
    int rejectEPOValidREF=0;
    int rejectEPOValidANC=0;
    int rejectERROR_REF=0;

	    

    bool stayLoop=true;
    while(stayLoop){
	//reached the end
	if( !lineLeftREF ||  //reference
	    !lineLeftSMP ||  //sample
	    !lineLeftEPO ){ //EPO
	    stayLoop=false;
	    break;
	}

	vector<string> fieldsEPO=allTokens(lineFromEPO,'\t');

	//current position in the VCF files
	unsigned int positionREF; 
	unsigned int positionSMP=smvcfSMP->getPosition();
	unsigned int positionEPO=string2uint(fieldsEPO[1]);

	if(useHumanGenomeAsReference){
	    positionREF =positionEPO;
	}else{
	    positionREF =smvcfREF->getPosition();
	}

#ifdef PRINTPROGRESS
	if(coordCurrent > milestonesCoord[currentMilestone]){
	    cerr<<"Progress "<<(10*(currentMilestone+1))<<"%"<<endl;
	    //	    milestonesPrint[currentMilestone]=true;
	    currentMilestone++;
	}
#endif		

#ifdef DEBUG2
	cout<<endl<<"######"<<endl;
	cout<<"coordCurrent\t"<<coordCurrent<<endl;
	cout<<"positionREF\t"<<positionREF<<endl;
	cout<<"positionSMP\t"<<positionSMP<<endl;
	cout<<"positionEPO\t"<<positionEPO<<endl;
	if(useBedFileRegions){
	    cout<<"current bed     "<<bedRegionsToFilter->size()<<endl;
	    cout<<"current bed idx "<<bedRegionsToFilterIndex<<endl;
	}
#endif		

	if(useBedFileRegions){//we seek only coordinates in the regions given

	    //if coordinate is before next allowed range, skip to that coordinate
	    if(coordCurrent<(*bedRegionsToFilter)[bedRegionsToFilterIndex].getStartCoord() &&
	       coordCurrent<(*bedRegionsToFilter)[bedRegionsToFilterIndex].getEndCoord()  ){//
		coordCurrent=(*bedRegionsToFilter)[bedRegionsToFilterIndex].getStartCoord();
		continue;
	    }

	    //if coordinate is before next allowed range, skip to that coordinate
	    if( (*bedRegionsToFilter)[bedRegionsToFilterIndex].getStartCoord() <= coordCurrent &&
		                                                                  coordCurrent <= (*bedRegionsToFilter)[bedRegionsToFilterIndex].getEndCoord() ){
		//fine
		
	    }

	    //if coordinate is ahead, need to increase counter
	    if((*bedRegionsToFilter)[bedRegionsToFilterIndex].getStartCoord() < coordCurrent &&
	       (*bedRegionsToFilter)[bedRegionsToFilterIndex].getEndCoord()   < coordCurrent){
		coordCurrent=(*bedRegionsToFilter)[bedRegionsToFilterIndex].getStartCoord();
		if(bedRegionsToFilterIndex == (bedRegionsToFilter->size()-1))//last index
		    break;
		else{
		   bedRegionsToFilterIndex++;
		   continue;
		}
	    }
	    
	}


	//all at the same genomic coord
	if(coordCurrent == positionREF &&
	   coordCurrent == positionSMP &&
	   coordCurrent == positionEPO ){ 
#ifdef DEBUG2
	    cout<<"entered"<<endl;
	    if(!useHumanGenomeAsReference)
		cout<<*smvcfREF<<endl;
	    cout<<*smvcfSMP<<endl;
	    cout<<lineFromEPO<<endl;
#endif
	    //increase for next iteration
	    coordCurrent++;
	    


	    //continue;


	    //////////////////////////////////////
	    //                                  //
	    //   BEGIN  VCF ACTUAL PROCESSING   //
	    //                                  //
	    //////////////////////////////////////
	    //Checking EPO nucleotides
	    if(!validOneBP(fieldsEPO[2])){ rejectEPOValidREF++; continue; } //REF in the EPO   
	    if(!validOneBP(fieldsEPO[4])){ rejectEPOValidANC++; continue; } //chimp/human ancestor in EPO   

			
	    //DETERMINING CPG
	    bool isCpG=false;
	    //First check CpG as flagged in the ancestral states
	    //then check in reference or sample (if cpgInData is used)
	    if(fieldsEPO[9]     == "0"){
		if(cpgInData){
		    if(useHumanGenomeAsReference){
			if(smvcfSMP->isCpg() ){
			    isCpG=true;
			}

		    }else{
			if(smvcfREF->isCpg() || smvcfSMP->isCpg() ){
			    isCpG=true;
			}
		    }
		}
	    }else{
		if(fieldsEPO[9] == "1"){ 
		    isCpG=true;
		}else{
		    cerr<<"ERROR: wrong CpG flag in EPO"<<endl;
		    exit(1);
		}
	    }



	    /////////////////////////////////////////////////////////////////
	    //           FILTERING INDELS AND UNDEFINED REF ALLELE         //
	    /////////////////////////////////////////////////////////////////

	    bool refPassingFilter;
	    if(useHumanGenomeAsReference){
		refPassingFilter=true;
	    }else{
		if(refIfVCF){ //we only call passedFilters if the referencewas a VCF
		    //refPassingFilter = passedFilters( ((SimpleVCF *)smvcfREF),  minCovcutoffREF, maxCovcutoffREF, minMapabilitycutoff, minMQcutoff, minGQcutoff);
		    refPassingFilter = passedFilters( ((SimpleVCF *)smvcfREF),  filtersVCFREF);//minCovcutoffREF, maxCovcutoffREF, minMapabilitycutoff, minMQcutoff, minGQcutoff);
		}else{
		    refPassingFilter=true;
		}		   
	    }

	    bool smpPassingFilter;
	    if(smpIfVCF){ //we only call passedFilters if the sample was a VCF
		smpPassingFilter     = passedFilters( ((SimpleVCF *)smvcfSMP),  filtersVCFSMP);//minCovcutoffSMP, maxCovcutoffSMP, minMapabilitycutoff, minMQcutoff, minGQcutoff);
	    }else{
		smpPassingFilter=true;
	    }

	    if( !refPassingFilter ||  !smpPassingFilter){
		continue;
	    }





	    if(useHumanGenomeAsReference){

		//if this is not the case, this is not good since we have no way of ascertaining the reference allele
		//ideally, we should have the reference allele in BAMTable output
		if(smpIfVCF){ 

		    if( ( (SimpleVCF *)smvcfSMP )->getRef()  != fieldsEPO[2]) {
			cerr<<"ERROR: Ref not equal to EPO"<<endl;
			//cerr<<*smvcfREF<<endl;
			cerr<<*smvcfSMP<<endl;
			cerr<<lineFromEPO<<endl;
			rejectERROR_REF++;
			continue;
		    }
		}
		
	    }else{
		if(smpIfVCF && refIfVCF){ //both are VCF

		    if( (( (SimpleVCF *)smvcfREF )->getRef() != ( (SimpleVCF *)smvcfSMP )->getRef()) ){
			cerr<<"ERROR: Reference ref allele not equal to sample's"<<endl;
			cerr<<*smvcfREF<<endl;
			cerr<<*smvcfSMP<<endl;
			cerr<<lineFromEPO<<endl;
			rejectERROR_REF++;
			continue;
		    }
		}

		if(refIfVCF) //only the ref is VCF
		    if((( (SimpleVCF *)smvcfREF )->getRef()  != fieldsEPO[2]) ){
			cerr<<"ERROR: Ref not equal to EPO"<<endl;
			cerr<<*smvcfREF<<endl;
			cerr<<*smvcfSMP<<endl;
			cerr<<lineFromEPO<<endl;
			rejectERROR_REF++;
			continue;
		    }
	    }

	    /////////////////////////////////
	    //  PASSED QUALITY FILTERS     //
	    /////////////////////////////////

	    //previousCoord=positionREF; //they should all be equal

	    char allel_chimpHumanAncestor = fieldsEPO[3][0];
	    char allel_reference          ;

	    if(useHumanGenomeAsReference){
		allel_reference          = fieldsEPO[2][0]; //we use the reference allele from the EPO (NOT the one from a VCF) to avoid biasing towards unresolved regions in the EPO
	    }else{
		//if(specifiedPL){
		if(refIfVCF){
		    allel_reference          = ( (SimpleVCF *)smvcfREF )->getRandomAlleleUsingPL(minPLdiffind);
		    if(allel_reference == 'X'){ continue; } //PL values not high enough to generate an allele
		}else{
		    allel_reference          = smvcfREF->getRandomAllele();
		}		
	    }

	    
	    char allel_sample;
	    // if(specifiedPL){
	    if(smpIfVCF){
		allel_sample= ( (SimpleVCF *)smvcfSMP )->getRandomAlleleUsingPL(minPLdiffind);
		if(allel_sample == 'X'){ continue; } //PL values not high enough to generate an allele
	    }else{
		allel_sample= smvcfSMP->getRandomAllele();
	    }

	    if(maximizeDiv){
		if(useHumanGenomeAsReference){
		    //check for het in smp
		    // this double check is done because sometimes, the most likely genotype is het but no alt is produced
		    if( ((SimpleVCF *)smvcfSMP )->isHeterozygous() && 
			((SimpleVCF *)smvcfSMP )->isHeterozygousUsingPL(minPLdiffind) ) {
			if(allel_reference == allel_sample){
			    allel_sample = ((SimpleVCF *)smvcfSMP )->getOtherAllele(allel_sample);
			}
		    }

		}else{
		    //check for both het
		    if( ((SimpleVCF *)smvcfREF )->isHeterozygous()                    && 
			((SimpleVCF *)smvcfREF )->isHeterozygousUsingPL(minPLdiffind) &&		      
			((SimpleVCF *)smvcfSMP )->isHeterozygous()                    && 
			((SimpleVCF *)smvcfSMP )->isHeterozygousUsingPL(minPLdiffind) ) {
			if(allel_reference == allel_sample){
			    if(randomBool())
				allel_reference = ((SimpleVCF *)smvcfREF )->getOtherAllele(allel_reference);
			    else
				allel_sample    = ((SimpleVCF *)smvcfSMP )->getOtherAllele(allel_sample);
			}
		    }else{
			if( ((SimpleVCF *)smvcfREF )->isHeterozygous()                    && 
			    ((SimpleVCF *)smvcfREF )->isHeterozygousUsingPL(minPLdiffind) ){
			    if(allel_reference == allel_sample)				
				allel_reference = ((SimpleVCF *)smvcfREF )->getOtherAllele(allel_reference);
			}else{
			    
			    if( ((SimpleVCF *)smvcfSMP )->isHeterozygous()                    && 
			        ((SimpleVCF *)smvcfSMP )->isHeterozygousUsingPL(minPLdiffind) ) {
				if(allel_reference == allel_sample){
				    allel_sample    = ((SimpleVCF *)smvcfSMP )->getOtherAllele(allel_sample);
				}
			    }

			}

		    }


		}

	    }

	    numberOfGoodSites++; //increase the counter for # of sites passed quality filters
#ifdef  DEBUGCOORD
	    cout<<"coordCurrent "<<(coordCurrent-1)<<endl;
#endif
	    computeDiv(allel_chimpHumanAncestor,allel_reference,allel_sample,isCpG,&divr);

			     

#ifdef DEBUG3
	    if(!useHumanGenomeAsReference){	
		cout<<"lineFromREF "<<*smvcfREF<<endl;
	    }
	    cout<<"lineFromSMP "<<*smvcfSMP<<endl;
	    cout<<"lineFromEPO "<<lineFromEPO<<endl;
	    cout<<all<<endl<<endl;
#endif
	    

	    //////////////////////////////////////
	    //                                  //
	    //   END  VCF ACTUAL PROCESSING     //
	    //                                  //
	    //////////////////////////////////////





			  
	}else{
#ifdef DEBUG2

	    cout<<"needs line ajusting "<<endl;
#endif


	    //REPOSITING THE FILES AND COUNTER
	    if(!useHumanGenomeAsReference){

		if(coordCurrent < positionREF){ //overshot , repositioning there
		    coordCurrent = positionREF;
		    continue;
		}

		if(coordCurrent == positionREF){ //fine
		}

		if(coordCurrent > positionREF){ //running behind

		    if( (coordCurrent - positionREF ) >= limitToReOpenFP){ //seeking in the file
			vcfREF->repositionIterator( chrName, int(coordCurrent), endChrCoord);
		    }

		    lineLeftREF=vcfREF->hasData();
		    if(lineLeftREF){
			smvcfREF=vcfREF->getData();
		    }
		}
	    }

	    if(coordCurrent < positionSMP){ //overshot , repositioning there
		coordCurrent = positionSMP;
		continue;
	    }

	    if(coordCurrent == positionSMP){ //fine
	    }

	    if(coordCurrent > positionSMP){ //running behind
		if( (coordCurrent - positionSMP ) >= limitToReOpenFP){ //seeking in the file
		    vcfSMP->repositionIterator( chrName, int(coordCurrent),endChrCoord);
		}

		lineLeftSMP=vcfSMP->hasData();
		if(lineLeftSMP){
		    smvcfSMP=vcfSMP->getData();
		}

	    }
			    

	    if(coordCurrent < positionEPO){ //overshot , repositioning there
		coordCurrent = positionEPO;
		continue;
	    }

	    if(coordCurrent == positionEPO){ //fine
	    }

	    if(coordCurrent > positionEPO){ //running behind
		if( (coordCurrent - positionEPO ) >= limitToReOpenFP){ //seeking in the file
		    rtEPO.repositionIterator(chrName, int(coordCurrent),endChrCoord ); 		    
		}
		
		lineLeftEPO=(rtEPO.readLine( lineFromEPO ));
	    }

	}
			   
			
		
    }//end of data processing for a given random chromosome (closing while(stayInLoop))

#ifdef DEBUG


    cout<<"rejectEPOValidREF     "<<rejectEPOValidREF<<endl;
    cout<<"rejectEPOValidANC     "<<rejectEPOValidANC<<endl;
    cout<<"rejectERROR_REF       "<<rejectERROR_REF<<endl;

	    
#endif
	
    //destructor for EPO reader is implicitely called (allocation on the stack)   
    if(!useHumanGenomeAsReference)
	delete vcfREF;
    delete vcfSMP;
			
    //we found a satisfactory # of sites
    if(numberOfGoodSites >= minimumNumberOfGoodSites){
	cout<<chrName<<":"
	    <<startChrCoord<<"-"
	    <<endChrCoord<<"\t"
	    <<divr<<"\t"
	    <<endl;
	return 0;
    }else{
	cout<<chrName<<":"
	    <<startChrCoord<<"-"
	    <<endChrCoord<<"\tN/A"<<endl;
	cerr<<"Failed enough sites good sites for "<<chrName<<":"<<startChrCoord<<"-"<<endChrCoord<<"\t"<<endl;
	cerr<<"Needed "<<minimumNumberOfGoodSites<<" found "<<numberOfGoodSites<<endl;
	return 1;
    }

    
	    
}//end of if found random coordinate in genome


