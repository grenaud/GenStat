/*
 * Dstats
 * Date: Oct-15-2012 
 * Author : Gabriel Renaud gabriel.reno@gmail.com
 *
 */

#include "Dstats.h"
//#include "ComputeDivergence_core.h"

// #define DEBUG
// #define DEBUG2
// #define DEBUG3

static unsigned int limitToReOpenFP = 200; //if the coordinate is this far away, we will re-open the file pointer


// inline bool isTransition(char allel_derived,char allel_ancestral){
//     return (  
// 	    (allel_ancestral == 'C' && allel_derived == 'T' )  || 
// 	    (allel_ancestral == 'T' && allel_derived == 'C' )  || 
// 	    (allel_ancestral == 'A' && allel_derived == 'G' )  || 
// 	    (allel_ancestral == 'G' && allel_derived == 'A' )  
// 	      );
// }


int dstats(string refereVCF,string refereVCFidx,string sampleVCF,string sampleVCFidx,string conditionalVCF,string conditionalVCFidx,string epoFile,string epoFileidx,int minimumNumberOfGoodSites,GenomicRange grc,
	   SetVCFFilters * filtersVCFREF,
	   SetVCFFilters * filtersVCFSMP,
	   SetVCFFilters * filtersVCFCOND,
	   int bpForIndels,bool cpgInData,bool useHumanGenomeAsReference,bool useEPOInferedAsAncestor,bool requiredHomozygousCond){
	   //int minQCcutoff,int minCovcutoffREF,int maxCovcutoffREF,int minCovcutoffSMP,int maxCovcutoffSMP,int minCovcutoffCOND,int maxCovcutoffCOND,int minMQcutoff,double minMapabilitycutoff,int bpForIndels,bool cpgInData,bool useHumanGenomeAsReference,bool useEPOInferedAsAncestor,bool requiredHomozygousCond){
    return dstats( refereVCF, refereVCFidx, sampleVCF, sampleVCFidx, conditionalVCF, conditionalVCFidx,epoFile, epoFileidx, minimumNumberOfGoodSites, grc.getChrName(), grc.getStartCoord(), grc.getEndCoord(), filtersVCFREF, filtersVCFSMP,filtersVCFCOND, bpForIndels,cpgInData,useHumanGenomeAsReference,useEPOInferedAsAncestor,requiredHomozygousCond);
		   
		   //minQCcutoff, minCovcutoffREF, maxCovcutoffREF, minCovcutoffSMP, maxCovcutoffSMP, minCovcutoffCOND,maxCovcutoffCOND, minMQcutoff, minMapabilitycutoff, bpForIndels,cpgInData,useHumanGenomeAsReference,useEPOInferedAsAncestor,requiredHomozygousCond);
}


int dstats( string refereVCF,string refereVCFidx,string sampleVCF,string sampleVCFidx,string conditionalVCF,string conditionalVCFidx,string epoFile,string epoFileidx,int minimumNumberOfGoodSites,string chrName, int startChrCoord,  int endChrCoord,SetVCFFilters * filtersVCFREF,
	   SetVCFFilters * filtersVCFSMP,
	    SetVCFFilters * filtersVCFCOND,int bpForIndels,bool cpgInData,bool useHumanGenomeAsReference,bool useEPOInferedAsAncestor,bool requiredHomozygousCond){
	    //int minQCcutoff,int minCovcutoffREF,int maxCovcutoffREF,int minCovcutoffSMP,int maxCovcutoffSMP,int minCovcutoffCOND,int maxCovcutoffCOND,int minMQcutoff,double minMapabilitycutoff,int bpForIndels,bool cpgInData,bool useHumanGenomeAsReference,bool useEPOInferedAsAncestor,bool requiredHomozygousCond){

     int numberOfGoodSites=0;


     VCFreader * vcfREF;

     if(!useHumanGenomeAsReference){
	 //READING VCF DATA
	 if( !strEndsWith(refereVCF,"vcf.gz")){
	     cerr<<"ERROR: the reference "<< refereVCF <<" must be a tabix indexed vcf "<<endl;
	     exit(1);	
	 }

	 vcfREF=new VCFreader ( refereVCF.c_str(), refereVCFidx.c_str(), chrName, startChrCoord,endChrCoord,bpForIndels );    
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
	     vcfSMP = new BAMTABLEreader (sampleVCF.c_str(), sampleVCFidx.c_str(), chrName, startChrCoord,endChrCoord);
	 }else{
	     cerr<<"ERROR: the sample must be a tabix indexed vcf or bamtable file"<<endl;
	     exit(1);
	 }
     }

     // VCFreader * vcfSMP;   
     // //READING VCF DATA
     // if( !strEndsWith(sampleVCF,"vcf.gz")){
     // 	cerr<<"ERROR: the reference must be a tabix indexed vcf"<<endl;
     // 	exit(1);	
     // }
     // vcfSMP=new VCFreader ( sampleVCF.c_str(), sampleVCFidx.c_str(), chrName, startChrCoord,endChrCoord,bpForIndels );    


     VCFreader * vcfCOND;

     //READING VCF DATA
     if( !strEndsWith(conditionalVCF,"vcf.gz")){
	 cerr<<"ERROR: the conditional must be a tabix indexed vcf"<<endl;
	 exit(1);	
     }

     vcfCOND=new VCFreader ( conditionalVCF.c_str(), conditionalVCFidx.c_str(), chrName, startChrCoord,endChrCoord,bpForIndels );    








     SimpleVCF * smvcfREF;
     //    SimpleVCF * smvcfSMP;
     AlleleInfo * smvcfSMP;

     SimpleVCF * smvcfCOND;

     bool lineLeftREF;
     if(useHumanGenomeAsReference){
	 lineLeftREF=true;
     }else{
	 lineLeftREF=vcfREF->hasData();
     }
     bool lineLeftSMP =vcfSMP->hasData();
     bool lineLeftCOND=vcfCOND->hasData();

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
	 smvcfSMP  = vcfSMP->getData();  
     }else{
	 smvcfSMP  = 0;
     }
     if(lineLeftCOND){ 
	 smvcfCOND = vcfCOND->getData();  
     }else{
	 smvcfCOND = 0;
     }



     //READING EPO DATA CHECK IF FLAG FAIL
     ReadTabix rtEPO ( epoFile.c_str()  , epoFileidx.c_str()  , chrName, startChrCoord, endChrCoord ); //the destructor should be called automatically
     string lineFromEPO;
     lineLeftEPO=(rtEPO.readLine( lineFromEPO ));






     //Quit if there is no data to start with
     if(!( lineLeftREF &&  //reference
	   lineLeftSMP && //sample
	   lineLeftCOND && //sample
	   lineLeftEPO )){ //EPO
	 if(!useHumanGenomeAsReference)
	     delete vcfREF;
	 delete vcfSMP;		
	 delete vcfCOND;		

	 return -1;
     }

     unsigned int coordCurrent=startChrCoord;
     unsigned int previousCoord=0;

     // //This is everything
     // DstatCounter all;
     // //This is without the ones marked as CpG
     // DstatCounter noCpg;
     // //This is only with the ones marked as CpG
     // DstatCounter onlyCpg;
     // //This is only with the ones marked as CpG
     // DstatCounter noTransitions;
     DstatResult dSres;

     //DEBUG 

     int rejectEPOValidREF=0;
     int rejectEPOValidANC=0;
     int rejectERROR_REF=0;



     bool stayLoop=true;
     while(stayLoop){
	 //reached the end
	 if( !lineLeftREF  ||  //reference
	     !lineLeftSMP  ||  //sample
	     !lineLeftCOND ||  //sample
	     !lineLeftEPO ){ //EPO
	     stayLoop=false;
	     break;
	 }

	 vector<string> fieldsEPO=allTokens(lineFromEPO,'\t');

	 //current position in the VCF files
	 unsigned int positionREF; 
	 unsigned int positionSMP  = smvcfSMP->getPosition();
	 unsigned int positionCOND = smvcfCOND->getPosition();

	 unsigned int positionEPO=string2uint(fieldsEPO[1]);

	 if(useHumanGenomeAsReference){
	     positionREF =positionEPO;
	 }else{
	     positionREF =smvcfREF->getPosition();
	 }


 #ifdef DEBUG2
	 cout<<endl<<"######"<<endl;
	 cout<<"coordCurrent\t"<<coordCurrent<<endl;
	 cout<<"positionREF\t"<<positionREF<<endl;
	 cout<<"positionSMP\t"<<positionSMP<<endl;
	 cout<<"positionCOND\t"<<positionCOND<<endl;

	 cout<<"positionEPO\t"<<positionEPO<<endl;
 #endif		


	 //all at the same genomic coord
	 if(coordCurrent == positionREF  &&
	    coordCurrent == positionSMP  &&
	    coordCurrent == positionCOND &&
	    coordCurrent == positionEPO ){ 
 #ifdef DEBUG2
	     cout<<"entered"<<endl;
	     if(!useHumanGenomeAsReference)
		 cout<<"REF  "<<*smvcfREF<<endl;
	     cout<<"SMP  "<<*smvcfSMP<<endl;
	     cout<<"COND "<<*smvcfCOND<<endl;
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
		     if(smvcfREF->isCpg() || smvcfSMP->isCpg() ){
			 isCpG=true;
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
		 refPassingFilter=passedFilters(               smvcfREF  , filtersVCFREF);//minCovcutoffREF, maxCovcutoffREF, minMapabilitycutoff,minMQcutoff,minQCcutoff);
	     }

	     bool smpPassingFilter;
	     if(smpIfVCF){ //we only call passedFilters if the sample was a VCF
		 smpPassingFilter=passedFilters( ((SimpleVCF *)smvcfSMP)  , filtersVCFSMP);//minCovcutoffSMP, maxCovcutoffSMP, minMapabilitycutoff,minMQcutoff,minQCcutoff);
	     }else{
		 smpPassingFilter=true;
	     }

	     bool condPassingFilter;	    
	     condPassingFilter=passedFilters(                  smvcfCOND  , filtersVCFCOND);//minCovcutoffCOND,maxCovcutoffCOND,minMapabilitycutoff,minMQcutoff,minQCcutoff);


	     if( !refPassingFilter ||  !smpPassingFilter || !condPassingFilter){
		 continue;
	     }







	     if(useHumanGenomeAsReference){

		 if(smpIfVCF){ 
		     if(  ((SimpleVCF *)smvcfSMP )->getRef()  != fieldsEPO[2]) {
			 cerr<<"ERROR: Ref not equal to EPO"<<endl;
			 //cerr<<*smvcfREF<<endl;
			 cerr<<*smvcfSMP<<endl;
			 cerr<<lineFromEPO<<endl;
			 rejectERROR_REF++;
			 continue;
		     }
		 }


		 if(  smvcfCOND->getRef()  != fieldsEPO[2]) {
		     cerr<<"ERROR: Ref not equal to EPO"<<endl;
		     //cerr<<*smvcfREF<<endl;
		     cerr<<*smvcfSMP<<endl;
		     cerr<<lineFromEPO<<endl;
		     rejectERROR_REF++;
		     continue;
		 }




	     }else{

		 if(smpIfVCF){ 
		     if( smvcfREF->getRef() !=  ((SimpleVCF *)smvcfSMP )->getRef() ){
			 cerr<<"ERROR: Reference ref allele not equal to sample's"<<endl;
			 cerr<<*smvcfREF<<endl;
			 cerr<<*smvcfSMP<<endl;
			 cerr<<*smvcfCOND<<endl;
			 cerr<<lineFromEPO<<endl;
			 rejectERROR_REF++;
			 continue;
		     }
		 }

		 if( smvcfREF->getRef() !=  smvcfCOND->getRef() ){
		     cerr<<"ERROR: Reference ref allele not equal to sample's"<<endl;
		     cerr<<*smvcfREF<<endl;
		     cerr<<*smvcfSMP<<endl;
		     cerr<<*smvcfCOND<<endl;
		     cerr<<lineFromEPO<<endl;
		     rejectERROR_REF++;
		     continue;
		 }


		 if((smvcfREF->getRef()  != fieldsEPO[2]) ){
		     cerr<<"ERROR: Ref not equal to EPO"<<endl;
		     cerr<<*smvcfREF<<endl;
		     cerr<<*smvcfSMP<<endl;
		     cerr<<*smvcfCOND<<endl;
		     cerr<<lineFromEPO<<endl;
		     rejectERROR_REF++;
		     continue;
		 }


	     }





	     /////////////////////////////////
	     //  PASSED QUALITY FILTERS     //
	     /////////////////////////////////
	     numberOfGoodSites++; //increase the counter for # of sites passed quality filters

	     previousCoord=positionREF; //they should all be equal

	     char allel_chimpHumanAncestor;
	     //if we use the infered ancestral allele as ancestral allele
	     if(useEPOInferedAsAncestor){
		 allel_chimpHumanAncestor  = fieldsEPO[3][0];
	     }else{
		 allel_chimpHumanAncestor  = fieldsEPO[4][0];
	     }
	     char allel_condition          ;


	     //BEGIN USING THE CONDITIONAL FOR FILTERING FOR DERIVED SITES
	     if(requiredHomozygousCond){
		 if(smvcfCOND->isHomozygousREF() || smvcfCOND->isHomozygousALT() ){ //we require homozygosity

		     //if the conditional is homozygous reference and 
		     //derived (not ancestral)
		     //then the reference becomes the derived allele
		     if( smvcfCOND->isHomozygousREF() && 
			 (smvcfCOND->getRef()[0] !=  allel_chimpHumanAncestor) ){
			 allel_condition = smvcfCOND->getRef()[0];
		     }else{

			 //if the conditional is homozygous alternative and 
			 //derived (not ancestral)
			 //then the alternative becomes the derived allele
			 if( smvcfCOND->isHomozygousALT() && 
			     (smvcfCOND->getAlt()[0] !=  allel_chimpHumanAncestor) ){
			     allel_condition = smvcfCOND->getAlt()[0];

			 }else{ //neither condition was satisfied, skipping to next site
			     continue;
			 }

		     }

		 }else{  //else, move to next site
		     continue;
		 }
	     }else{
		 allel_condition          = smvcfCOND->getRandomAllele();
	     }
	     //END USING THE CONDITIONAL FOR FILTERING



	     char allel_reference          ;

	     if(useHumanGenomeAsReference){
		 allel_reference          = fieldsEPO[2][0]; //we use the reference allele from the EPO (NOT the one from a VCF) to avoid biasing towards unresolved regions in the EPO
	     }else{
		 allel_reference          = smvcfREF->getRandomAllele();
	     }
	     char allel_sample             = smvcfSMP->getRandomAllele();


	     computeDstat(allel_chimpHumanAncestor,allel_condition,allel_reference,allel_sample,isCpG,&dSres);


	     // //require the reference allele to be either derived or ancestral
	     // if( ( allel_reference == allel_chimpHumanAncestor) ||
	     // 	 ( allel_reference == allel_condition)   ){
	     // 	 //fine
	     // }else{
	     // 	 continue;
	     // }

	     // //require the sample allele to be either derived or ancestral
	     // if( ( allel_sample == allel_chimpHumanAncestor) ||
	     // 	 ( allel_sample == allel_condition)   ){
	     // 	 //fine
	     // }else{
	     // 	 continue;
	     // }


	     // //AA
	     //  if( ( allel_reference == allel_chimpHumanAncestor) &&
	     // 	  ( allel_sample    == allel_chimpHumanAncestor) ){


	     // 	  all.counterAncAnc++;
	     // 	  if(isCpG)
	     // 	     onlyCpg.counterAncAnc++;
	     // 	 else
	     // 	     noCpg.counterAncAnc++;

	     // 	  if(!isPotentialTransition(allel_condition,allel_chimpHumanAncestor))
	     // 	      noTransitions.counterAncAnc++;

	     //  }else{
	     // 	  //AD
	     // 	  if( ( allel_reference == allel_chimpHumanAncestor) &&
	     // 	      ( allel_sample    == allel_condition) ){



	     // 	      all.counterAncDer++;
	     // 	      if(isCpG)
	     // 		  onlyCpg.counterAncDer++;
	     // 	      else
	     // 		  noCpg.counterAncDer++;

	     // 	      if(!isPotentialTransition(allel_condition,allel_chimpHumanAncestor))
	     // 		  noTransitions.counterAncDer++;


	     // 	  }else{
	     // 	      //DA
	     // 	      if( ( allel_reference == allel_condition) &&
	     // 		  ( allel_sample    == allel_chimpHumanAncestor) ){


	     // 		  all.counterDerAnc++;
	     // 		  if(isCpG)
	     // 		      onlyCpg.counterDerAnc++;
	     // 		  else
	     // 		      noCpg.counterDerAnc++;

	     // 		  if(!isPotentialTransition(allel_condition,allel_chimpHumanAncestor))
	     // 		      noTransitions.counterDerAnc++;



	     // 	      }else{

	     // 		  //DD
	     // 		  if( ( allel_reference == allel_condition) &&
	     // 		      ( allel_sample    == allel_condition) ){



	     // 		      all.counterDerDer++;
	     // 		      if(isCpG)
	     // 			  onlyCpg.counterDerDer++;
	     // 		      else
	     // 			  noCpg.counterDerDer++;

	     // 		      if(!isTransition(allel_condition,allel_chimpHumanAncestor))
	     // 			  noTransitions.counterDerDer++;


	     // 		  }else{
	     // 		      cerr<<"Invalid state for "<<endl;	     
	     // 		      cerr<<"lineFromREF "<<*smvcfREF<<endl;
	     // 		      cerr<<"lineFromSMP "<<*smvcfSMP<<endl;
	     // 		      cerr<<"lineFromCOND "<<*smvcfCOND<<endl;
	     // 		  }


	     // 	      }

	     // 	  }


	     //  }






 #ifdef DEBUG3
	     cout<<"lineFromREF "<<*smvcfREF<<endl;
	     cout<<"lineFromSMP "<<*smvcfSMP<<endl;
	     cout<<"lineFromCOND "<<*smvcfCOND<<endl;
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
			 //vecOfAlleleReaders[indexAlleleFile]->repositionIterator( grc.getChrName(), int(coordCurrent), int(grc.getEndCoord()) );
			 vcfREF->repositionIterator( chrName, int(coordCurrent),endChrCoord);
		     }else{
			 //nothing to do
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
		     //vecOfAlleleReaders[indexAlleleFile]->repositionIterator( grc.getChrName(), int(coordCurrent), int(grc.getEndCoord()) );
		     vcfSMP->repositionIterator( chrName, int(coordCurrent),endChrCoord);
		 }else{
		     //nothing to do
		 }


		 lineLeftSMP=vcfSMP->hasData();
		 if(lineLeftSMP){
		     smvcfSMP=vcfSMP->getData();
		 }

	     }




	     if(coordCurrent < positionCOND){ //overshot , repositioning there
		 coordCurrent = positionCOND;
		 continue;
	     }

	     if(coordCurrent == positionCOND){ //fine
	     }

	     if(coordCurrent > positionCOND){ //running behind

		 if( (coordCurrent - positionCOND ) >= limitToReOpenFP){ //seeking in the file
		     //vecOfAlleleReaders[indexAlleleFile]->repositionIterator( grc.getChrName(), int(coordCurrent), int(grc.getEndCoord()) );
		     vcfCOND->repositionIterator( chrName, int(coordCurrent),endChrCoord);
		 }else{
		     //nothing to do
		 }

		 lineLeftCOND=vcfCOND->hasData();
		 if(lineLeftCOND){
		     smvcfCOND=vcfCOND->getData();
		 }

	     }




	     //reposition EPO
	     if(coordCurrent < positionEPO){ //overshot , repositioning there
		 coordCurrent = positionEPO;
		 continue;
	     }

	     if(coordCurrent == positionEPO){ //fine
	     }

	     if(coordCurrent > positionEPO){ //running behind

		 if( (coordCurrent - positionEPO ) >= limitToReOpenFP){ //seeking in the file
		    rtEPO.repositionIterator( chrName, int(coordCurrent),endChrCoord);
		}
		//retrieving new data
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
     delete vcfCOND;

     //we found a satisfactory # of sites
     if(numberOfGoodSites >= minimumNumberOfGoodSites){
	 cout<<chrName<<":"
	     <<startChrCoord<<"-"
	     <<endChrCoord<<"\t"
	     <<dSres
	     // <<all<<"\t"
	     // <<onlyCpg<<"\t"
	     // <<noCpg<<"\t"
	     // <<noTransitions<<"\t"
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


