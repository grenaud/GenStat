/*
 * MultiVCFParser
 * Date: Oct-04-2012 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include "MultiVCFParser.h"

static unsigned int limitToReOpenFP = 200; //if the coordinate is this far away, we will seek the file pointer
		
// #define DEBUG


inline void MultiVCFParser::addPair( pair<int,int> * destination, pair<int,int>  toAdd){
    destination->first   +=toAdd.first;
    destination->second  +=toAdd.second;
}


inline bool MultiVCFParser::isPairEmpty( pair<int,int>  & toCheck){
    if(toCheck.first == 0 && toCheck.second == 0)
	return true;
    return false;
}



inline int MultiVCFParser::refBp2Index(const char toCheck){
    if(toCheck == 'A')
	return 1;
    if(toCheck == 'C')
	return 2;
    if(toCheck == 'G')
	return 3;
    if(toCheck == 'T')
	return 4;

    return -1; //invalid
    // cerr<<"Invalid ref bp = "<<toCheck<<endl;
    // exit(1);
}




inline int MultiVCFParser::altBp2Index(const char toCheck){
    if(toCheck == '.')
	return 0;
    if(toCheck == 'A')
	return 1;
    if(toCheck == 'C')
	return 2;
    if(toCheck == 'G')
	return 3;
    if(toCheck == 'T')
	return 4;

    return -1;
    // cerr<<"Invalid alt bp = "<<toCheck<<endl;
    // exit(1);
}

void MultiVCFParser::printCutoffs(){

    for(unsigned int indexCut=0;indexCut<vcfCutoffs->size();indexCut++){
	cerr<<(*(*vcfCutoffs)[indexCut])<<endl;
    }
}

MultiVCFParser::MultiVCFParser(string popfile,bool printChimp,int minPLdiffind,bool printChrCoord,int minIndWithAlt,int allowedFailPop,
			       int    minGQcutoff         ,
			       int    minMQcutoff         ,
			       double minMapabilitycutoff ,
			       bool   filterIndelProx     ,
			       bool   repeatMasking       ,
			       bool   systemError   ,
			       bool   allowall,
			       bool   allowallMQ){  
    

    //int minGQcutoff,int minMQcutoff,double minMapabilitycutoff){
    this->printChrCoord = printChrCoord;
    this->printChimp    = printChimp;
    this->minPLdiffind  = minPLdiffind;
    this->minIndWithAlt = minIndWithAlt;
    this->allowedFailPop= allowedFailPop;

    
    pop2PrefixFiles   = new map< string,vector<string> * >    ();
    pop2SuffixFiles   = new map< string,vector<string> * >    ();

    populations = new vector<string> ();
    //coverageCutoffs = new vector< pair<int,int> >; //vector of cutoffs for every vcf file
    vcfCutoffs  = new vector<SetVCFFilters *> ();


    ifstream myFile;	
    string line;

    if(!isFile(popfile)){
	cerr << "File "<<popfile<<" is not a file"<<endl;
	exit(1);
    }

    myFile.open(popfile.c_str(), ios::in);

    string popName="";
    popName.erase();
    if (myFile.is_open()){
	while ( getline (myFile,line)){
	    //line.erase(remove_if( line.begin(), line.end(), ::isspace), line.end());
	    trimWhiteSpacesBothEnds(&line);
	    if(line.length() == 0)
		continue;
	    if(line.substr(0,1) == "#")
		continue;

	    if(line.substr(line.length()-1,1) == ":"){ //new population
		popName=line.substr(0,line.length()-1);
		(*pop2PrefixFiles)[popName]   = new vector<string> ();	
		(*pop2SuffixFiles)[popName]   = new vector<string> ();	
		
		populations->push_back(popName);
		cerr<<"new population\t"<<popName<<endl;
	    }else{
				
		if(popName == ""){ //unset yet
		    cerr << "Population name for "<<line<<" was not defined"<<endl;
		    exit(1);
		}


		vector<string> tempLineSplit=allTokens(line,' ');
		if(tempLineSplit.size() == 2 ){ //no coverage bounds information
		    cerr<<"new file for "<<popName<< " "<<tempLineSplit[0] <<" "<<tempLineSplit[1] << " no coverage cutoffs"<<endl;
		    ((*pop2PrefixFiles)[popName])->push_back( tempLineSplit[0] );
		    ((*pop2SuffixFiles)[popName])->push_back( tempLineSplit[1] );
		    
		    SetVCFFilters * filterToAdd =new SetVCFFilters(minGQcutoff         ,
								   minMQcutoff         ,
								   minMapabilitycutoff ,
								   filterIndelProx     ,
								   repeatMasking       ,
								   systemError         ,
								   0                   ,
								   10000       ,
								   allowall,
								   allowallMQ);
		    filterToAdd->setName( "pop name "+popName+" file : "+tempLineSplit[0]+" "+tempLineSplit[1]);
		    vcfCutoffs->push_back( filterToAdd );
		    

		}else{
		    if(tempLineSplit.size() == 4 ){ //with coverage cutoffs information
			cerr<<"new file for "<<popName<< " "<<tempLineSplit[0] <<" "<<tempLineSplit[1] << " coverage cutoffs:"<<tempLineSplit[2]<<"-"<<tempLineSplit[3]<<endl;
			((*pop2PrefixFiles)[popName])->push_back( tempLineSplit[0] );
			((*pop2SuffixFiles)[popName])->push_back( tempLineSplit[1] );
			
			SetVCFFilters * filterToAdd = new SetVCFFilters(minGQcutoff         ,
									minMQcutoff         ,
									minMapabilitycutoff ,
									filterIndelProx     ,
									repeatMasking       ,
									systemError         ,
									destringify<int>(tempLineSplit[2]) , 
									destringify<int>(tempLineSplit[3])  ,
									allowall,
									allowallMQ);
			
			filterToAdd->setName( "pop name "+popName+" file : "+tempLineSplit[0]+" "+tempLineSplit[1] );
			vcfCutoffs->push_back( filterToAdd );
			//vcfCutoffs->back()->setName( "pop name "+popName+" file : "+tempLineSplit[0]+" "+tempLineSplit[1] );
		    }else{
			cerr << "Invalid number of fields for line "<< line<<endl<<"we allow 2 or 4 (with coverage cutoffs) only"<<endl;
			exit(1);
		    }
		}

	    }

	}
	myFile.close();
    }else{ 
	cerr << "Unable to open file "<<popfile<<endl;
	exit(1);
    }



    cerr<<"Done parsing pop file file"<<endl;

}

MultiVCFParser::~MultiVCFParser(){

    for(unsigned int i=0;i<populations->size();i++){
	delete( (*pop2PrefixFiles)[ populations->at(i) ] ) ;    
	delete( (*pop2SuffixFiles)[ populations->at(i) ] ) ;    
    }

    delete(pop2PrefixFiles);
    delete(pop2SuffixFiles);
    //delete(coverageCutoffs);
    //delete(vcfCutoffs);
    for(unsigned int indexCut=0;indexCut<vcfCutoffs->size();indexCut++){
	delete( (*vcfCutoffs)[indexCut] );
    }
    delete( vcfCutoffs );
    delete(populations);
}

int MultiVCFParser::produceOutput(string epoFile,string epoFileidx,GenomicRange grc,int outputType,bool onlySegSites,bool useEPOInferedAsAncestor,string programLine,string gitHubVersion){


    //EPO FILES

    ReadTabix * rtEPO = new ReadTabix( epoFile.c_str()  , epoFileidx.c_str()  , grc.getChrName(), grc.getStartCoord(), grc.getEndCoord() ); 
    string lineFromEPO="";
    bool lineLeftEPO=(rtEPO->readLine( lineFromEPO ));


    //  opening VCF or BAMtable files
    vector<AlleleInfoReader *> vecOfAlleleReaders;
    vector<string> filesNames;
    vector<string> filesNamestbi;
    vector<int>    filesNamesType;  // 1 = VCF, 2 = BAMTable
    vector<int>    popIndexForFile;
    vector<int>    popFileTypes (populations->size(),0);    // 1 = VCF, 2 = BAMTable cannot mix VCF and BAMTable in the same population
    string  fSeparator=" ";
    if(outputType == 1 )
	fSeparator=" ";
    if(outputType == 2 )
	fSeparator="\t";

    //opening the files
    //for each population
    for(unsigned int i=0;i<populations->size();i++){
	//for each file within each population
    	for(unsigned int j=0;j< ( (*pop2PrefixFiles)[ populations->at(i) ] )->size();j++){

    	    string fileToOpen    = ( (*pop2PrefixFiles)[ populations->at(i) ] )->at(j) + grc.getChrName() + ( (*pop2SuffixFiles)[ populations->at(i) ] )->at(j);
    	    string fileToOpenTbi = fileToOpen+".tbi";

	    if(!isFile(fileToOpen)){
		cerr << "File "<<fileToOpen<<" does not exist"<<endl;
		exit(1);
	    }

	    if(!isFile(fileToOpenTbi)){
		cerr << "Tabix "<<fileToOpenTbi<<" does not exist, files must be bgzipped and tabix indexed"<<endl;
		exit(1);
	    }
	    
	    if( strEndsWith(fileToOpen,"vcf.gz")){
	    	AlleleInfoReader * vcfSMP     = new VCFreader (      fileToOpen, fileToOpenTbi, grc.getChrName(), int(grc.getStartCoord()), int(grc.getEndCoord()) );
		vecOfAlleleReaders.push_back(vcfSMP);
		filesNames.   push_back(fileToOpen);
		filesNamestbi.push_back(fileToOpenTbi);
		filesNamesType.push_back(1);
		popIndexForFile.push_back(i);
	    }else{
	    	if( strEndsWith(fileToOpen,"bed.gz")){
	    	    AlleleInfoReader * vcfSMP = new BAMTABLEreader ( fileToOpen, fileToOpenTbi, grc.getChrName(), int(grc.getStartCoord()), int(grc.getEndCoord()) );
		    vecOfAlleleReaders.push_back(vcfSMP);
		    filesNames.   push_back(fileToOpen);
		    filesNamestbi.push_back(fileToOpenTbi);		    
		    filesNamesType.push_back(2);
		    popIndexForFile.push_back(i);
	    	}else{
	    	    cerr<<"ERROR: the sample must be a tabix indexed vcf or bamtable file"<<endl;
	    	    exit(1);
	    	}
	    }    
    	}
    }

    
    //detect file types, disable mixture of file types for a given population
    for(unsigned int indexAlleleFile=0;indexAlleleFile<vecOfAlleleReaders.size();indexAlleleFile++){
	if( popFileTypes[ popIndexForFile[indexAlleleFile] ]  == 0){
	    popFileTypes[ popIndexForFile[indexAlleleFile] ]  = filesNamesType[ indexAlleleFile ];
	}else{
	    if(popFileTypes[ popIndexForFile[indexAlleleFile] ]  != filesNamesType[ indexAlleleFile ] ){
		cerr<<"ERROR: cannot have a mixture of file types for "<< populations->at( popIndexForFile[indexAlleleFile] ) <<endl;
		exit(1);
	    }
	}
    }

    bool hasDataInAll=true;
    for(unsigned int indexAlleleFile=0;indexAlleleFile<vecOfAlleleReaders.size();indexAlleleFile++){
	hasDataInAll &= (vecOfAlleleReaders[indexAlleleFile]->hasData());//if one is false, the whole thing is false	
    }


    //no data to begin with for this loci, delete allocated variables
    if(!hasDataInAll){
	for(unsigned int indexAlleleFile=0;indexAlleleFile<vecOfAlleleReaders.size();indexAlleleFile++){
	    delete(vecOfAlleleReaders[indexAlleleFile]);
	}
	delete(rtEPO);
	return -1;
    }

   

    unsigned int coordCurrent=grc.getStartCoord();
    vector<AlleleInfo * > currentAlleleInfo;
    //All have data to return at this point 
    currentAlleleInfo.clear();
    for(unsigned int indexAlleleFile=0;indexAlleleFile<vecOfAlleleReaders.size();indexAlleleFile++){
	currentAlleleInfo.push_back( vecOfAlleleReaders[indexAlleleFile]->getData() );
    }


    if(currentAlleleInfo.size() != vecOfAlleleReaders.size()){
	cerr<<"Internal error, cannot get the same amount of data as in the input"<<endl;
	exit(1);
    }


    if(outputType == 2){
	printChrCoord=true; //we need the coordinate for mistar
	cout<<"#MISTAR"<<endl;    	
	cout<<"#PG:"<<programLine<<endl;
	cout<<"#GITVERSION: "<<gitHubVersion<<endl;
	cout<<"#DATE: "<<getDateString()<<endl;
	


	cout<<"#";
    }

	

    if(printChrCoord  ){ 
	cout<<"chr"<<fSeparator<<"coord"<<fSeparator;
    }

    if(outputType == 2){ 
	cout<<"REF,ALT\troot\tanc"<<fSeparator;
    }else{
	if(printChimp)    
	    cout<<"root"<<fSeparator;

    }
    cout<<vectorToString(*populations,fSeparator)<<endl;

    // for(int indexPopul=0; indexPopul < populations->size() ;indexPopul++){
    // 	cout<<populations->at(indexPopul)<<fSeparator;
    // }
    // cout<<endl;
    


    //MAIN LOOP OVER COORDINATES
    bool stayLoop=true;
    while(stayLoop){
	// cout<<endl<<"coordCurrent "<<coordCurrent<<endl;




	//determining position
	vector<string> fieldsEPO=allTokens(lineFromEPO,'\t');			
	unsigned int positionEPO=string2uint(fieldsEPO[1]);

	bool allSameCoord=true;
	allSameCoord &= (positionEPO ==  coordCurrent);
	for(unsigned int indexAlleleFile=0;indexAlleleFile<currentAlleleInfo.size();indexAlleleFile++){
	    allSameCoord &= (currentAlleleInfo[indexAlleleFile]->getPosition() ==  coordCurrent);
	}


	if(allSameCoord){	//we can safely enter

#ifdef DEBUG
	    // cout<<"allsame"<<endl;
	    cout<<lineFromEPO<<endl;
#endif
	    
	    bool  altAllelesFlags [5] = {false,false,false,false,false};  //will be true once we encounter certain alt alleles
	    int countPositiveAlt=0; //number of alternative alleles
	    int numberOfIndsWithAlternativeAllele=0; //number of individuals with the alternative allele
	    //BEGIN FILTERING VCF

	    //1,2,3,4 for A,C,G,T
	    int indexRefAllele;
	    int indexAltAllele=0;
	    int indexAltAlleleTemp;

	    //actual allele frequencies
	    vector<  pair<int,int>   > alleleCountByPop  (populations->size(), pair<int,int>(0,0)  );
	    vector<  bool   >  popPassed  (populations->size(), false  ); //set of flags to say at least one passed in that population
	    vector<  bool   >  cpGinPop   (populations->size(), false  ); //set of flags to flag cpg

	    string chimpInfoToprint     = "";
	    string ancestralInfoToprint = "";

	    bool   ancestralCpG=false;
	    int popFailedFilter; //quality filters
	    
	    int chimpAlleleIndex = -1;
	    int ancAlleleIndex   = -1;

	    indexRefAllele=  refBp2Index( ((SimpleVCF *)currentAlleleInfo[0])->getRef()[0] );
	    if(indexRefAllele == -1)
		goto nextiteration;
#ifdef DEBUG

	    cout<<"Ref fine "<<lineFromEPO<<endl;
#endif





	    //BEGIN COUNTING POTENTIAL ALTERNATIVE ALLELES AND VCF FILTERING
	    //filter VCF and detect the number of alternative alleles
	    //we only look at VCF for alternative alleles
	    for(unsigned int indexAlleleFile=0;indexAlleleFile<currentAlleleInfo.size();indexAlleleFile++){

		cpGinPop[  popIndexForFile[indexAlleleFile] ]  =  cpGinPop[ popIndexForFile[indexAlleleFile] ]  || currentAlleleInfo[indexAlleleFile]->isCpg();
	
		if(filesNamesType[indexAlleleFile] == 1){ //only for VCFs
		    //                       ,
		    bool thisVCFPassed= passedFilters((SimpleVCF *)currentAlleleInfo[indexAlleleFile] ,  //vcf
						      (*vcfCutoffs)[indexAlleleFile]);
		    // coverageCutoffs->at(indexAlleleFile).first,   //int minCovcutoff
		    // coverageCutoffs->at(indexAlleleFile).second,  //int maxCovcutoff
		    // minMapabilitycutoff,                          // double minMapabilitycutoff
		    // minMQcutoff,                                  //int minMQcutoff
		    // minGQcutoff);                                 //int minGQcutoff


		    //This next line does not work because vector<bool> is a packed representation and it seems the C++ guys have not overloaded
		    //every operator.
		    //(popPassed[ popIndexForFile[indexAlleleFile] ])  ||= thisVCFPassed; //if one is true, all are true for this pop
		    popPassed[ popIndexForFile[indexAlleleFile] ]  =  popPassed[ popIndexForFile[indexAlleleFile] ] || thisVCFPassed;

		    if(!thisVCFPassed) //don't bother adding those that did not pass
			continue;

		    //add the count for VCF files in alleleCountByPop, for bamtable files, they are added below
		    addPair ( &(alleleCountByPop[ popIndexForFile[indexAlleleFile] ]) ,
			      (((SimpleVCF *)currentAlleleInfo[indexAlleleFile])->returnLikelyAlleleCountForRefAlt(minPLdiffind)) );

		    //check alt allele
		    indexAltAlleleTemp=altBp2Index(  ((SimpleVCF *)currentAlleleInfo[indexAlleleFile])->getAlt()[0]  ) ;
		    if(indexAltAlleleTemp == -1) //-1 is for not A,C,G,T or .
			goto nextiteration;


		    altAllelesFlags[ indexAltAlleleTemp   ]  = true ; //should be safe to use [0]		    
		}else{ //BAM table records

		    //BAM table records are assumed to have been quality filtered
		    popPassed[ popIndexForFile[indexAlleleFile] ]  =  true;

		    //check which potential alternative alleles are present, at this point we ignore if the record has the reference
		    for(int indexAltNuc=1;indexAltNuc<=4;indexAltNuc++){
			if(indexRefAllele != indexAltNuc){ //putative alternative allele, non-reference
			    //cout<<"indexAltNuc "<<indexAltNuc<<endl;
			    if(  ((BAMTableObj *)currentAlleleInfo[indexAlleleFile])->hasAllele(indexAltNuc) ){
				altAllelesFlags[ indexAltNuc   ]  = true ; 
			    }
			}
		    }

		}

	    }// end for(unsigned int indexAlleleFile=0;indexAlleleFile<currentAlleleInfo.size();indexAlleleFile++){

	    //END COUNTING POTENTIAL ALTERNATIVE ALLELES AND VCF FILTERING









	    //BEGIN COUNT ALTERNATIVE ALLELES
	    for(int i=1;i<=4;i++){ //don't care about '.' as Alt
		if( altAllelesFlags[i] ){ //detected in either file
		    indexAltAllele=i;
		    countPositiveAlt++;
		}
	    }
	    

	    //if we do not require only segregating sites and the number of 
	    //alternative allele is 0, we let it go through
	    if(!onlySegSites &&
	       countPositiveAlt == 0){
		//fine
		indexAltAllele= 0; //dummy value to produce 'N'
	    }else{//we want to produce only segragating sites

		if(countPositiveAlt != 1) //skip all sites with more than one alternative allele
		    goto nextiteration;

		//checking if alternative was observed in more than one individual
		//this can be done to minimize the noise in a low coverage sample to 
		//skip private mutations
		for(unsigned int indexAlleleFile=0;indexAlleleFile<currentAlleleInfo.size();indexAlleleFile++){
		    if( currentAlleleInfo[indexAlleleFile]->hasAllele( indexAltAllele ) )
			numberOfIndsWithAlternativeAllele++;
		}

		if(numberOfIndsWithAlternativeAllele == 0) {
		    cerr<<"The number of individuals with the alternative allele is zero, problem in the program"<<endl;
		    exit(1);
		}

		//skip all sites where alternative was observed less than minIndWithAlt individual
		//with minIndWithAlt == 2, we require at least two individuals
		//with minIndWithAlt == 1, we allow private mutations for an individual
		if(numberOfIndsWithAlternativeAllele < minIndWithAlt) 
		    goto nextiteration;  
	    }
	    //END COUNT ALTERNATIVE ALLELES	    






	    //BEGIN CHECKING BAMTABLE RECORDS
	    //if the BAMTable has more alleles than the reference or the alternative, we skip this site
	    for(unsigned int indexAlleleFile=0;indexAlleleFile<currentAlleleInfo.size();indexAlleleFile++){	
		if(filesNamesType[indexAlleleFile] == 2){ //bamtable obj


		    if(indexAltAllele == 0){ //if there is an no indexAltAllele to check
			if(!(  ((BAMTableObj *)currentAlleleInfo[indexAlleleFile])->hasOnlyThisAlleles(indexRefAllele) )){
			    cerr<<"The following bamtable object should not have more alleles than the reference, exiting "<<(*((BAMTableObj *)currentAlleleInfo[indexAlleleFile]))<<endl;
			    exit(1);
			}

			//double check...
			if( ((BAMTableObj *)currentAlleleInfo[indexAlleleFile])->hasAllele(indexRefAllele) )
			    alleleCountByPop[ popIndexForFile[indexAlleleFile] ].first  ++;
		    }else{
			//no triallelic sites allowed
			if(!(  ((BAMTableObj *)currentAlleleInfo[indexAlleleFile])->hasOnly2Alleles(indexRefAllele,indexAltAllele) ))		   
			    goto nextiteration;

			if( ((BAMTableObj *)currentAlleleInfo[indexAlleleFile])->hasAllele(indexRefAllele) )
			    alleleCountByPop[ popIndexForFile[indexAlleleFile] ].first  ++;
			if( ((BAMTableObj *)currentAlleleInfo[indexAlleleFile])->hasAllele(indexAltAllele) )
			    alleleCountByPop[ popIndexForFile[indexAlleleFile] ].second ++;
		    }


		}
	    }
	    //END CHECKING BAMTABLE RECORDS	    






	    
	    //BEGIN COUNT HOW MANY POPULATION FAILED THE FILTERS
	    popFailedFilter=0;

	    for(unsigned int indexPopul=0; indexPopul < populations->size() ;indexPopul++){
		//if no individual within the population passed or if the pop is empty (due to bad quality data and low difference in the PL fields)
		if(  (!popPassed[ indexPopul ]) || isPairEmpty(alleleCountByPop[ indexPopul ]) ){
#ifdef DEBUG
		    cout<<"pop failed : "<<populations->at(indexPopul)<<endl;
#endif

		    popFailedFilter++;
		}
	    }
#ifdef DEBUG

	    cout<<"popFailedFilter "<<popFailedFilter<<endl;
#endif
	    //END COUNT HOW MANY POPULATION FAILED THE FILTERS

	    //if more than allowedFailPop failed, continue to next site
	    if(popFailedFilter > allowedFailPop)
	    	goto nextiteration;
	    //END FILTERING VCF






	    //PASSED EVERYTHING, BEGIN PRINTING THE ALLELE COUNT



	    //ancestral allele
	    if(printChimp){
		char ancAlleleRaw;
		char chimpAlleleRaw;

		//if(useEPOInferedAsAncestor){
		ancAlleleRaw   = fieldsEPO[3][0]; //use the Ortheus inference
		chimpAlleleRaw = fieldsEPO[4][0]; //use the chimp
		//		}
		


		if(fieldsEPO[9][0] == '1') //Ancestral is CpG
		    ancestralCpG=true;


		//look at ancestor
		if(ancAlleleRaw == 'N') //skip if the ancestral allele is not determined
		    goto ancestorfail;
		else{
		    ancAlleleIndex=refBp2Index( ancAlleleRaw );
		}

		if(ancAlleleIndex == -1 )//gaps and such
		    goto ancestorfail;
		

		if(ancAlleleIndex == indexRefAllele ){  //the reference is ancestral
		    ancestralInfoToprint="1,0";
		    goto epochimp; //fine
		}else{// the reference is not the (ancestral/chimp)

		    if(indexAltAllele == 0 ){ //this site is not a segregating site because there is not alternative allele, take the ancestor as the alternative allele
			//set the alternative allele to be the ancestral
			indexAltAllele=ancAlleleIndex;
		    }

		    if(ancAlleleIndex == indexAltAllele ){ //the alternative is ancestral
			ancestralInfoToprint="0,1";
			goto epochimp; //fine
		    }else{
			goto nextiteration; //skip if the ancestral allele is not the reference nor the alt allele, triallelic 
		    }
		    
		}

	    epochimp:

		//look at chimp
		if(chimpAlleleRaw == 'N') //skip if the chimp allele is not determined
		    goto ancestorfail;
		else{
		    chimpAlleleIndex=refBp2Index( chimpAlleleRaw );
		}

		if(chimpAlleleIndex == -1 ) //gaps and such
		    goto ancestorfail;		


		if(chimpAlleleIndex == indexRefAllele ){  //the reference is ancestral
		    chimpInfoToprint="1,0";
		    goto printremaining; //fine
		}else{// the reference is not the (ancestral/chimp)

		    if(indexAltAllele == 0 ){ //this site is not a segregating site because there is not alternative allele, take the chimp as the alternative allele
			//set the alternative allele to be the chimp
			indexAltAllele=chimpAlleleIndex;
		    }

		    if(chimpAlleleIndex == indexAltAllele ){ //the alternative is ancestral
			chimpInfoToprint="0,1";
			goto printremaining; //fine
		    }else{
			goto nextiteration; //skip if the chimp allele is not the reference nor the alt allele, triallelic
		    }
		    
		}



	    ancestorfail:
		if(outputType == 2){ //do not care for mistar but for other modes we do		  
		    chimpInfoToprint     = "0,0";
		    ancestralInfoToprint = "0,0";
		    goto printremaining; //fine
		}else{
		    goto nextiteration;
		}

	    }else{// do not print chimp/ancestor

		if(printChrCoord){ 	    //we print the chr and coordinate for tabix indexing if need be
		    cout<<grc.getChrName()<<fSeparator<<coordCurrent<<fSeparator;
		}

	    }	    

	printremaining:
	    if(printChrCoord){ 	    //we print the chr and coordinate for tabix indexing if need be
		cout<<grc.getChrName()<<fSeparator<<coordCurrent<<fSeparator;
	    }

	    //TREEMIX output
	    if(outputType == 1){ //

		if(useEPOInferedAsAncestor)
		    cout<<ancestralInfoToprint<<fSeparator;
		else
		    cout<<chimpInfoToprint<<fSeparator;

		for(unsigned int indexPopul=0; indexPopul < populations->size() ;indexPopul++){		    
		    if(indexPopul == (populations->size()-1))
			cout<<alleleCountByPop[indexPopul ].first<<","<<alleleCountByPop[indexPopul ].second;
		    else
			cout<<alleleCountByPop[indexPopul ].first<<","<<alleleCountByPop[indexPopul ].second<<fSeparator;
		}
		cout<<endl;
		goto nextiteration;
	    }

	    //mistar output
	    if(outputType == 2){ 

		cout<<("NACGT"[indexRefAllele])<<","<<("NACGT"[indexAltAllele])<<fSeparator;
		//if(printChimp)
		//cpg for both
		cout<<chimpInfoToprint    <<":"<<ancestralCpG<<fSeparator;
		cout<<ancestralInfoToprint<<":"<<ancestralCpG<<fSeparator;

		for(unsigned int indexPopul=0; indexPopul < populations->size() ;indexPopul++){
		    if(indexPopul == (populations->size()-1))
			cout<<alleleCountByPop[ indexPopul ].first<<","<<alleleCountByPop[indexPopul ].second<<":"<<cpGinPop[ indexPopul ];
		    else
			cout<<alleleCountByPop[ indexPopul ].first<<","<<alleleCountByPop[indexPopul ].second<<":"<<cpGinPop[ indexPopul ]<<fSeparator;
		}
		
		cout<<endl;
		goto nextiteration;
	    }



	    cerr<<"Invalid output type = "<<outputType<<endl;
	    exit(1);
	    //END printing


	    
	    

	nextiteration:
	    
	    //increase for next iteration
	    coordCurrent++;

	    lineLeftEPO=(rtEPO->readLine( lineFromEPO ));

	    if(!lineLeftEPO ){//reached end chr
		stayLoop=false;
		break;
	    }
		
	    for(unsigned int indexAlleleFile=0;indexAlleleFile<vecOfAlleleReaders.size();indexAlleleFile++){
		if(vecOfAlleleReaders[indexAlleleFile]->hasData()){
		    currentAlleleInfo[indexAlleleFile]=vecOfAlleleReaders[indexAlleleFile]->getData();
		}else{ //reached the end
		    stayLoop=false;
		    break;
		}
	    }



	}
	//if not all same coord, repositing the files to the correct position
	else{ 
	    
	    //BEGIN REPOSITIONING FILE POINTERS 

	    if(coordCurrent < positionEPO){ //overshot , repositioning there
		coordCurrent = positionEPO;
		continue;
	    }

	    if(coordCurrent == positionEPO){ //fine
	    }

	    if(coordCurrent > positionEPO){ //running behind
		if( (coordCurrent - positionEPO ) >= limitToReOpenFP){ //seeking in the file
		    // delete(rtEPO);
		    // rtEPO = new ReadTabix( epoFile.c_str()  , epoFileidx.c_str()  , grc.getChrName(), coordCurrent, grc.getEndCoord() ); 		    
		    rtEPO->repositionIterator(grc.getChrName(), coordCurrent, grc.getEndCoord() ); 		    
		}
		//retrieving new data
		lineLeftEPO=(rtEPO->readLine( lineFromEPO ));	    
	    }

	    if(!lineLeftEPO ){//reached end chr
		stayLoop=false;
		break;
	    }

	    // cout<<"positionEPO2 "<<positionEPO<<endl;

	    bool continueMain=false;

	    for(unsigned int indexAlleleFile=0;indexAlleleFile<vecOfAlleleReaders.size();indexAlleleFile++){
		// cout<<indexAlleleFile<<":"<<currentAlleleInfo[indexAlleleFile]->getPosition();

		if( coordCurrent == currentAlleleInfo[indexAlleleFile]->getPosition()   ){ //fine, do nothing
		    continue;
		}

		if( coordCurrent < currentAlleleInfo[indexAlleleFile]->getPosition()   ){//overshot , repositioning there
		    coordCurrent = currentAlleleInfo[indexAlleleFile]->getPosition();
		    continueMain=true;
		    continue;
		}


		if(coordCurrent > currentAlleleInfo[indexAlleleFile]->getPosition() ){ //running behind

		    if( (coordCurrent - currentAlleleInfo[indexAlleleFile]->getPosition() ) >= limitToReOpenFP){ //seeking in the file
			vecOfAlleleReaders[indexAlleleFile]->repositionIterator( grc.getChrName(), int(coordCurrent), int(grc.getEndCoord()) );
		    }else{
			//nothing to do
		    }
		    
		    //retrieving new data
		    if(vecOfAlleleReaders[indexAlleleFile]->hasData()){
			currentAlleleInfo[indexAlleleFile]=vecOfAlleleReaders[indexAlleleFile]->getData();
		    }else{ //reached the end for one
			stayLoop=false;
			break;
		    }
		    
		}
		// cout<<indexAlleleFile<<":"<<currentAlleleInfo[indexAlleleFile]->getPosition()<<"\t";
	    }//end for each vector allele
	    	    
	    //END REPOSITIONING FILE POINTERS 
	    if(continueMain)
		continue;

	}// else of if coordinates not equal
	
	
    } //end main loop


    

    //Clean up
    for(unsigned int indexAlleleFile=0;indexAlleleFile<vecOfAlleleReaders.size();indexAlleleFile++){
	delete(vecOfAlleleReaders[indexAlleleFile]);
    }

    delete(rtEPO);
    return 0;

}
