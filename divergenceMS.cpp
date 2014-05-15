/*
 * DivergenceMS
 * Date: Nov-05-2012 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include <iostream>
#include <fstream>
#include <cstring>

#include "MSobject.h"
#include "MSParser.h"
#include "utils.h"

using namespace std;

int main (int argc, char *argv[]) {

   // string line;
   // ifstream myFile;
   // string filename = string(argv[1]);
   // myFile.open(filename.c_str(), ios::in);

   // if (myFile.is_open()){
   //   while ( getline (myFile,line)){

   //   }
   //   myFile.close();
   // }else{
   //     cerr << "Unable to open file "<<filename<<endl;
   //     return 1;
   //  }

    int referenceIndex = -1;
    int sampleIndex    = -1;
    int ancestorIndex  = -1;

    string usage=string(""+string(argv[0])+" [options] [ms file]\n\n"+
			"Please specify the index (1-based) of the following individual\n"+
			"the simulation:\n"+
			"\t-r reference index\n"+
			"\t-s sample index\n"+
			"\t-a ancestor index (or none for 0=ancestral)\n");

    if(argc == 1 ||
       (argc == 2 && (string(argv[2]) == "-h" || string(argv[2]) == "--help") )
       ){
	cerr << "Usage "<<usage<<endl;
	return 1;       
    }


    for(int i=1;i<(argc);i++){ 
	if(strcmp(argv[i],"-r") == 0 ){
	    referenceIndex=destringify<int>(argv[i+1]);
            continue;
        }
	if(strcmp(argv[i],"-s") == 0 ){
	    sampleIndex=destringify<int>(argv[i+1]);
            continue;
        }
	if(strcmp(argv[i],"-a") == 0 ){
	    ancestorIndex=destringify<int>(argv[i+1]);
            continue;
        }

    }

    
    if(referenceIndex == -1 ){
	cerr << "Error: the index of the reference must be specified using the -r option"<<endl;
	return 1;       
    }

    if(sampleIndex == -1 ){
	cerr << "Error: the index of the sample must be specified using the -s option"<<endl;
	return 1;       
    }

    // if(ancestorIndex == -1 ){
    // 	cerr << "Error: the index of the ancestor must be specified using the -a option"<<endl;
    // 	return 1;       
    // }

    MSParser ms (argv[ argc -1 ]); //last argument
    
    //    cout<<ms.numberOfRecords()<<endl;
    for(int indexRecord = 0 ;
	indexRecord<ms.numberOfRecords();
	indexRecord++){
	const MSobject * foundMSO=ms.getMSObj( indexRecord );
	const bool * refallele=foundMSO->getBoolArrayIndividual(referenceIndex);
	const bool * smpallele=foundMSO->getBoolArrayIndividual(   sampleIndex);
	const bool * ancallele;
	if(ancestorIndex != -1 ){
	    ancallele=foundMSO->getBoolArrayIndividual( ancestorIndex);
	}
	
	int noMutation       =0;
	int commonMutation   =0;
	int referenceMutation=0;
	int sampleMutation   =0;

	for(unsigned int indexSegsite=0;
	    indexSegsite<foundMSO->getNumberSegSites();
	    indexSegsite++){
	    //cout<<refallele[ indexSegsite ]<<smpallele[ indexSegsite ]<<ancallele[ indexSegsite ]<<endl;
	    bool ancestralAllele;
	    if(ancestorIndex == -1 ){
		ancestralAllele = false;
	    }else{
		ancestralAllele = ancallele[ indexSegsite ];
	    }

	    //no mutation
	    if(refallele[ indexSegsite ] == smpallele[ indexSegsite ]  &&
	       refallele[ indexSegsite ] == ancestralAllele ){
		noMutation++;
		continue;
	    }

	    //common branch
	    if(refallele[ indexSegsite ] == smpallele[ indexSegsite ] &&
	       refallele[ indexSegsite ] != ancestralAllele &&
	       smpallele[ indexSegsite ] != ancestralAllele ){
		commonMutation++;
		continue;
	    }
	    
	    //reference branch
	    if(refallele[ indexSegsite ] != smpallele[ indexSegsite ] &&
	       refallele[ indexSegsite ] != ancestralAllele &&
	       smpallele[ indexSegsite ] == ancestralAllele ){
		referenceMutation++;
		continue;
	    }

	    //sample branch
	    if(refallele[ indexSegsite ] != smpallele[ indexSegsite ] &&
	       refallele[ indexSegsite ] == ancestralAllele &&
	       smpallele[ indexSegsite ] != ancestralAllele ){
		sampleMutation++;
		continue;
	    }

	    cerr << "Error: invalid state for :"<<endl;
	    cerr<<refallele[ indexSegsite ]<<smpallele[ indexSegsite ]<<ancestralAllele<<endl;
	    return 1;      
	    
	}//end for each seg site


	cout<<indexRecord<<"\t"<<noMutation<<"\t"<<commonMutation<<"\t"<<referenceMutation<<"\t"<<sampleMutation<<"\t"<<(double(referenceMutation)/double(referenceMutation+commonMutation))<<"\t"<<(double(sampleMutation)/double(sampleMutation+commonMutation))<<endl;
    }//end for each ms record


    return 0;
}

