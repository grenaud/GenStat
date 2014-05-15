/*
 * PairwiseDivergence
 * Date: Jan-24-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include "PairwiseDivergence.h"

void pairwiseDivergence(string filename){

    ifstream myFile;
    myFile.open(filetoread.c_str(), ios::in);
    string line;
    int numberPopulations=0;
    vector<string> populationNames;

    if (myFile.is_open()){
	//parse header
	if(!getline (myFile,line)){
	    cerr << "Unable to open file "<<filetoread<<endl;
	    exit(1);	    
	}else{
	    vector<string> fields=allTokens(line,'\t');
	    if(fields[0] != "chr")    { cerr<<"Field #1 of header must be chr ";     exit(1); }
	    if(fields[1] != "coord")  { cerr<<"Field #2 of header must be coord ";   exit(1); }
	    if(fields[2] != "REF,ALT"){ cerr<<"Field #3 of header must be REF,ALT "; exit(1); }
	    if(fields[3] != "root")   { cerr<<"Field #4 of header must be root ";    exit(1); }
	    for(int i=4;i<fields.size();i++){
		populationNames.push_back(fields[i]);
		numberPopulations++;
	    }
	}

	while ( getline (myFile,line)){	    
	    string chrName;
	    unsigned int startCoord;
	    unsigned int endCoord;
	    vector<string> temp=allTokens(line,'\t');
	    chrName     = destringify<string>(temp[0]);
	    startCoord  = destringify<unsigned int>(temp[1]);
	    endCoord    = destringify<unsigned int>(temp[2]);
	    GenomicRange toadd (chrName,startCoord,endCoord);
	    toReturn.push_back(toadd);
	}
	myFile.close();
    }else{
	cerr << "Unable to open file "<<filetoread<<endl;
	exit(1);
    }

    
}
