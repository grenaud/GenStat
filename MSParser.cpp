/*
 * MSParser
 * Date: Nov-02-2012 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include "MSParser.h"

MSParser::MSParser(string filename){
    string line;
    igzstream myFile;
    mso=new vector<MSobject *>();

    myFile.open(filename.c_str(), ios::in);
    bool inABlock=false;
    vector<string> linesBlock;
    if (myFile.good()){
	while ( getline (myFile,line)){
	   // cout<<inABlock <<"\t"<<line<<endl;
	    if(inABlock ){
		if(line == ""){
		    inABlock=false;
		    //MSobject mso ( &linesBlock);
		    mso->push_back( new MSobject ( &linesBlock) );
		}else{
		    linesBlock.push_back(line);
		}
	    }

	    if(line == "//"){
		inABlock=true;
		linesBlock.clear();
	    }
	}
	myFile.close();
	//last one
	//MSobject mso ( &linesBlock);
	mso->push_back( new MSobject ( &linesBlock) );
    }else{
	cerr << "Unable to open file "<<filename<<endl;
	exit(1);
    }

}

MSParser::~MSParser(){
    for(unsigned int i =0;i<mso->size();i++){
	delete( mso->at(i) );
    }
    delete(mso);
}

int MSParser::numberOfRecords() const{
    return int(mso->size());
}


const MSobject * MSParser::getMSObj(int indexObj) const{
    if(indexObj<0 || (indexObj >= int(mso->size()) ) ){
	cerr<<"Msparser cannot return object #"<< indexObj <<" must be between 0 and  "<<(mso->size()-1)<<" objects"<<endl;	    
	exit(1);
    }
    return mso->at(indexObj);
}

const vector<MSobject *> * MSParser::getVectorMSObj() const{
    return mso;
}
