/*
 * AllelePairCounter
 * Date: Jan-26-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include "AllelePairCounter.h"


AllelePairCounter::AllelePairCounter(){
    reinitializedCounters();
}




AllelePairCounter::AllelePairCounter(const string stringForm){
    vector<string> allLines=allTokens(stringForm,'\t');
    for(int i=0;i<16;i++){
	count[i] = destringify<unsigned int>(allLines[i+1]);
    }
}




AllelePairCounter::~AllelePairCounter(){

}

void AllelePairCounter::reinitializedCounters(){
    for(int i=0;i<16;i++)
	count[i]=0;
}

void AllelePairCounter::addAllelePair(int indexDimer){
    count[indexDimer]++;
}


string AllelePairCounter::headerForCount() const{
    string dnaAlphabet = "ACGT";
    string toReturn="";

    for(int i=0;i<4;i++){
	for(int j=0;j<4;j++){
	    toReturn+=dnaAlphabet.substr(i,1)+dnaAlphabet.substr(j,1);
	    if(j!=3)
		toReturn+="\t";
	}
	if(i!=3)
	    toReturn+="\t";	
    }

    toReturn+="\tident\tmutat\trate";	

    return toReturn;
}

unsigned int AllelePairCounter::getIndent(){
    //unsigned int mutations=0;
    unsigned int ident    =0;
    for(int i=0;i<16;i++){
	if( i == 0 || i == 5 || i == 10 || i == 15 ){
	    ident+=     count[i] ;
	}else{
	    //mutations+= apc.count[i];	
	}
    }
    return ident;
}

unsigned int AllelePairCounter::getMutations(){
    unsigned int mutations=0;
    for(int i=0;i<16;i++){
	if( i == 0 || i == 5 || i == 10 || i == 15 ){
	    //ident+=     apc.count[i] ;
	}else{
	    mutations+= count[i];	
	}
    }
    return mutations;
}

std::ostream & operator << (std::ostream & os, const AllelePairCounter & apc){
    unsigned int mutations=0;
    unsigned int ident    =0;

    for(int i=0;i<16;i++){
	if( i == 0 || i == 5 || i == 10 || i == 15 ){
	    ident+=     apc.count[i] ;
	}else{
	    mutations+= apc.count[i];	
	}
    }

    for(int i=0;i<15;i++){
	os<<apc.count[i]<<"\t";
    }
    os<<apc.count[15]<<"\t"<<ident<<"\t"<<mutations<<"\t"<<double(mutations) /double(mutations+ident);

    return os;
}


AllelePairCounter &  AllelePairCounter::operator+=(const AllelePairCounter & other){

    for(int i=0;i<16;i++){
	this->count[i] += other.count[i];
    }

    return *this;

}
