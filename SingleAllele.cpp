/*
 * SingleAllele
 * Date: Mar-18-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include "SingleAllele.h"

SingleAllele::SingleAllele(){
    refCount=0;
    altCount=0;
    isCpg=false;
}


SingleAllele::SingleAllele(int refCount,int altCount,bool isCpg):
    refCount(refCount),
    altCount(altCount),
    isCpg(isCpg){

}


SingleAllele::SingleAllele(const SingleAllele & other):
    refCount(other.refCount),
    altCount(other.altCount),
    isCpg(other.isCpg){

}

SingleAllele::~SingleAllele(){

}



int SingleAllele::getRefCount(){
    return refCount;
}

int SingleAllele::getAltCount(){
    return altCount;
}

bool SingleAllele::getIsCpg(){
    return isCpg;
}


void SingleAllele::setRefCount(int refCount_){
    refCount = refCount_;
}

void SingleAllele::setAltCount(int altCount_){
    altCount = altCount_;
}

void SingleAllele::setIsCpg(bool isCpg_){
    isCpg    = isCpg_;
}

string  SingleAllele::toString(){
    string toReturn =""+ 
	stringify( refCount)+","+
	stringify( altCount)+":"+
	stringify( isCpg);
    return toReturn;
}

SingleAllele operator+(const SingleAllele & first,const SingleAllele & second){

    return SingleAllele( first.refCount + second.refCount, 
    			 first.altCount + second.altCount,
    			 (first.isCpg   ||  second.isCpg)
    			 );
}


SingleAllele & SingleAllele::operator+=(const SingleAllele & other){
    this->refCount += other.refCount;
    this->altCount += other.altCount;
    this->isCpg    =  (this->isCpg   ||  other.isCpg);
    return *this;
}


bool operator== (const SingleAllele & first,const SingleAllele & second){
    return (  (first.refCount ==  second.refCount) &&
	      (first.altCount ==  second.altCount) &&
	      (first.isCpg    ==  second.isCpg)  );	      
}


bool operator!= (const SingleAllele & first,const SingleAllele & second){
    return !(first==second);
}
