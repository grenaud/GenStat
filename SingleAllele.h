/*
 * SingleAllele
 * Date: Mar-18-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef SingleAllele_h
#define SingleAllele_h

#include <iostream>
#include "utils.h"

using namespace std;

class SingleAllele{
private:
    int refCount;
    int altCount;
    bool isCpg;
public:
    SingleAllele();
    SingleAllele(int refCount,int altCount,bool isCpg);

    SingleAllele(const SingleAllele & other);
    ~SingleAllele();
    SingleAllele & operator= (const SingleAllele & other);

    int getRefCount();
    int getAltCount();
    bool getIsCpg();


    void setRefCount(int refCount);
    void setAltCount(int altCount);
    void setIsCpg(bool isCpg);
    string toString();


    friend bool operator== (const SingleAllele & first,const SingleAllele & second);
    friend bool operator!= (const SingleAllele & first,const SingleAllele & second);

    friend ostream& operator<<(ostream& os, const SingleAllele & at){
	os<<at.refCount<<","<<at.altCount<<":"<<at.isCpg<<endl;
	return os;
    }


    friend SingleAllele   operator+(const SingleAllele & first,const SingleAllele & second);
     SingleAllele & operator+=(const SingleAllele & other);

    //friend SingleAllele operator+(const SingleAllele & other);/* { */
    /* 	return SingleAllele( refCount + other.refCount,  */
    /* 			     altCount + other.altCount, */
    /* 			     (isCpg ||  other.isCpg) */
    /* 			     ); */
    /* } */
    
};
#endif
