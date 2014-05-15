/*
 * DivergenceResult
 * Date: Jan-24-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef DivergenceResult_h
#define DivergenceResult_h

using namespace std;

#include "AlleleCounter.h"

class DivergenceResult{
 private:
 public:
    //This is everything
    AlleleCounter all;
    //This is without the ones marked as CpG
    AlleleCounter noCpg;
    //This is only with the ones marked as CpG
    AlleleCounter onlyCpg;
    //This excludes the following cases:
    // S = C, R or A = T
    // S = T, R or A = C
    // S = A, R or A = G
    // S = G, R or A = A    
    AlleleCounter transversions;


    AlleleCounter transitions;


    //This excludes the following cases:
    // S = T, R or A = C
    // S = A, R or A = G
    AlleleCounter noDamage;

    DivergenceResult();
    ~DivergenceResult();
    string getHeader();

    friend ostream& operator<<(ostream& os, const DivergenceResult & ct){
	os<<ct.all<<"\t"
	  <<ct.onlyCpg<<"\t"
	  <<ct.noCpg<<"\t"
	  <<ct.transitions<<"\t"
	  <<ct.transversions<<"\t"
	  <<ct.noDamage;
	return os;
    }
};
#endif
