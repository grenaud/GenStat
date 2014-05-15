/*
 * DistanceResult
 * Date: Jan-26-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef DistanceResult_h
#define DistanceResult_h

#include <ostream> 

#include "DistanceResult.h"
#include "AllelePairCounter.h"

using namespace std;

class DistanceResult{
private:

public:
    //This is everything
    AllelePairCounter all;
    //This is without the ones marked as CpG
    AllelePairCounter noCpg;
    //This is only with the ones marked as CpG
    AllelePairCounter onlyCpg;
    //This is only with the ones marked as transitions
    AllelePairCounter transitions;
    //This is only with the ones marked as transversions
    AllelePairCounter transversions;


    DistanceResult();
    DistanceResult(const string stringForm);

    ~DistanceResult();
    friend std::ostream & operator<<(std::ostream & os, const DistanceResult & dr);
    DistanceResult &  operator+=(const DistanceResult & other);



};
#endif
