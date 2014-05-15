/*
 * AllelePairCounter
 * Date: Jan-26-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef AllelePairCounter_h
#define AllelePairCounter_h

#include <vector> 
#include <ostream> 
#include <iostream> 
#include <stdlib.h>

#include "utils.h"

using namespace std;


class AllelePairCounter{

private:
    unsigned int count[16];
public:
    void reinitializedCounters();
    void addAllelePair(int indexDimer);
    string headerForCount() const;
    unsigned int getIndent();
    unsigned int getMutations();

    AllelePairCounter();
    AllelePairCounter(const string stringForm);

    ~AllelePairCounter();
    friend std::ostream & operator<<(std::ostream & os, const AllelePairCounter & apc);
    AllelePairCounter &  operator+=(const AllelePairCounter & other);

};
#endif
