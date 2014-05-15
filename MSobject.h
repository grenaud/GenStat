/*
 * MSobject
 * Date: Nov-02-2012 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef MSobject_h
#define MSobject_h

#include <vector>
#include <string>

#include "utils.h"

using namespace std;

class MSobject{

private:
    int numberSegSites;
    int numberIndividual;
    //vector<double> positions; //works with custom ms
    vector<int> positions;
    bool ** alleleMatrix;


public:
    MSobject(vector<string> * inputString);
    ~MSobject();
    const bool * getBoolArrayIndividual(int indexInd) const;
    //const vector<double> * getPositions() const; //works with custom ms
    const vector<int> * getPositions() const;


    //bool * getBoolArraySite(); //to implement ? to retun a slice according to a site, will require memory allocation...
    unsigned int getNumberSegSites() const;
    int getNumberIndividuals() const;

};
#endif
