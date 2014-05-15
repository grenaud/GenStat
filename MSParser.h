/*
 * MSParser
 * Date: Nov-02-2012 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef MSParser_h
#define MSParser_h

#include <iostream>
/* #include <fstream> */
#include <gzstream.h>

#include "MSobject.h"

using namespace std;

class MSParser{
private:
    vector<MSobject *> * mso;

public:
    MSParser(string filename);
    ~MSParser();
    const vector<MSobject *> * getVectorMSObj() const;
    const MSobject * getMSObj(int indexObj) const;
    int numberOfRecords() const;
};
#endif
