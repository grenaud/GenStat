/*
 * testMS
 * Date: Nov-02-2012 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include <iostream>
#include <fstream>


#include "MSParser.h"
#include "utils.h"

using namespace std;

int main (int argc, char *argv[]) {

    MSParser ms (argv[1]);
    cout<<ms.numberOfRecords()<<endl;
    
    return 0;
}

