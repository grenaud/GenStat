/*
 * testRandomCoordGenerator
 * Date: Aug-17-2012 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include <iostream>
#include <fstream>

#include "RandomGenomicCoord.h"
#include "GenomicWindows.h"

using namespace std;

int main (int argc, char *argv[]) {


    RandomGenomicCoord rgc (string("index.hg19.fai"));
    GenomicWindows rw (string("index.hg19.fai"));

    if(0){
	for(int i=0;i<100;i++){
	    GenomicRange c=rgc.getRandomGenomic(10000);
	    cout<<c<<endl;
	    // randomCoord c=rgc.getRandomGenomic(10000);
	    // cout<<c.chrName<<"\t"<<c.chrCoord<<endl;
	}


    vector<GenomicRange> v=rw.getGenomicWindowsChr("22",100000,10000);
    for(unsigned int i=0;i<v.size();i++){    
	cout<<v[i]<<endl;
    }

    // return 1;
    vector<GenomicRange> v2=rw.getGenomicWindows(100000,10000);
    for(unsigned int i=0;i<v2.size();i++){    
	//cout<<v2[i]<<endl;
    }
    }

    vector<GenomicRange> v3=rw.getChr("2");
    for(unsigned int i=0;i<v3.size();i++){    
	cout<<v3[i]<<endl;
    }

    vector<GenomicRange> v4=rw.getGenomeWide();
    for(unsigned int i=0;i<v4.size();i++){    
	cout<<v4[i]<<endl;
    }



    return 0;
}

