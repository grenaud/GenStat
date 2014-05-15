/*
 * MultiVCFParser
 * Date: Oct-04-2012 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef MultiVCFParser_h
#define MultiVCFParser_h

#include <stdlib.h>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <map>

#include "SetVCFFilters.h"
#include "FilterVCF.h"
#include "SimpleVCF.h"
#include "ReadTabix.h"
#include "VCFreader.h"
#include "BAMTABLEreader.h"
#include "AlleleCounter.h"
#include "GenomicRange.h"

using namespace std;


class MultiVCFParser{
private:

    inline void addPair( pair<int,int> * destination, pair<int,int>  toAdd);
    inline bool isPairEmpty( pair<int,int>  & toCheck);

    //map name of populations to files
    map< string,vector<string> * >         * pop2PrefixFiles;   //map population name to vector of string representing the prefix of the files
    map< string,vector<string> * >         * pop2SuffixFiles;   //map population name to vector of string representing the prefix of the files
    
    vector<string>                       * populations; //vector of string for every population
    vector<SetVCFFilters *>              * vcfCutoffs;  //vector of SetVCFFilters for every vcf file

    // vector< pair<int,int> >              * coverageCutoffs; //vector of cutoffs for every vcf file
    inline int refBp2Index(const char toCheck);
    inline int altBp2Index(const char toCheck);

    bool printChrCoord; //print the chr[tab]coordinate ?
    bool printChimp   ; //print the chimp allele ?
    int  minPLdiffind ; //minimum pl difference 
    int  minIndWithAlt ;       //minimum # of individuals for private mutation
    int allowedFailPop;
public:
    //MultiVCFParser(string popfile);
    //MultiVCFParser(string popfile,bool printChimp,int minPLdiffind,bool printChrCoord,int minind,int allowedFailPop,int minGQcutoff,int minMQcutoff,double minMapabilitycutoff);
    MultiVCFParser(string popfile,bool printChimp,int minPLdiffind,bool printChrCoord,int minind,int allowedFailPop,
		   int    minGQcutoff         ,
		   int    minMQcutoff         ,
		   double minMapabilitycutoff ,
		   bool   filterIndelProx     ,
		   bool   repeatMasking       ,
		   bool   systemError     ,
		   bool   allowall,
		   bool   allowallMQ);
    ~MultiVCFParser();


    int produceOutput(string epoFile,string epoFileidx,GenomicRange grc,int outputType,bool onlySegSites,bool useEPOInferedAsAncestor,string programLine="",string gitHubVersion="NA");
    void printCutoffs();

};
#endif
