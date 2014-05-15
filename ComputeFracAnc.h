/*
 * ComputeDivergence
 * Date: Aug-20-2012 
 * Author : Gabriel Renaud gabriel.reno@gmail.com
 *
 */

#ifndef ComputeFracAnc_h
#define ComputeFracAnc_h

#include "SetVCFFilters.h"
#include "FilterVCF.h"
#include "SimpleVCF.h"
#include "ReadTabix.h"
#include "VCFreader.h"
#include "GenomicRange.h"

using namespace std;

int computeFracAnc(string refereVCF,string refereVCFidx,string epoFile,string epoFileidx,int minimumNumberOfGoodSites,GenomicRange  grg,SetVCFFilters * filtersVCF, int bpForIndels);
		   //int minGQcutoff,int minCovcutoffREF,int maxCovcutoffREF,int minMQcutoff,double minMapabilitycutoff,int bpForIndels);

int computeFracAnc(string refereVCF,string refereVCFidx,string epoFile,string epoFileidx,int minimumNumberOfGoodSites,string chrName, int startChrCoord,  int endChrCoord,SetVCFFilters * filtersVCF, int bpForIndels);
		   //int minGQcutoff,int minCovcutoffREF,int maxCovcutoffREF,int minMQcutoff,double minMapabilitycutoff,int bpForIndels);
//void computeDivergence( string refereVCF,string refereVCFidx,string sampleVCF,string sampleVCFidx,string epoFile,string epoFileidx,int minimumNumberOfGoodSites,string chrName, int startChrCoord,  int endChrCoord, int bpForIndels);

#endif
