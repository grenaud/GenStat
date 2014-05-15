/*
 * ComputeDivergence
 * Date: Aug-15-2012 
 * Author : Gabriel Renaud gabriel.reno@gmail.com
 *
 */

#ifndef ComputeDivergence_h
#define ComputeDivergence_h

#include "FilterVCF.h"
#include "DivergenceResult.h"
#include "SimpleVCF.h"
#include "ReadTabix.h"
#include "VCFreader.h"
#include "BAMTABLEreader.h"
//#include "AlleleCounter.h"
#include "GenomicRange.h"
#include "ComputeDivergence_core.h"


using namespace std;

int computeDivergence(string refereVCF,string refereVCFidx,string sampleVCF,string sampleVCFidx,bool   useBedFileRegions,vector<GenomicRange> *  bedRegionsToFilter,string epoFile,string epoFileidx,int minimumNumberOfGoodSites,string chrName, int startChrCoord,  int endChrCoord,SetVCFFilters * filtersVCFREF,SetVCFFilters * filtersVCFSMP,int bpForIndels,bool cpgInData,bool useHumanGenomeAsReference,int minPLdiffind,bool maximizeDiv);


int computeDivergence(string refereVCF,string refereVCFidx,string sampleVCF,string sampleVCFidx,bool   useBedFileRegions,vector<GenomicRange> *  bedRegionsToFilter,string epoFile,string epoFileidx,int minimumNumberOfGoodSites,GenomicRange grc,SetVCFFilters * filtersVCFREF,SetVCFFilters * filtersVCFSMP,int bpForIndels,bool cpgInData,bool useHumanGenomeAsReference,int minPLdiffind,bool maximizeDiv);

#endif
