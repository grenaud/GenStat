/*
 * Dstats
 * Date: Oct-15-2012 
 * Author : Gabriel Renaud gabriel.reno@gmail.com
 *
 */

#ifndef Dstats_h
#define Dstats_h

#include "FilterVCF.h"
#include "SimpleVCF.h"
#include "ReadTabix.h"
#include "VCFreader.h"
#include "BAMTABLEreader.h"
#include "DstatResult.h"
#include "GenomicRange.h"
#include "AlleleInfoReader.h"
#include "Dstat_core.h"

using namespace std;


int dstats(string refereVCF,
	   string refereVCFidx,
	   string sampleVCF,
	   string sampleVCFidx,
	   string conditionalVCF,
	   string conditionalVCFidx,
	   string epoFile,
	   string epoFileidx,
	   int minimumNumberOfGoodSites,
	   string chrName, 
	   int startChrCoord,  int endChrCoord,
	   //int minQCcutoff,int minCovcutoffREF,int maxCovcutoffREF,int minCovcutoffSMP,int maxCovcutoffSMP,
	   //int minCovcutoffCOND,int maxCovcutoffCOND,
	   //int minMQcutoff,double minMapabilitycutoff,
	   SetVCFFilters * filtersVCFREF,
	   SetVCFFilters * filtersVCFSMP,
	   SetVCFFilters * filtersVCFCOND,
	   int bpForIndels,bool cpgInData,bool useHumanGenomeAsReference,bool useEPOInferedAsAncestor,bool requiredHomozygousCond);

int dstats(string refereVCF,string refereVCFidx,string sampleVCF,string sampleVCFidx,
	   string conditionalVCF,
	   string conditionalVCFidx,string epoFile,string epoFileidx,int minimumNumberOfGoodSites,GenomicRange grc,
	   //int minQCcutoff,int minCovcutoffREF,int maxCovcutoffREF,int minCovcutoffSMP,int maxCovcutoffSMP,int minCovcutoffCOND,int maxCovcutoffCOND,int minMQcutoff,double minMapabilitycutoff,
	   SetVCFFilters * filtersVCFREF,
	   SetVCFFilters * filtersVCFSMP,
	   SetVCFFilters * filtersVCFCOND,
	   int bpForIndels,bool cpgInData,bool useHumanGenomeAsReference,bool useEPOInferedAsAncestor,bool requiredHomozygousCond);


#endif
