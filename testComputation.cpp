/*
 * testRandomCoordGenerator
 * Date: Aug-17-2012 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include <iostream>
#include <fstream>

#include "ComputeFracAnc.h"
#include "ComputeFracHetero.h"
#include "SetVCFFilters.h"

using namespace std;

int main (int argc, char *argv[]) {



    // returncode= computeFracAnc("/mnt/454/HCNDCAM/1_Extended_VCF/HGDP00521/HGDP00521.hg19_1000g.1.mod.vcf.gz",
    // 				   "/mnt/454/HCNDCAM/1_Extended_VCF/HGDP00521/HGDP00521.hg19_1000g.1.mod.vcf.gz.tbi",
    // 				   "/mnt/454/Altaiensis/users/gabriel/EPOindices/chr1.epo.gz",
    // 				   "/mnt/454/Altaiensis/users/gabriel/EPOindices/chr1.epo.gz.tbi",
    // 				   30,
    // 				   "1",
    // 				   222217540,
    // 				   222227540,
    // 				   40,
    // 				   10,
    // 				   80,
    // 				   20,
    // 				   1,
    // 				   5);

    SetVCFFilters * filterToUse=new SetVCFFilters(40          ,
						  30          ,
						  1.0  ,
						  true ,
						  true  ,
						  true         ,
						  10,
						  80  ,
						  false,
						  false); 

    int returncode= computeFracHetero("/mnt/454/HCNDCAM/1_Extended_VCF/HGDP00521/HGDP00521.hg19_1000g.1.mod.vcf.gz",
				  "/mnt/454/HCNDCAM/1_Extended_VCF/HGDP00521/HGDP00521.hg19_1000g.1.mod.vcf.gz.tbi",
				  // "/mnt/454/Altaiensis/users/gabriel/EPOindices/chr1.epo.gz",
				  // "/mnt/454/Altaiensis/users/gabriel/EPOindices/chr1.epo.gz.tbi",			  
				  30,
				   "1",
				  222217540,
				  222297540,
				  filterToUse,
				  // 40,
				  // 10,
				  // 80,
				  // 20,
				  // 1,
				  5);

    cout<<returncode<<endl;
    delete(filterToUse);
    //string refereVCF,string refereVCFidx,string epoFile,string epoFileidx,int minimumNumberOfGoodSites,string chrName, int startChrCoord,  int endChrCoord,int minQCcutoff,int minCovcutoffREF,int maxCovcutoffREF,int minMQcutoff,double minMapabilitycutoff,int bpForIndels){

    
    // RandomGenomicCoord rgc (string("index.hg19.fai"));
    // GenomicWindows rw (string("index.hg19.fai"));

    // if(0)
    // for(int i=0;i<100;i++){
    // 	GenomicRange c=rgc.getRandomGenomic(10000);
    // 	cout<<c<<endl;
    // 	// randomCoord c=rgc.getRandomGenomic(10000);
    // 	// cout<<c.chrName<<"\t"<<c.chrCoord<<endl;
    // }

    // vector<GenomicRange> v=rw.getGenomicWindowsChr("22",100000,10000);
    // for(int i=0;i<v.size();i++){    
    // 	//cout<<v[i]<<endl;
    // }
    // // return 1;
    // vector<GenomicRange> v2=rw.getGenomicWindows(100000,10000);
    // for(int i=0;i<v2.size();i++){    
    // 	//cout<<v2[i]<<endl;
    // }

    return 0;
}

