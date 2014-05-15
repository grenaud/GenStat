#include <iostream>
#include <fstream>
#include <sstream> 
#include <string> 
#include <vector> 
#include <map>
#include <algorithm>
#include <time.h>
#include <stdlib.h>


#include "utils.h"
#include "ComputeDivergence.h"
#include "ComputeFracAnc.h"
#include "ComputeFracHetero.h"
#include "RandomGenomicCoord.h"
#include "GenomicWindows.h"
#include "GenomicRange.h"

#include "MultiVCFParser.h"
#include "Dstats.h"
#include "SetVCFFilters.h"
#include "DivergenceResult.h"

//#define DEBUG
// #define DEBUG2
//#define DEBUG3

//TODO
//sum up the counters for genome-wide 
//using an overloaded sum
//Multiple coverage for D-stats ?





using namespace std;

//! Method to read a sorted bed file 
/*!
 *
 * This method checks for the records being ordered

  \param filetoread : String with the full path to the file to read
  \return           : Return(head) a pointer to a map where the key is the chromosome name and the value a vector of genomic ranges
  \sa  readBEDSortedfile()
*/

map< string, vector<GenomicRange> * > * readBEDSortedfile(string filetoread){
    //vector<GenomicRange> toReturn;
    map< string, vector<GenomicRange> * > * toReturn= new map< string, vector<GenomicRange> *>();
    ifstream myFile;
    myFile.open(filetoread.c_str(), ios::in);
    string line;
    unsigned int     lastEndCoord = 0;
    string           lastChrname  = "###";

    if (myFile.is_open()){
	while ( getline (myFile,line)){	    
	    //cout<<line<<endl;

	    string       chrName;
	    unsigned int startCoord;
	    unsigned int endCoord;
	    vector<string> temp=allTokens(line,'\t');
	    if(temp.size() != 3){
		cerr << "Error in vcfcompute in readBEDSortedfile(): following line does not have 3 fields"<< line<<endl;
		exit(1);		
	    }
	    chrName     = destringify<string>(temp[0]);
	    startCoord  = destringify<unsigned int>(temp[1])+1; //the left coordinate is zero based
	    endCoord    = destringify<unsigned int>(temp[2]);

	    if(lastChrname != chrName){//new chr found
		//check for previously found		
		lastChrname  = chrName;
		lastEndCoord = endCoord;
		if(toReturn->find(chrName) == toReturn->end() ){//not previously found
		    //cout<<"new chr "<<chrName<<endl;
		    (*toReturn)[chrName]=new vector<GenomicRange>();
		}else{
		    cerr << "Cannot have multiple records on the same chromosome at different parts of the file "<< line<<endl;
		    exit(1);
		}
	    }else{//stay on same chr
		if(startCoord <= lastEndCoord ){
		    cerr << "Problem with line =  "<<line<<" the start of the coordinate lesser than the end of the previous record "<<lastEndCoord<<endl;
		    exit(1);
		}
		if(endCoord   <= lastEndCoord ){
		    cerr << "Problem with line =  "<<line<<" the end of the coordinate lesser than the end of the previous record "<<lastEndCoord<<endl;
		    exit(1);		    
		}
	    }

	    GenomicRange toadd (chrName,startCoord,endCoord);	    
	    (*toReturn)[chrName]->push_back(toadd);
	    //cout<<(*toReturn)[chrName]->size()<<endl;
	}
	myFile.close();
    }else{
	cerr << "Unable to open file "<<filetoread<<endl;
	exit(1);
    }

    return toReturn;
}


//! Simple method to read a bed file 
/*!
 *
 * This method does not check for the records being ordered

  \param filetoread : String with the full path to the file to read
  \return           : Return(stack) the  value of a vector of GenomicRange objects
  \sa  readBEDSortedfile()
*/

vector<GenomicRange> readBEDfile(string filetoread){
    vector<GenomicRange> toReturn;
    ifstream myFile;
    myFile.open(filetoread.c_str(), ios::in);
    string line;

    if (myFile.is_open()){
	while ( getline (myFile,line)){	    
	    string chrName;
	    unsigned int startCoord;
	    unsigned int endCoord;
	    vector<string> temp=allTokens(line,'\t');
	    if(temp.size() != 3){
		cerr << "Error in vcfcompute in readBEDSortedfile(): following line does not have 3 fields"<< line<<endl;
		exit(1);		
	    }

	    chrName     = destringify<string>(temp[0]);
	    startCoord  = destringify<unsigned int>(temp[1])+1; //the left coordinate is zero based
	    endCoord    = destringify<unsigned int>(temp[2]);
	    GenomicRange toadd (chrName,startCoord,endCoord);
	    toReturn.push_back(toadd);
	}
	myFile.close();
    }else{
	cerr << "Unable to open file "<<filetoread<<endl;
	exit(1);
    }

    return toReturn;
}

// string booleanAsString(bool toprint){
//     if(toprint)
// 	return string("currently turned on/used");
//     else
// 	return string("not on/not used");
// }


void setVectorUsedAsString(string arg,vector<int> * vectorToset){
    vector<string> tempStr=allTokens( arg, ',' );
    vectorToset->clear();
    for(unsigned int ind=0;ind<tempStr.size();ind++)
	vectorToset->push_back( destringify<int>(tempStr[ind]) );
}

int main (int argc, char *argv[]) {
 
    ///////////////////////////////
    //
    //   BEGIN DEFINE VARIABLES
    //
    ///////////////////////////////


    //quality cutoffs
    vector<SetVCFFilters *> * filtersToUse = new vector<SetVCFFilters *>();
    int minGQcutoff=0;
    int minMQcutoff=30;
    double minMapabilitycutoff=1.0;
    int bpForIndels=5;
    bool allowCloseIndelProx = false;
    bool allowRepeatMasking  = false;
    bool allowSysErr         = false;
    bool allowall            = false;
    bool allowallMQ          = false;

    
    int minCovcutoffREF=0;
    int maxCovcutoffREF=1000;
    int minCovcutoffSMP=0;
    int maxCovcutoffSMP=1000;
    int minCovcutoffCOND=0;
    int maxCovcutoffCOND=1000;
    

    char mode; //d = divergence, a=ancall, h = hetero, t=treemix,  s=d-stats

    bool allowSexChr=false;
    bool useHumanGenomeAsReference=false;
    bool useEPOInferedAsAncestor=false;
    int overlapBetweenWindows=0;

    //treemix options
    bool printChrCoord=false;
    bool printChimp   =true;

    int    minPLdiffind=33;
    string popfile="";
    int outtype=-1;
    bool requiredHomozygousCond=true;
    int minIndWithAlt=1;
    int allowedFailPop=0;
    bool onlySegSites=true;

    // case 'x':  //undefined
    // case 'f':  //reading regions from a file
    // case 'r':  //random windows
    // case 'o':  //entire genome per window
    // case 'c':  //entire chr per window	    
    // case 'w':  //genome wide 
    // case 'v':  //entire chr (only one window)
    char lociSelection='x';
    string chrToUse;
    string bedFileWithRegions;

    string epoIndex      ;
    string fastaIndex    ;

    string sampleVCFSuf  ;
    string sampleVCFPre  ;

    string refereVCFSuf  ;
    string refereVCFPre  ;

    string condVCFSuf  ;
    string condVCFPre  ;


    unsigned int amountOfGoodSitesTARGET=0;  //number of contiguous chunks to extract for random

    int bpToExtract =-1      ; //amount of bp per chunk
    double fractionGood=-1.0;
    string bedFileRegions;
    string bedFileRegionsTabix;

    bool   useBedFileRegions     =false;
    bool   useBedFileRegionsTabix=false;

    map< string, vector<GenomicRange> * > * bedRegionsToFilter;


    int minimumNumberOfGoodSites; //minimum number of sites that pass qual filters
    bool cpgData=false;
    bool maximizeDiv=false;

    MultiVCFParser * tmim;

    ReadTabix * bedFileTabix=0;


    ///////////////////////////////
    //
    //    END DEFINE VARIABLES
    //
    ///////////////////////////////




    ///////////////////////////////
    //
    //   BEGIN DEFINE OPTIONS
    //
    ///////////////////////////////

    string usage=string(""+string(argv[0])+"  [mode]"+
			"\nThis program has the following modes:\n\n"+
			//"\t\t"+"pairdivcompute  To compute divergence pairwise for mistar input\n"+
			"\t\t"+"divcompute      To identify in which lineage the mutations occured for two files\n"+
			"\t\t"+"ancall          To compute the percentage of ancestral alleles\n"+
			"\t\t"+"hetero          To compute the percentage of heterozygous sites\n"+
			"\t\t"+"dstat           To compute D-statistics    \n"+
			// "\t\t"+"dstatcond   Conditional d-stats (on homozygous derived sites in the conditional VCF)\n"+
			"\t\t"+"multivcf        Parse multiple vcfs To produce treemix/mistar input \n");

    string selection=string("")+
	"\tLoci selection mode, select either:\n"+
	"\t\t"+"--region  [file]"   +"\t\t\t"     +"Read [file] in BED format to select which loci to use\n"+
	"\t\t"+"--random  [number]" +"\t\t\t"     +"Select [number] random loci from the genome\n"+
	"\t\t"+"--loci"             +"\t\t\t\t\t" +"Select multiple loci from the genome\n"+
	"\t\t"+"--locichr [chr]"    +"\t\t\t\t"   +"Select multiple loci from a chromosome\n"+
	"\t\t"+"--chr     [chr]"    +"\t\t\t\t"   +"Run on an entire chromosome\n"+
	"\t\t"+"--wide"             +"\t\t\t\t\t" +"Run on the entire genome\n";

    string multivcf=string("")+
	"\tMultivcf options:\n"+
	"\t\t"+"--outtype [type]" +"\t\t\t"+"Output type:\n"+

	"\t\t\t1:"+"\t"+"treemix\n"+	
	"\t\t\t2:"+"\t"+"mistar\n"+	
	"\t\t"+"--popfile [path]" +"\t\t\t"+"File containing the populations to use (default: none)\n"+
	"\n\t\t\tThis file must have the following format:\n"+
	"\n\t\t\tPopulation1:\n"+
	"\t\t\t/path/to/vcffile1 .vcf.gz\n"+
	"\t\t\t/path/to/vcffile2 .vcf.gz\n"+
	"\t\t\tPopulation2:\n"+
	"\t\t\t/path/to/vcffile3 .vcf.gz\n\n"+
	"\t\t\tWhere the first part is the prefix for each file and the remaining part the suffix\n"+
	"\t\t\tCoverage boundaries can be specified as two integer at the end of a VCF file as such:\n"+
	"\t\t\t/path/to/vcffile3 .vcf.gz 10 40\n\n"+
	"\t\t\tHence, any site in /path/to/vcffile3 .vcf.gz not having coverage between 10 and 40 will be discarded\n\n"+       
	"\t\t\tThis file is space separated, no tabs\n\n"+       

	"\t\t"+"--nochimp"          +"\t\t\t\t" +"Do not print chimp allele                                                 (default= "+booleanAsString(printChimp)+")"+"\n"+ 
	// "\t\t"+"--anc"              +"\t\t\t\t\t" +"Do not print chimp allele                                                 (default= "+booleanAsString(printChimp)+")"+"\n"+ 

	//"\t\t"+"--pldiffind [value]"+"\t\t\t"     +"Minimum difference in PL value for a confident allele score               (default= "+stringify(minPLdiffind)+")"+"\n"+
	"\t\t"+"--chrcoord  [value]"+"\t\t\t"     +"Print the chr and coordinate, useful for tabix indexing                   (default= "+booleanAsString(printChrCoord)+")"+"\n"+
	"\t\t"+"--minind    [value]"+"\t\t\t"     +"Require this minimum number of individuals to have the alternative allele (default= "+stringify(minIndWithAlt)+")"+"\n"+
	"\t\t"+""                   +"\t\t\t\t\t" +"At 1, it will allow mutations private to an individual "+"\n"+
	"\t\t"+"--failpop   [value]"+"\t\t\t"     +"Allow this many # of populations to fail the quality test, must still be in the VCF file (default= "+stringify(allowedFailPop)+")"+"\n"+
	"\t\t"+"--nonseg             "+"\t\t\t"   +"Print non-segrating sites as well                                         (default= "+booleanAsString(!onlySegSites)+")"+"\n";

    
 

    string refvcffiles=string("")+
	"\t\t"+"--refpre [path]" +"\t\t\t\t"+ "Prefix for reference vcf files           (default: none)\n"+
	"\t\t"+"--refsuf [suf]"  +"\t\t\t\t"+ "Suffix for reference vcf files           (default: none)\n";
    string refOptions=string("")+
	"\t\t"+"--refhum"        +"\t\t\t\t"+"Use human genome as the reference allele (default: "+booleanAsString(useHumanGenomeAsReference)+")\n";
    string ancOptions=string("")+
	"\t\t"+"--useanc"        +"\t\t\t\t"+"Use inferred ancestral, not chimp, from EPO as ancestral allele    (default: "+booleanAsString(useEPOInferedAsAncestor)+")\n";
    string dstatsOptions=string("")+
	"\t\t"+"--nohomocond"    +"\t\t\t\t"+"Do not require the conditional to be homozygous, sample at random  (default: "+booleanAsString(requiredHomozygousCond)+")\n";





    string smpvcffiles=string("")+
	"\t\t"+"--smppre [path]" +"\t\t\t\t"+ "Prefix for sample vcf files              (default: none)\n"+
	"\t\t"+"--smpsuf [suf]" +"\t\t\t\t"+  "Suffix for sample vcf files              (default: none)\n";

   string condvcf=string("")+
	"\t\t"+"--condpre [path]" +"\t\t\t"+  "Prefix for conditional vcf files         (default: none)\n"+
	"\t\t"+"--condsuf [suf]"  +"\t\t\t\t"+"Suffix for conditional vcf files         (default: none)\n";


    // string refcutoff=string("")+
    // 	"\t\t\t"+"--minCovREF [comma separated cov cutoffs]" +"\t\t\t"+"Minimal coverage reference (default: "+stringify(minCovcutoffREF)+")\n"+
    // 	"\t\t\t"+"--maxCovREF [comma separated cov cutoffs]" +"\t\t\t"+"Maximal coverage reference (default: "+stringify(maxCovcutoffREF)+")\n";

    // string smcutoff=string("")+
    // 	"\t\t\t"+"--minCovSMP [comma separated cov cutoffs]" +"\t\t\t"+"Minimal coverage sample (default: "+stringify(minCovcutoffSMP)+")\n"+
    // 	"\t\t\t"+"--maxCovSMP [comma separated cov cutoffs]" +"\t\t\t"+"Maximal coverage sample (default: "+stringify(maxCovcutoffSMP)+")\n";

    // string condcutoff=string("")+
    // 	"\t\t\t"+"--minCovCOND [comma separated cov cutoffs]" +"\t\t\t"+"Minimal coverage conditional (default: "+stringify(minCovcutoffCOND)+")\n"+
    // 	"\t\t\t"+"--maxCovCOND [comma separated cov cutoffs]" +"\t\t\t"+"Maximal coverage conditional (default: "+stringify(maxCovcutoffCOND)+")\n";


    string refcutoff=string("")+
    	"\t\t\t"+"--minCovREF " +"\t\t\t"+"Minimal coverage reference (default: "+stringify(minCovcutoffREF)+")\n"+
    	"\t\t\t"+"--maxCovREF " +"\t\t\t"+"Maximal coverage reference (default: "+stringify(maxCovcutoffREF)+")\n";

    string smcutoff=string("")+
    	"\t\t\t"+"--minCovSMP " +"\t\t\t"+"Minimal coverage sample (default: "+stringify(minCovcutoffSMP)+")\n"+
    	"\t\t\t"+"--maxCovSMP " +"\t\t\t"+"Maximal coverage sample (default: "+stringify(maxCovcutoffSMP)+")\n";

    string condcutoff=string("")+
    	"\t\t\t"+"--minCovCOND " +"\t\t\t"+"Minimal coverage conditional (default: "+stringify(minCovcutoffCOND)+")\n"+
    	"\t\t\t"+"--maxCovCOND " +"\t\t\t"+"Maximal coverage conditional (default: "+stringify(maxCovcutoffCOND)+")\n";


    string sizeRegions=string("")+
	"\tOnly consider genomic loci in these regions:\n"+
	"\t\t"+"--bed    [file]"     +"\t\t\t\t"+"BED file with the sorted non-overlapping genomic intervals to consider\n"+
	"\t\t"+"--bedtx  [file]"     +"\t\t\t\t"+"Genome-wide tabix indexed file of sorted non-overlapping genomic intervals to consider\n"+
	"\n"+
	"\tLoci size option:\n"+
	"\t\t"+"--chunk   [size]"     +"\t\t\t"+"Size of contiguous genomic region to take                         (default: none)\n"+
	"\t\t"+"--overlap [size]"     +"\t\t\t"  +"Size of overlap between windows (for overlap mode)                (default: "+stringify(overlapBetweenWindows)+")\n"+
	"\t\t"+"--frac    [fraction]" +"\t\t\t"  +"Fraction of sites passing quality filters e.g  --frac 0.3 for 30% (default: none)\n";
	
    string qualoptions=string("")+
	"\t\t"+"--minGQ [gq]"       +"\t\t\t\t" +"Minimal genotype quality                                         (default: "+stringify(minGQcutoff)+")\n"+
	"\t\t"+"--minMQ [mq]"       +"\t\t\t\t" +"Minimal mapping quality                                          (default: "+stringify(minMQcutoff)+")\n"+
	"\t\t"+"--minPL [pl]"       +"\t\t\t\t" +"Use this as the minimum difference of PL values for alleles      (default: "+stringify(minPLdiffind)+")\n"+
	"\t\t"+"--minMap20 [map20]" +"\t\t\t"   +"Mapability score                                                 (default: "+stringify(minMapabilitycutoff)+")\n"+
	"\t\t"+"--allowsexchr"      +"\t\t\t\t" +"Allow sites to be generated on X and Y                           (default: "+booleanAsString(allowSexChr)+")\n"+
	"\t\t"+"--allowindel"       +"\t\t\t\t" +"Allow sites considered within 5bp of an indel                    (default: "+booleanAsString(allowCloseIndelProx)+")\n"+
	"\t\t"+"--allowrm"          +"\t\t\t\t" +"Allow sites labeled repeat masked                                (default: "+booleanAsString(allowRepeatMasking)+")\n"+
	"\t\t"+"--allowSysErr"      +"\t\t\t\t" +"Allow sites labeled as systematic error                          (default: "+booleanAsString(allowSysErr)+")\n"+
	"\n\t\t"+"--allowall"       +"\t\t\t\t" +"Allow all sites, good for pre-filtered VCF                       (default: "+booleanAsString(allowall)+")\n"+
	"\n\t\t"+"--allowallMQ"     +"\t\t\t\t" +"Allow all sites, but still filter on MQ                          (default: "+booleanAsString(allowallMQ)+")\n";
    // bool allowCloseIndelProx = false;
    // bool allowRepeatMasking  = false;
    // bool allowSysErr         = false;



    string cpgoptions=string("")+
	"\t\t"+"--cpg "            +"\t\t\t\t\t"+"Consider CpGs in the reference or sample as CpGs                 (default: "+booleanAsString(cpgData)+")\n"+
	                        "\t\t\t\t\t\t\t"+"By default, the program only consider CpGs in the EPO flags\n";

    string maxdivoptions=string("")+
	"\t\t"+"--maxdiv "            +"\t\t\t\t"+"Maximize divergence at het. sites              (default: "+booleanAsString(maximizeDiv)+")\n\n";

    string externalDatabase=string("")+
	"\t\t"+"--fai [file]"      +"\t\t\t\t"+"Fasta index for genome (produced by \"samtools faidx\")          (default: none)\n"+
	"\t\t"+"--epo [directory]" +"\t\t\t"  +"Directory of the parsed EPO alignment (produces by EPOparser)  (default: none)\n";




    
    if(argc == 1 ||
       (argc == 2 && (string(argv[1]) == "-h" || string(argv[1]) == "--help") )
       ){
	cerr << "Usage "<<usage<<endl;
	return 1;       
    }


    //DEFINING OPTIONS FOR EACH MODE
    if( string(argv[1]) == "divcompute"){
	mode='d';	
	usage=string(string(argv[0])+" "+string(argv[1])+" [options]\n\n"+
		     selection+"\n"+
		     sizeRegions+
		     "\n\tLocation of database files:\n"+		     
		     externalDatabase+
		     "\n\tLocation of input files:\n"+		     
		     refvcffiles+
		     refOptions+
		     smpvcffiles+
		     "\n\tQuality filters:\n"+
		     "\t\tReference:\n"+
		     refcutoff+
		     "\t\tSample:\n"+
		     smcutoff+
		     "\t\tVCF filters:\n"+
		     qualoptions+
		     "\t\tMisc:\n"+
		     maxdivoptions+
		     cpgoptions);

    }else{
	if( string(argv[1]) == "ancall" ){
	    mode='a';
	    usage=string(string(argv[0])+" "+string(argv[1])+" [options]\n\n"+
			 selection+"\n"+
			 sizeRegions+
			 "\n\tLocation of database files:\n"+		     
			 externalDatabase+
			 "\n\tLocation of input files:\n"+		     
			 refvcffiles+
			 "\n\tQuality filters:\n"+
			 "\t\tReference:\n"+
			 refcutoff+
			 "\t\tGeneral:\n"+
			 qualoptions);	  
	}else{
	    if(string(argv[1]) == "hetero" ){
		mode='h';
		usage=string(string(argv[0])+" "+string(argv[1])+" [options]\n\n"+
			     selection+"\n"+
			     sizeRegions+
			     "\n\tLocation of database files:\n"+		     
			     externalDatabase+
			     "\n\tLocation of input files:\n"+		     
			     refvcffiles+
			     "\n\tQuality filters:\n"+
			     "\t\tReference:\n"+
			     "\t"+refcutoff+
			     "\t\tGeneral:\n"+
			     qualoptions);	  
	    }else{
		if(string(argv[1]) == "multivcf" ){
		    mode='t';
		    usage=string(string(argv[0])+" "+string(argv[1])+" [options]\n\n"+
				 selection+"\n"+
				 sizeRegions+"\n"+			 
				 multivcf+
				 ancOptions+
				 "\n\tLocation of database files:\n"+		     
				 externalDatabase+
				 "\n\tQuality filters:\n"+
				 "\t\tGeneral:\n"+
				 qualoptions);	  		   
		}else{
		    // if( string(argv[1]) == "dstatcond"){
		    // 	mode='c';	
		    // 	usage=string(string(argv[0])+" "+string(argv[1])+" [options]\n\n"+
		    // 		     selection+"\n"+
		    // 		     sizeRegions+
		    // 		     "\n\tLocation of database files:\n"+		     
		    // 		     externalDatabase+
		    // 		     "\n\tLocation of input files:\n"+		     
		    // 		     refvcffiles+
		    // 		     refOptions+
		    // 		     smpvcffiles+
		    // 		     condvcf+
		    // 		     "\n\tQuality filters:\n"+
		    // 		     "\t\tReference:\n"+
		    // 		     refcutoff+
		    // 		     "\t\tSample:\n"+
		    // 		     smcutoff+
		    // 		     "\t\tCondition:\n"+
		    // 		     condcutoff+
		    // 		     "\t\tGeneral:\n"+
		    // 		     ancOptions+
		    // 		     qualoptions+
		    // 		     cpgoptions);	
		    // }else{
		    if( string(argv[1]) == "dstat"){
			mode='s';	
			usage=string(string(argv[0])+" "+string(argv[1])+" [options]\n\n"+
				     selection+"\n"+
				     sizeRegions+
				     "\n\tLocation of database files:\n"+		     
				     externalDatabase+
				     "\n\tLocation of input files:\n"+		     
				     refvcffiles+
				     refOptions+
				     smpvcffiles+
				     condvcf+
				     "\n\tQuality filters:\n"+
				     "\t\tReference:\n"+
				     refcutoff+
				     "\t\tSample:\n"+
				     smcutoff+
				     "\t\tCondition:\n"+
				     condcutoff+
				     "\t\tGeneral:\n"+
				     ancOptions+
				     dstatsOptions+
				     qualoptions+
				     cpgoptions);	
		
		    }else{
			cerr << "Invalid mode\nUsage:\n"<<usage<<endl;
			return 1; 
		    }
	         
		}
	    }
	}
    }
    
    ////////////////////////////////////
    //                                //
    //      END DEFINE OPTIONS        //
    //                                //
    ////////////////////////////////////
       









































    ///////////////////////////////////////////
    //                                       //
    //      BEGIN PARSING CMD LINE OPTIONS   //
    //                                       //
    ///////////////////////////////////////////
       




    if(argc == 2 ||
       (argc == 3 && (string(argv[2]) == "-h" || string(argv[2]) == "--help") )
       ){
	cerr << "Usage "<<usage<<endl;
	return 1;       
    }


    bool specifiedPL  =false;
    bool specifiedGQ=false;

#ifdef DEBUG
    cerr<<"Parsing options"<<endl;
#endif


    //starts at two for the sub-program name
    for(int i=2;i<(argc);i++){ 

	////////////////////
	//  COMMON TO ALL //
	////////////////////

	//BEGIN loci selection
        if(strcmp(argv[i],"--bed") == 0 ){
	     bedFileRegions    = string(argv[i+1]);
	     useBedFileRegions = true;  
	     i++;
	     continue;
	}


	// bool   useBedFileRegions     =false;
	// bool   useBedFileRegionsTabix=false;

        if(strcmp(argv[i],"--bedtx") == 0 ){
	     bedFileRegionsTabix    = string(argv[i+1]);
	     useBedFileRegionsTabix = true;  
	     i++;
	     continue;
	}


	
        if(strcmp(argv[i],"--random") == 0 ){
	    if(lociSelection != 'x'){
		cerr << "Choose only one loci selection option"<<endl;
		return 1;       
	    }
	    amountOfGoodSitesTARGET=destringify<int>(argv[i+1]);
	    lociSelection='r';
	    i++;
            continue;
        }

        if(strcmp(argv[i],"--loci") == 0 ){
	    if(lociSelection != 'x'){
		cerr << "Choose only one loci selection option"<<endl;
		return 1;       
	    }

	    lociSelection='o';
            continue;
        }

        if(strcmp(argv[i],"--locichr") == 0 ){
	    if(lociSelection != 'x'){
		cerr << "Choose only one loci selection option"<<endl;
		return 1;       
	    }
	    chrToUse=string(argv[i+1]);
	    i++;	    
	    lociSelection='c';
            continue;
        }

        if(strcmp(argv[i],"--region") == 0 ){
	    if(lociSelection != 'x'){
		cerr << "Choose only one loci selection option"<<endl;
		return 1;       
	    }
	    bedFileWithRegions=string(argv[i+1]);
	    i++;	    
	    lociSelection='f';
            continue;
        }


        if(strcmp(argv[i],"--chr") == 0 ){
	    if(lociSelection != 'x'){
		cerr << "Choose only one loci selection option"<<endl;
		return 1;       
	    }
	    chrToUse=string(argv[i+1]);
	    i++;	    
	    lociSelection='v'; 
            continue;
        }

	if(strcmp(argv[i],"--wide") == 0 ){
	    if(lociSelection != 'x'){
		cerr << "Choose only one loci selection option"<<endl;
		return 1;       
	    }
	    lociSelection='w';
            continue;
        }
	//END loci selection


	//loci size
        if(strcmp(argv[i],"--chunk") == 0 ){
	    bpToExtract=destringify<int>(argv[i+1]);
	    i++;
            continue;
        }

	if(strcmp(argv[i],"--overlap") == 0 ){
	    overlapBetweenWindows=destringify<int>(argv[i+1]);
	    i++;
            continue;
        }

	if(strcmp(argv[i],"--frac") == 0 ){
	    fractionGood=destringify<double>(argv[i+1]);
	    i++;
            continue;
        }





	//location of database
        if(strcmp(argv[i],"--fai") == 0 ){
	    fastaIndex=string(argv[i+1]);
	    i++;
            continue;
        }

        if(strcmp(argv[i],"--epo") == 0 ){
	    epoIndex=string(argv[i+1]);
	    i++;
            continue;
        }



	//Quality filters

	if(strcmp(argv[i],"--minPL") == 0 ){
	    minPLdiffind=destringify<int>(argv[i+1]);
	    specifiedPL  =true;
	    // minGQcutoff=0;
	    // cerr<<"vcfcompute: Warning: using minPL sets the GQ cutoffs to 0"<<endl;
	    i++;
	    continue;
	}


        if(strcmp(argv[i],"--minGQ") == 0 ){
	    minGQcutoff=destringify<int>(argv[i+1]);
	    specifiedGQ=true;    
	    i++;
            continue;
        }

	if(strcmp(argv[i],"--minMQ") == 0 ){
	    minMQcutoff=destringify<int>(argv[i+1]);
	    i++;
            continue;
        }

	if(strcmp(argv[i],"--minMap20") == 0 ){
	    minMapabilitycutoff=destringify<double>(argv[i+1]);
	    i++;
	    continue;
	}
	 
	if(strcmp(argv[i],"--allowsexchr") == 0 ){
	    allowSexChr=true;
	    continue;
	}

	if(strcmp(argv[i],"--allowindel") == 0 ){
	    allowCloseIndelProx =true;
	    continue;
	}

	if(strcmp(argv[i],"--allowrm") == 0 ){
	    allowRepeatMasking   =true;
	    continue;
	}

	if(strcmp(argv[i],"--allowSysErr") == 0 ){
	    allowSysErr     =true;
	    continue;
	}

	if(strcmp(argv[i],"--allowall") == 0 ){
	    allowall     =true;
	    continue;
	}

	if(strcmp(argv[i],"--allowallMQ") == 0 ){
	    allowallMQ    =true;
	    continue;
	}


	if(mode=='t'){ //options unique to treemix
	    


	    if(strcmp(argv[i],"--nonseg") == 0 ){
		onlySegSites=false;    
		continue;
	    }

	    if(strcmp(argv[i],"--failpop") == 0 ){
		allowedFailPop=destringify<int>(argv[i+1]);
		i++;
		continue;
	    }


	    if(strcmp(argv[i],"--outtype") == 0 ){
		outtype=destringify<int>(argv[i+1]);
		i++;
		continue;
	    }

	    if(strcmp(argv[i],"--popfile") == 0 ){
		popfile=string(argv[i+1]);
		i++;
		continue;
	    }

	    if(strcmp(argv[i],"--nochimp") == 0 ){
		printChimp=false;
		continue;
	    }



	    // if(strcmp(argv[i],"--pldiffind") == 0 ){
	    // 	minPLdiffind=destringify<int>(argv[i+1]);
	    // 	i++;
	    // 	continue;
	    // }

	    if(strcmp(argv[i],"--chrcoord") == 0 ){
		printChrCoord=true;
		continue;
	    }


	    if(strcmp(argv[i],"--minind") == 0 ){
		minIndWithAlt=destringify<int>(argv[i+1]);
		i++;
		continue;
	    }


	    ///////////////////////
	    // END COMMON TO ALL //
	    ///////////////////////
	    
	}else{ //no treemix/mistar

	    //Reference filters
	    if(strcmp(argv[i],"--minCovREF") == 0 ){
		minCovcutoffREF=destringify<int>(argv[i+1]);
		//setVectorUsedAsString( string(argv[i+1]), &minCovcutoffREF);
		// vector<string> tempStr=allTokens(string(argv[i+1]),",");
		// minCovcutoffREF.clear();
		// for(int ind=0;ind<tempStr.size();ind++)
		//     minCovcutoffREF.push_back( destringify<int>(tempStr[ind]) );
		i++;
		continue;
	    }

       
	    if(strcmp(argv[i],"--maxCovREF") == 0 ){
		maxCovcutoffREF=destringify<int>(argv[i+1]);
		//maxCovcutoffREF=allTokens(string(argv[i+1]),",");
		//setVectorUsedAsString( string(argv[i+1]), &maxCovcutoffREF);
		i++;
		continue;
	    }

	    if(strcmp(argv[i],"--refpre") == 0 ){
		refereVCFPre=string(argv[i+1]);
		i++;
		continue;
	    }

	    if(strcmp(argv[i],"--refsuf") == 0 ){
		refereVCFSuf=string(argv[i+1]);
		i++;
		continue;
	    }

	}//end no treemix/mistar

	if( mode=='s'){ //options unique to  d-stats
	    
	    if(strcmp(argv[i],"--condpre") == 0 ){
		condVCFPre=string(argv[i+1]);
		i++;
		continue;
	    }

	    if(strcmp(argv[i],"--condsuf") == 0 ){
		condVCFSuf=string(argv[i+1]);
		i++;
		continue;
	    }

	

	    if(strcmp(argv[i],"--nohomocond") == 0 ){
		requiredHomozygousCond=false;
		continue;
	    }


	    if(strcmp(argv[i],"--minCovCOND") == 0 ){
		minCovcutoffCOND=destringify<int>(argv[i+1]);
		//setVectorUsedAsString( string(argv[i+1]), &minCovcutoffCOND);
		i++;
		continue;
	    }

       
	    if(strcmp(argv[i],"--maxCovCOND") == 0 ){
		maxCovcutoffCOND=destringify<int>(argv[i+1]);
		//setVectorUsedAsString( string(argv[i+1]), &maxCovcutoffCOND);
		i++;
		continue;
	    }

	} //end options unique to  d-stats

	if( mode=='d' ){ //options unique to divergence or  d-stats

	    if(strcmp(argv[i],"--maxdiv") == 0 ){
		maximizeDiv=true;
		continue;
	    }

	}

	if( mode=='s' || mode=='d' ){ //options unique to divergence or  d-stats


	    if(strcmp(argv[i],"--refhum") == 0 ){
		useHumanGenomeAsReference=true;
		continue;
	    }

	    if(strcmp(argv[i],"--cpg") == 0 ){
		cpgData=true;
		continue;
	    }

	    if(strcmp(argv[i],"--minCovSMP") == 0 ){
		minCovcutoffSMP=destringify<int>(argv[i+1]);
		//setVectorUsedAsString( string(argv[i+1]), &minCovcutoffSMP);
		i++;
		continue;
	    }

       
	    if(strcmp(argv[i],"--maxCovSMP") == 0 ){
		maxCovcutoffSMP=destringify<int>(argv[i+1]);
		//setVectorUsedAsString( string(argv[i+1]), &maxCovcutoffSMP);
		i++;
		continue;
	    }

	    if(strcmp(argv[i],"--smppre") == 0 ){
		sampleVCFPre=string(argv[i+1]);
		i++;
		continue;
	    }

	    if(strcmp(argv[i],"--smpsuf") == 0 ){
		sampleVCFSuf=string(argv[i+1]);
		i++;
		continue;
	    }
	}// end options unique to divergence or  d-stats

	if( mode=='t' || mode=='d' ){ //options unique to multivcf or  d-stats

	    if(strcmp(argv[i],"--useanc") == 0 ){
		useEPOInferedAsAncestor=true;
		continue;
	    }

	}

        cerr<<"Unknown option "<<argv[i] <<" exiting"<<endl;
        return 1;           

    }//end argc


    //Double checking options
    // case 'x':  //undefined

    // case 'r':  //random windows
    // case 'o':  //entire genome per window
    // case 'c':  //entire chr per window	    

    // case 'f':  //reading regions from a file
    // case 'w':  //genome wide 
    // case 'v':  //entire chr (only one window)


    if(useBedFileRegions && useBedFileRegionsTabix){
	cerr << "Cannot specify a set of regions if not using genome wide or entire chr"<<endl;
	return 1;   	    	    
    }

    if(useBedFileRegions){

	//if it's anything but 'w' and 'v', error
	if(!(lociSelection == 'w' || //genome wide
	     lociSelection == 'v' )){ //entire chr
	     cerr << "Cannot specify a set of regions if not using genome wide or entire chr"<<endl;
	     return 1;   	    
	}
	cerr<<"Reading bed file ...";
	bedRegionsToFilter = readBEDSortedfile(bedFileRegions);
	cerr<<"..done"<<endl;
    }

    if(useBedFileRegionsTabix){
	if(isFile(bedFileRegionsTabix+".tbi")){
	    cerr << "The file "<<bedFileRegionsTabix<<" should have a tabix index"<<endl;
	    return 1;   	    	    	
	}
    }



     if(lociSelection == 'x'){
	 cerr << "Choose a loci selection option"<<endl;
	 return 1;       
     }

     if(lociSelection == 'f' || //file
	lociSelection == 'r' || //random windows
	lociSelection == 'v' || //entire chr
	lociSelection == 'w'    //genome wide
	){
	 //if it's not a specific chr or genome wide, the user must specify the # of bp to extract and the fraction passing qual filters
	 if(overlapBetweenWindows != 0){
	     cerr << "Cannot specify an overlap for this loci selection mode"<<endl;
	     return 1;   
	 }
     }


     //if it's not a specific chr or genome wide, the user must specify the # of bp to extract and the fraction passing qual filters
     if(lociSelection != 'v' &&  //entire chr
	lociSelection != 'w' &&	 //genome wide
	lociSelection != 'f' 	 //reading a file
	){ 
	 if(bpToExtract == -1 ){
	     cerr << "Need to specify the size of the size of contiguous genomic region to take"<<endl;
	     return 1;       
	 }
     
	 if(fractionGood == -1.0){
	     cerr << "Need to specify the fraction of sites passing quality filters"<<endl;
	     return 1;       
	 }

	 minimumNumberOfGoodSites=int( fractionGood*bpToExtract ); //amount of bp passing filters that are needed per chunk
     }else{
	 minimumNumberOfGoodSites=0;

	 if(bpToExtract != -1 ){
	     cerr << "Error: no need to specify the size of the size of contiguous genomic region to take"<<endl;
	     return 1;       
	 }
     
	 if(fractionGood != -1.0){
	     cerr << "Error: no need to specify the fraction of sites passing quality filters"<<endl;
	     return 1;       
	 }


     }

    
    if(epoIndex.empty()){
	cerr << "Need to specify the directory of the parsed EPO alignment"<<endl;
	return 1;       
    }

    if(fastaIndex.empty()){
	cerr << "Need to specify the directory of the fasta index for genome"<<endl;
	return 1;       
    }


    if(mode != 't'){ //no need for reference for multivcf

	if(!useHumanGenomeAsReference){
	    if(refereVCFPre.empty()){
		cerr << "Need to specify the prefix for reference vcf files"<<endl;
		return 1;       
	    }
	

	    if(refereVCFSuf.empty()){
		cerr << "Need to specify the suffix for reference vcf files"<<endl;
		return 1;       
	    }
	}

    }else{

	if( popfile == ""){
	    cerr << "Need to specify the population file with --popfile"<<endl;
	    return 1;       
	}    
	if( outtype == -1){
	    cerr << "Need to specify the output type with --outtype"<<endl;
	    return 1;       
	}    

    }

    if(mode=='d'){ //unique to divergence 

	if(useHumanGenomeAsReference){
	    if( minCovcutoffREF!=0){
		cerr << "No need to specify the min coverage for the human reference "<<endl;
		return 1;  
	    }
	    if( maxCovcutoffREF!=1000){
		cerr << "No need to specify the max coverage for the human reference "<<endl;
		return 1;  
	    }
	    if(!refereVCFPre.empty()){
		cerr << "No need to specify the prefix for the human reference"<<endl;
		return 1;       
	    }
	    

	    if(!refereVCFSuf.empty()){
		cerr << "No need to specify the suffix for the human reference"<<endl;
		return 1;       
	    }
	}

	if(sampleVCFPre.empty()){
	    cerr << "Need to specify the prefix for sample vcf files"<<endl;
	    return 1;       
	}
	

	if(sampleVCFSuf.empty()){
	    cerr << "Need to specify the suffix for sample vcf files"<<endl;
	    return 1;       
	}

    }else{
	if(useHumanGenomeAsReference){
	    cerr << "Cannot use the human reference for this mode"<<endl;
	    return 1;       
	}
    }

    if( specifiedPL  &&   specifiedGQ ){
	cerr << "Cannot specify both PL filtering and GQ filtering"<<endl;
	return 1;       
    }



    //////////////////////////////////////////
    //
    //      END PARSING CMD LINE OPTIONS
    //
    /////////////////////////////////////////
       







    







    //BEGIN loci selection

    if(mode == 't'){ 
	//tmim = new MultiVCFParser (popfile,printChimp,minPLdiffind,printChrCoord,minIndWithAlt,allowedFailPop) ;
	tmim = new MultiVCFParser (popfile,printChimp,minPLdiffind,printChrCoord,minIndWithAlt,allowedFailPop,
				   minGQcutoff          ,
				   minMQcutoff          ,
				   minMapabilitycutoff  ,
				   !allowCloseIndelProx ,
				   !allowRepeatMasking  ,
				   !allowSysErr         ,
				   allowall,
				   allowallMQ) ;
    }else{
	tmim=0;
    }



    RandomGenomicCoord rgc (fastaIndex,allowSexChr);
    GenomicWindows     rw  (fastaIndex,allowSexChr);
    GenomicRange       rangeFound;
    vector<GenomicRange> v;

    //BEGIN generate the random coord
    unsigned int indexOfLoci=0;
    bool loopingCondition;



    switch(lociSelection){
    case 'r': //random loci
	loopingCondition=(indexOfLoci<amountOfGoodSitesTARGET);
	rangeFound=rgc.getRandomGenomic(bpToExtract);
	break;
    case 'c':  //windows on entire chr
	v=rw.getGenomicWindowsChr(chrToUse,bpToExtract,overlapBetweenWindows);
	if( v.size() == 0 ){    cerr<<"No range found using these chr/loci settings"<<endl; return 1;}
	loopingCondition=(indexOfLoci<v.size());   
	rangeFound=v[indexOfLoci];
	break;
    case 'f':    //read a bed file
	v=readBEDfile(bedFileWithRegions);	//a bed file provides the regions to use
	if( v.size() == 0 ){    cerr<<"No range found using these chr/loci settings"<<endl; return 1;}
	loopingCondition=(indexOfLoci<v.size());   
	rangeFound=v[indexOfLoci];
	break;	
    case 'o':    //windows on entire genome
	v=rw.getGenomicWindows(bpToExtract,overlapBetweenWindows);
	if( v.size() == 0 ){    cerr<<"No range found using these chr/loci settings"<<endl; return 1;}
	loopingCondition=(indexOfLoci<v.size());   
	rangeFound=v[indexOfLoci];
	break;
    case 'v':  //entire chr
	v=rw.getChr(chrToUse);
	if( v.size() == 0 ){    cerr<<"No range found using these chr/loci settings"<<endl; return 1;}
	loopingCondition=(indexOfLoci<v.size());   
	rangeFound=v[indexOfLoci];
	break;
    case 'w':  //entire genome
	v=rw.getGenomeWide();
	if( v.size() == 0 ){    cerr<<"No range found using these chr/loci settings"<<endl; return 1;}
	loopingCondition=(indexOfLoci<v.size());   
	rangeFound=v[indexOfLoci];
	break;


    default:
	cerr << "Unknown loci selection"<<endl;
	return 1;       
	break;
    }


    if(useBedFileRegionsTabix){
	bedFileTabix=new ReadTabix(bedFileRegionsTabix,
				   bedFileRegionsTabix+".tbi",
				   rangeFound.getChrName(),
				   rangeFound.getStartCoord(),
				   rangeFound.getEndCoord());
    }
    //END loci selection

#ifdef DEBUG
    cerr<<"done loci selection"<<endl;
#endif




    //we negate  allowCloseIndelProx,allowRepeatMasking,allowSysErr because
    //allow = false means filter =true in setVCF
    DivergenceResult justToPrint;
    switch(mode){

    case 'd': //divergence
	
	filtersToUse->push_back( new SetVCFFilters(minGQcutoff          ,
						   minMQcutoff          ,
						   minMapabilitycutoff  ,
						   !allowCloseIndelProx ,
						   !allowRepeatMasking  ,
						   !allowSysErr         ,
						   minCovcutoffREF      ,
						   maxCovcutoffREF      ,
						   allowall,
						   allowallMQ) );
	(*filtersToUse)[0]->setName("Reference");

	filtersToUse->push_back(new SetVCFFilters(minGQcutoff          ,
						  minMQcutoff          ,
						  minMapabilitycutoff  ,
						  !allowCloseIndelProx ,
						  !allowRepeatMasking  ,
						  !allowSysErr         ,
						  minCovcutoffSMP      ,
						  maxCovcutoffSMP      ,
						  allowall,
						  allowallMQ) );
	(*filtersToUse)[1]->setName("Sample");
	// //printing header

	cout<<"#loc\t"<<justToPrint.getHeader()<<endl;

	break;

    case 's': // d-stats

	filtersToUse->push_back(new SetVCFFilters(minGQcutoff          ,
						  minMQcutoff          ,
						  minMapabilitycutoff  ,
						  !allowCloseIndelProx ,
						  !allowRepeatMasking  ,
						  !allowSysErr         ,
						  minCovcutoffREF      ,
						  maxCovcutoffREF      ,
						  allowall,
						  allowallMQ ) );
	(*filtersToUse)[0]->setName("Reference");

	filtersToUse->push_back(new SetVCFFilters(minGQcutoff          ,
						  minMQcutoff          ,
						  minMapabilitycutoff  ,
						  !allowCloseIndelProx ,
						  !allowRepeatMasking  ,
						  !allowSysErr         ,
						  minCovcutoffSMP      ,
						  maxCovcutoffSMP      ,
						  allowall,
						  allowallMQ) );
	(*filtersToUse)[1]->setName("Sample");

	filtersToUse->push_back(new SetVCFFilters(minGQcutoff          ,
						  minMQcutoff          ,
						  minMapabilitycutoff  ,
						  !allowCloseIndelProx ,
						  !allowRepeatMasking  ,
						  !allowSysErr         ,
						  minCovcutoffCOND     ,
						  maxCovcutoffCOND     ,
						  allowall,
						  allowallMQ) );
	(*filtersToUse)[2]->setName("Conditional");

	break;

    case 'a': //ancestral stretches
	filtersToUse->push_back(new SetVCFFilters(minGQcutoff          ,
						  minMQcutoff          ,
						  minMapabilitycutoff  ,
						  !allowCloseIndelProx ,
						  !allowRepeatMasking  ,
						  !allowSysErr         ,
						  minCovcutoffREF      ,
						  maxCovcutoffREF      ,
						  allowall,
						  allowallMQ) );
	(*filtersToUse)[0]->setName("Reference");
	
	break;
    case 'h': //heterozygosity

	filtersToUse->push_back(new SetVCFFilters(minGQcutoff          ,
						  minMQcutoff          ,
						  minMapabilitycutoff  ,
						  !allowCloseIndelProx ,
						  !allowRepeatMasking  ,
						  !allowSysErr         ,
						  minCovcutoffREF      ,
						  maxCovcutoffREF      ,
						  allowall,
						  allowallMQ) );
	(*filtersToUse)[0]->setName("Reference");

	break;
    case 't': //multi vcf

	//should be handled internally
	
	break;

    default:
	cerr << "Unknown mode "<<mode<<endl;
	return 1;       
	break;	    

    }







    string programLine;
    for(int i=0;i<(argc);i++){ 
	programLine+=(string(argv[i])+" ");
    }
    string gitHubVersion=returnGitHubVersion(argv[0],"");

    //BEGIN summary of parameters used
    cerr<<"Parameters used:"<<endl;
    cerr<<"Date: "<<getDateString()<<"\t"<<getTimeString()<<endl;
    cerr<<"Program line : "<<programLine<<endl;
    cerr<<"Github version: "<<gitHubVersion<<endl;

    cerr<<"Quality filters:"<<endl;

    if(mode == 't'){
	tmim->printCutoffs();
    }else{
      for(unsigned int indexCut=0;
	  indexCut<filtersToUse->size();
	  indexCut++){
	  cerr<<(*(*filtersToUse)[indexCut])<<endl;
      }
    }
    cerr<<"PL difference\t\t"<<minPLdiffind<<endl;
    cerr<<"cpg in Ref/Sample\t\t"<<boolStringify(cpgData)<<endl;
    cerr<<"Do not apply any quality filter \t\t"<<boolStringify(allowall)<<endl;
    cerr<<"Do not apply any quality filter just MQ\t\t"<<boolStringify(allowallMQ)<<endl;

    cerr<<"bpForIndels\t\t"<<    bpForIndels<<endl;
    cerr<<"Loci options:"<<endl;

    cerr<<"mode\t\t\t"     <<    mode<<endl;
    cerr<<"amountOfGoodSitesTARGET\t\t"<<amountOfGoodSitesTARGET<<endl;
    cerr<<"allowSexChr\t\t"<<    allowSexChr<<endl;
    cerr<<"overlapBetweenWindows\t\t"<<    overlapBetweenWindows<<endl;
    cerr<<"lociSelection\t\t"<<    lociSelection<<endl;


    cerr<<"Misc options:"<<endl;

    cerr<<"epoIndex\t\t"<<    epoIndex<<endl;
    cerr<<"fastaIndex\t\t"<<    fastaIndex<<endl;
    cerr<<"useBedFileRegions\t"     <<boolStringify(useBedFileRegions)<<endl;
    cerr<<"useBedFileRegionsTabix\t"<<boolStringify(useBedFileRegionsTabix)<<endl;
    cerr<<"bedFileRegions\t"     <<bedFileRegions<<endl;
    cerr<<"bedFileRegionsTabix\t"<<bedFileRegionsTabix<<endl;

    cerr<<"bpToExtract\t\t"<<    bpToExtract<<endl;
    cerr<<"fractionGood\t\t"<<    fractionGood<<endl;
    cerr<<"minimumNumberOfGoodSites\t\t"<<minimumNumberOfGoodSites<<endl;

    cerr<<"Input options:"<<endl;

    cerr<<"refereVCFPre\t\t"<<  (refereVCFPre) <<endl;
    cerr<<"refereVCFSuf\t\t"<<  (refereVCFSuf) <<endl;
    cerr<<"sampleVCFPre\t\t"<<  (sampleVCFPre) <<endl;
    cerr<<"sampleVCFSuf\t\t"<<  (sampleVCFSuf) <<endl;
    cerr<<"condVCFPre\t\t"  <<  (condVCFPre)   <<endl;
    cerr<<"condVCFSuf\t\t"  <<  (condVCFSuf)   <<endl;
    cerr<<"use human genome as ref allele\t\t"<<boolStringify(useHumanGenomeAsReference)<<endl;




    cerr<<"Treemix/mistar options:"<<endl;
    cerr<<"use EPO as ancestor\t\t"<<boolStringify(useEPOInferedAsAncestor)<<endl;
    cerr<<"print the chr/coord\t\t"<<boolStringify(printChrCoord)<<endl;
    cerr<<"print the chimp\t\t"<<boolStringify(printChimp)<<endl;
    cerr<<"Population file\t\t"<<popfile<<endl;
    cerr<<"Minimum number of inds for mutation\t\t"<<minIndWithAlt<<endl;
    cerr<<"Allow this # of populations to fail the quality checks\t\t"<<allowedFailPop<<endl;
    cerr<<"Only print segregating sites\t\t"<<boolStringify(onlySegSites)<<endl;

    cerr<<"D-stats:"<<endl;
    cerr<<"Required homozygosity on conditional\t\t"<<boolStringify(requiredHomozygousCond)<<endl;
    cerr<<endl<<endl;

    //END summary of parameters used
    
    


    do{ //closed by while(loopingCondition)
	
	//cout<<"rangeFound "<<rangeFound<<endl;
    	// string refereVCF   =refereVCFPre+rangeFound.getChrName()+refereVCFSuf;
    	// string refereVCFidx=refereVCF+".tbi";

    	string refereVCF     = "";
    	string refereVCFidx  = "";

    	if(!useHumanGenomeAsReference){
    	    refereVCF    = refereVCFPre+rangeFound.getChrName()+refereVCFSuf;
    	    refereVCFidx = refereVCF+".tbi";
    	}	    
    

    	string sampleVCF   =sampleVCFPre+rangeFound.getChrName()+sampleVCFSuf;
    	string sampleVCFidx=sampleVCF+".tbi";


    	string condVCF   =condVCFPre+rangeFound.getChrName()+condVCFSuf;
    	string condVCFidx=condVCF+".tbi";


    	string epoFile     =epoIndex+"chr"+rangeFound.getChrName()+".epo.gz";
    	string epoFileidx  =epoFile+".tbi";

    	int returnCode=-1;

    	switch(mode){

    	case 'd': //divergence


	    //CODE FOR REGIONS, 
	    //this needs to be taken out of the switch
	    vector<GenomicRange> *  genomicRangesFilter;

	    if(useBedFileRegions){
		//no GenomicRange found, skip
		if( bedRegionsToFilter->find(rangeFound.getChrName() ) == bedRegionsToFilter->end() ){		    
		    genomicRangesFilter=0;
		    goto ENDOFIFCOND;		    
		}else{
		    genomicRangesFilter=(*bedRegionsToFilter)[ rangeFound.getChrName() ];		    
		}
	    }else{
		if(useBedFileRegionsTabix){
		       bedFileTabix->repositionIterator(rangeFound.getChrName(),
							rangeFound.getStartCoord(),
							rangeFound.getEndCoord());

		       genomicRangesFilter = new vector<GenomicRange> (); //if we use tabix, we actually allocate space
		       

		       string           buffer;
		       unsigned int     lastEndCoord = 0;

		       while(bedFileTabix->readLine(buffer)){

			   vector<string> temp=allTokens(buffer,'\t');
			   if(temp.size() != 3){
			       cerr << "Error in vcfcompute in readBEDSortedfile(): following line does not have 3 fields"<< buffer<<endl;
			       exit(1);		
			   }
			   string       chrName     = destringify<string>(temp[0]);
			   unsigned int startCoord  = destringify<unsigned int>(temp[1]);
			   unsigned int endCoord    = destringify<unsigned int>(temp[2]);

			   if(startCoord <= lastEndCoord ){
			       cerr << "Problem with line =  "<<buffer<<" the start of the coordinate lesser than the end of the previous record "<<lastEndCoord<<endl;
			       exit(1);
			   }
			   if(endCoord   <= lastEndCoord ){
			       cerr << "Problem with line =  "<<buffer<<" the end of the coordinate lesser than the end of the previous record "<<lastEndCoord<<endl;
			       exit(1);		    
			   }
	    
		       
			   GenomicRange toadd (chrName,startCoord,endCoord);	    
			   genomicRangesFilter->push_back(toadd);
		       }

		}else{
		    genomicRangesFilter=0;
		}
	    }
	    //END CODE FOR REGIONS
	    //returnCode=computeDivergence( refereVCF, refereVCFidx, sampleVCF,sampleVCFidx, useBedFileRegions,genomicRangesFilter , epoFile, epoFileidx, minimumNumberOfGoodSites, rangeFound, minGQcutoff, minCovcutoffREF, maxCovcutoffREF, minCovcutoffSMP, maxCovcutoffSMP, minMQcutoff, minMapabilitycutoff, bpForIndels,cpgData,useHumanGenomeAsReference,specifiedPL,minPLdiffind);

	    returnCode=computeDivergence( refereVCF, refereVCFidx, sampleVCF,sampleVCFidx, useBedFileRegions,genomicRangesFilter , epoFile, epoFileidx, minimumNumberOfGoodSites, rangeFound, 
					  (*filtersToUse)[0],
					  (*filtersToUse)[1],
					  //minGQcutoff, minCovcutoffREF, maxCovcutoffREF, minCovcutoffSMP, maxCovcutoffSMP, minMQcutoff, minMapabilitycutoff, 
					  bpForIndels,cpgData,useHumanGenomeAsReference,minPLdiffind,maximizeDiv);

	    //space was allocated for the genomicRangesFilter
	    if(useBedFileRegionsTabix){
		delete(genomicRangesFilter);
	    }
		       
	    
    	    break;


    	case 's': // d-stats
    	    returnCode=dstats(refereVCF, refereVCFidx, sampleVCF,sampleVCFidx,  condVCF, condVCFidx, epoFile, epoFileidx, minimumNumberOfGoodSites, rangeFound, 
			      (*filtersToUse)[0],
			      (*filtersToUse)[1],
			      (*filtersToUse)[2],
			      //minGQcutoff, minCovcutoffREF, maxCovcutoffREF, minCovcutoffSMP, maxCovcutoffSMP, minCovcutoffCOND, maxCovcutoffCOND, minMQcutoff, minMapabilitycutoff, 
			      bpForIndels,cpgData,useHumanGenomeAsReference,useEPOInferedAsAncestor,requiredHomozygousCond);
    	    break;

    	case 'a': //ancestral stretches
    	    returnCode=computeFracAnc( refereVCF, refereVCFidx,  epoFile, epoFileidx, minimumNumberOfGoodSites, rangeFound, (*filtersToUse)[0],  bpForIndels);
				       //minGQcutoff, minCovcutoffREF, maxCovcutoffREF, minMQcutoff, minMapabilitycutoff, bpForIndels);
    	    break;

    	case 'h': //heterozygosity
    	    returnCode=computeFracHetero( refereVCF, refereVCFidx,  minimumNumberOfGoodSites, rangeFound, (*filtersToUse)[0],  bpForIndels);
					  //minGQcutoff, minCovcutoffREF, maxCovcutoffREF, minMQcutoff, minMapabilitycutoff, bpForIndels);
    	    break;

    	case 't': //multi vcf
    	    //returnCode=computeDivergence( refereVCF, refereVCFidx, sampleVCF,sampleVCFidx, epoFile, epoFileidx, minimumNumberOfGoodSites, rangeFound, minGQcutoff, minCovcutoffREF, maxCovcutoffREF, minCovcutoffSMP, maxCovcutoffSMP, minMQcutoff, minMapabilitycutoff, bpForIndels,cpgData);
	    //returnCode =  tmim->produceOutput(epoFile,epoFileidx,rangeFound,minGQcutoff,minMQcutoff,minMapabilitycutoff,outtype,onlySegSites);
	    returnCode =  tmim->produceOutput(epoFile,epoFileidx,rangeFound,outtype,onlySegSites,useEPOInferedAsAncestor,programLine,gitHubVersion);
    	    break;


    	default:
    	    cerr << "Unknown mode"<<endl;
    	    return 1;       
    	    break;	    

    	}

    ENDOFIFCOND:


    	//increasing the index conditionally if random loci
    	if(lociSelection == 'r'){
    	    if(returnCode == 0 ){ //this loci worked
    		indexOfLoci++;
    	    }
    	}else{
    	    indexOfLoci++;
    	}
	    
    

    	    //for random windows, we ask for another ranom window
    	if(lociSelection == 'r'){
    	    loopingCondition=(indexOfLoci<amountOfGoodSitesTARGET);
	    if(loopingCondition)
		rangeFound=rgc.getRandomGenomic(bpToExtract);
	}else{

    	    //For the remaining case, we read from vector "v"
	    if(lociSelection == 'o' ||  //entire genome per window
	       lociSelection == 'c' ||  //entire chr per window	    
	       lociSelection == 'w' ||  //genome wide 
	       lociSelection == 'f' ){  //regions from a file 

		loopingCondition=(indexOfLoci<(v.size()));   
		if(loopingCondition)
		    rangeFound=v[indexOfLoci];    	    
	    }else{
		
		//for the entire genome, we do not need to loop
		if(lociSelection == 'v'){  //entire chr (only one window)
		    loopingCondition=false;
		}else{


		    cerr << "Unknown loci selection"<<endl;
		    return 1;       
		    break;
		}
	    }
	}

    
    }
    while(loopingCondition); // end of main while loop 



    if(mode == 't'){ 
    	delete(tmim);
    }

    cerr<<"Program terminated gracefully"<<endl;

    if(useBedFileRegions){
	map< string, vector<GenomicRange> * >::iterator iterMapBed;
	for (iterMapBed  = bedRegionsToFilter->begin(); 
	     iterMapBed != bedRegionsToFilter->end(); 
	     ++iterMapBed) {
	    delete(iterMapBed->second);	    
	}
	delete(bedRegionsToFilter);
    }

    if(mode == 't'){
	for(unsigned int indexCut=0;indexCut<filtersToUse->size();indexCut++){
	    delete( (*filtersToUse)[indexCut] );
	}
	delete( filtersToUse );
    }

if(useBedFileRegionsTabix){
    delete(bedFileTabix);
 }

    return 0;
    }
// END main

