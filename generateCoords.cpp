/*
 * generateCoords
 * Date: Oct-30-2012 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include <iostream>
#include <fstream>
#include <string.h>
#include <gzstream.h>

#include "RandomGenomicCoord.h"
#include "GenomicWindows.h"

#include "utils.h"

using namespace std;

int main (int argc, char *argv[]) {

    // case 'x':  //undefined
    // case 'f':  //reading regions from a file
    // case 'r':  //random windows
    // case 'o':  //entire genome per window
    // case 'c':  //entire chr per window	    
    // case 'w':  //genome wide 
    // case 'v':  //entire chr (only one window)
    char lociSelection='x';
    string fastaIndex  =  "/mnt/454/Altaiensis/users/gabriel/faidx/index.hg19.fai"  ;
    int bpToExtract =-1      ; //amount of bp per chunk
    int overlapBetweenWindows=0;
    string chrToUse;
    unsigned int amountOfGoodSitesTARGET=0;  //number of contiguous chunks to extract for random
    bool allowSexChr=false;

    bool specifiedChunk=false;
    bool specifiedSizeFile=false;
    string sizeFile;
    bool inBedFormat=false;

    string externalDatabase=string("")+
	"\t\t"+"--fai [file]" +"\t\t\t\t"+"Fasta index for genome (produced by \"samtools faidx\") (default: "+fastaIndex+")\n"+
	"\t\t"+"--allowsexchr" +"\t\t\t\t"+"Allow sites to be generated on X and Y (default: "+boolStringify(allowSexChr)+")\n";

    
    string sizeRegions=string("")+
	"\tLoci size option:\n"+
	"\t\t"+"--chunk [size]"     +"\t\t\t\t"+"Size of contiguous genomic region to take (default: none)\n"+
	"\t\t"+"--size  [file]"     +"\t\t\t\t"+"Instead of specifying --chunk for --random, read this file with the sizes (one per line, no header, just integers)\n"+
	"\t\t"+"--overlap [size]"   +"\t\t\t"+  "Size of overlap between windows (for overlap mode) (default: "+stringify(overlapBetweenWindows)+")\n";
    

    string selection=string("")+
	"\tLoci selection mode, select either:\n"+
	// "\t\t"+"--region [file]"    +"\t\t\t\t"   +"Read [file] in BED format to select which loci to use\n"+
	"\t\t"+"--random [number]"  +"\t\t\t"     +"Select [number] random loci from the genome\n"+
	"\t\t"+"--loci"             +"\t\t\t\t\t" +"Select multiple loci from the genome\n"+
	"\t\t"+"--locichr [chr]"    +"\t\t\t\t"   +"Select multiple loci from a chromosome\n"+
	"\t\t"+"--chr [chr]"        +"\t\t\t\t"   +"Run on an entire chromosome\n"+
	"\t\t"+"--wide"             +"\t\t\t\t\t" +"Run on the entire genome\n";

  string outputOpt=string("")+
      "\tOutput options:\n"+
      "\t\t"+"--bed"  +"\t\t\t"     +"Produce output in bed format (default: "+boolStringify(inBedFormat)+")\n";


    if(argc == 1 ||
       (argc == 3 && (string(argv[2]) == "-h" || string(argv[2]) == "--help") )
       ){
	cerr << "Usage:  "<<string(argv[0])<<" [options]"<<endl<<endl
	     <<"This program can generate coordinates for a given genome"<<endl
	     <<"Options:"<<endl
	     <<"\n\tLocation of database files:\n"	    
	     <<externalDatabase<<endl
	     <<sizeRegions<<endl
	     <<selection<<endl;
	return 1;       
    }


    for(int i=1;i<(argc);i++){ 

        if(strcmp(argv[i],"--bed") == 0 ){
	    inBedFormat=true;
            continue;
        }

        if(strcmp(argv[i],"--fai") == 0 ){
	    fastaIndex=string(argv[i+1]);
	    i++;
            continue;
        }

	//loci size
        if(strcmp(argv[i],"--chunk") == 0 ){
	    bpToExtract=destringify<int>(argv[i+1]);
	    i++;
	    specifiedChunk=true;
            continue;
        }

        if(strcmp(argv[i],"--size") == 0 ){
	    sizeFile=string(argv[i+1]);
	    specifiedSizeFile=true;
	    i++;
            continue;
        }

	if(strcmp(argv[i],"--overlap") == 0 ){
	    overlapBetweenWindows=destringify<int>(argv[i+1]);
	    i++;
            continue;
        }


	if(strcmp(argv[i],"--allowsexchr") == 0 ){
	    allowSexChr=true;
	    continue;
	}

	//BEGIN loci selection
        if(strcmp(argv[i],"--random") == 0 ){
	    if(lociSelection != 'x'){
		cerr << "Choose only one loci selection option"<<endl;
		return 1;       
	    }
	    amountOfGoodSitesTARGET=destringify<unsigned int>(argv[i+1]);
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

        // if(strcmp(argv[i],"--region") == 0 ){
	//     if(lociSelection != 'x'){
	// 	cerr << "Choose only one loci selection option"<<endl;
	// 	return 1;       
	//     }
	//     bedFileWithRegions=string(argv[i+1]);
	//     i++;	    
	//     lociSelection='f';
        //     continue;
        // }


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


	cerr<<"Unknown option"<<endl;
	return 1;
    }

    if(specifiedChunk && specifiedSizeFile){
	cerr << "Error, cannot specify both --size and --chunk"<<endl;
	return 1;       
    }

    if(lociSelection != 'r' &&
       specifiedSizeFile ){ //random
	cerr << "The option --size is only for random genomic region"<<endl;
	return 1;       	
    }

    if(lociSelection != 'v' &&  //entire chr
       lociSelection != 'w' ){	 //genome wide
	if(!specifiedChunk && !(lociSelection=='r' && specifiedSizeFile) ){
	    cerr << "Need to specify the size (with --chunk) of the size of contiguous genomic region to take"<<endl;
	    return 1;       
	}

    }


    vector<int> sizeOfBp;

    if(specifiedSizeFile){

	string line;
	igzstream myFile;
	string filename = string(sizeFile);
	myFile.open(filename.c_str(), ios::in);
	
	if (myFile.good()){
	    while ( getline (myFile,line)){
		int toadd=destringify<int>(line);
		if(toadd<0){
		    cerr << "Cannot have negative values in file for --size line="<<line<<endl;
		    return 1;     
		}
		    
		sizeOfBp.push_back(toadd);
	    }
	    myFile.close();
	}else{
	    cerr << "Unable to open file "<<filename<<endl;
	    return 1;
	}

	if(sizeOfBp.empty()){
	    cerr << "Error:  file "<<filename<<" does not contain any size"<<endl;
	    return 1;
	}
	
    }


    RandomGenomicCoord rgc (fastaIndex,allowSexChr);
    GenomicWindows     rw  (fastaIndex,allowSexChr);
    vector<GenomicRange> v;
    unsigned int indexOfLoci=0;
    switch(lociSelection){
    case 'r': //random loci
	
	while(indexOfLoci<amountOfGoodSitesTARGET){
	    if(specifiedSizeFile){
		//srand is called in RandomGenomicCoord
		unsigned int toReturn = (unsigned int)(mrand48());
		bpToExtract=sizeOfBp[toReturn%sizeOfBp.size()];
	    }
	    // v.push_back(rgc.getRandomGenomic(bpToExtract));
	    indexOfLoci++;
	    if(!inBedFormat)
		cout<<rgc.getRandomGenomic(bpToExtract)<<endl;
	    else
		cout<<rgc.getRandomGenomic(bpToExtract).asBed()<<endl;
	}
	cerr<<"Program terminated successfully"<<endl;
	return 0;
	break;
    case 'c':  //windows on entire chr
	v=rw.getGenomicWindowsChr(chrToUse,bpToExtract,overlapBetweenWindows);
	if( v.size() == 0 ){    cerr<<"No range found using these chr/loci settings"<<endl; return 1;}
	break;
    case 'o':    //windows on entire genome
	v=rw.getGenomicWindows(bpToExtract,overlapBetweenWindows);
	if( v.size() == 0 ){    cerr<<"No range found using these chr/loci settings"<<endl; return 1;}
	break;
    case 'v':  //entire chr
	v=rw.getChr(chrToUse);
	if( v.size() == 0 ){    cerr<<"No range found using these chr/loci settings"<<endl; return 1;}
	break;
    case 'w':  //entire genome
	v=rw.getGenomeWide();
	if( v.size() == 0 ){    cerr<<"No range found using these chr/loci settings"<<endl; return 1;}
	break;
    default:
	cerr << "Unknown loci selection"<<endl;
	return 1;       
	break;
    }


    for(unsigned int i=0;i<v.size();i++){
	if(!inBedFormat)
	    cout<<v[i]<<endl;
	else
	    cout<<v[i].asBed()<<endl;
    }



    return 0;
}

