/*
 * DistanceResult
 * Date: Jan-26-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include "DistanceResult.h"

DistanceResult::DistanceResult(){
    all.reinitializedCounters();
    noCpg.reinitializedCounters();
    onlyCpg.reinitializedCounters();
    transitions.reinitializedCounters();
    transversions.reinitializedCounters();

}

DistanceResult::DistanceResult(const string stringForm){
    vector<string> allLines=allTokens(stringForm,'\n');

    all           = AllelePairCounter(allLines[1]);
    noCpg         = AllelePairCounter(allLines[2]);
    onlyCpg       = AllelePairCounter(allLines[3]);
    transitions   = AllelePairCounter(allLines[4]);
    transversions = AllelePairCounter(allLines[5]);

}

DistanceResult::~DistanceResult(){

}



std::ostream & operator << (std::ostream & os, const DistanceResult & dr){
    os<<"Category:\t"       <<dr.all.headerForCount()<<"\n";
    os<<"all_sites:\t"      <<dr.all          <<"\n";
    os<<"nocpg_sites:\t"    <<dr.noCpg        <<"\n";
    os<<"onlycpg_sites:\t"  <<dr.onlyCpg      <<"\n";
    os<<"transitions:\t"    <<dr.transitions   <<"\n";
    os<<"transversions:\t"  <<dr.transversions<<"\n";
    return os;
}

DistanceResult &  DistanceResult::operator+=(const DistanceResult & other){
    // cout<<"DR before "<<this->all<<endl;
    // cout<<"DR other  "<<other.all<<endl;

    this->all           += other.all;
    // cout<<"DR after "<<this->all<<endl;
    this->noCpg         += other.noCpg;
    this->onlyCpg       += other.onlyCpg;
    this->transitions   += other.transitions;
    this->transversions += other.transversions;
    
    return *this;
}
