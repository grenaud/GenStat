/*
 * DstatResult
 * Date: Mar-15-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include "DstatResult.h"

DstatResult::DstatResult(){
    all.reinitializedCounters();
    noCpg.reinitializedCounters();
    onlyCpg.reinitializedCounters();
    transitions.reinitializedCounters();
    transversions.reinitializedCounters();
    // noDamage.reinitializedCounters();
}

DstatResult::DstatResult(const DstatResult & other){
     all             = other.all;
     noCpg           = other.noCpg;
     onlyCpg         = other.onlyCpg;
     transversions   = other.transversions;
     transitions     = other.transitions;
}

DstatResult::~DstatResult(){

}


string DstatResult::printWithBootstrap(list<vector< vector< vector<DstatResult> >  > >  & boostraps,unsigned int i,unsigned int j,unsigned int k,unsigned int numberOfBootstraps){
    stringstream toreturn;

    vector<double> allBoot;
    vector<double> nocpgBoot;
    vector<double> onlycpgBoot;
    vector<double> transitionsBoot;
    vector<double> transversionsBoot;

    //all
    for (list< vector< vector< vector<DstatResult> >  >   >::iterator it=boostraps.begin(); 
	 it != boostraps.end(); it++){
	allBoot.push_back(              (*it)[i][j][k].all.computeDST()           );
	nocpgBoot.push_back(            (*it)[i][j][k].noCpg.computeDST()         );
	onlycpgBoot.push_back(          (*it)[i][j][k].onlyCpg.computeDST()       );
	transitionsBoot.push_back(      (*it)[i][j][k].transitions.computeDST()   );
	transversionsBoot.push_back(    (*it)[i][j][k].transversions.computeDST() );
    }

    //cout<<vectorToString(allBoot)<<endl;
    pair<double,double> allDEV         = computeMeanSTDDEV(allBoot);
    pair<double,double> noCpgDEV           = computeMeanSTDDEV(nocpgBoot);
    pair<double,double> onlyCpgDEV         = computeMeanSTDDEV(onlycpgBoot);
    pair<double,double> transitionsDEV     = computeMeanSTDDEV(transitionsBoot);
    pair<double,double> transversionsDEV   = computeMeanSTDDEV(transversionsBoot);


    //cout<<allBootDEV.first<<"\t"<<allBootDEV.second<<endl;
    toreturn<<":\t"       <<all.headerForCount()  <<"\n";
    toreturn<<"all_sites:\t"      <<all           <<"\t"<<allDEV.first           <<"\t"<<allDEV.second           <<"\t"<<all.computeDST()           / allDEV.second <<"\n";
    toreturn<<"nocpg_sites:\t"    <<noCpg         <<"\t"<<noCpgDEV.first         <<"\t"<<noCpgDEV.second         <<"\t"<<noCpg.computeDST()         / noCpgDEV.second <<"\n";
    toreturn<<"onlycpg_sites:\t"  <<onlyCpg       <<"\t"<<onlyCpgDEV.first       <<"\t"<<onlyCpgDEV.second       <<"\t"<<onlyCpg.computeDST()       / onlyCpgDEV.second <<"\n";
    toreturn<<"transitions:\t"    <<transitions   <<"\t"<<transitionsDEV.first   <<"\t"<<transitionsDEV.second   <<"\t"<<transitions.computeDST()   / transitionsDEV.second <<"\n";
    toreturn<<"transversions:\t"  <<transversions <<"\t"<<transversionsDEV.first <<"\t"<<transversionsDEV.second <<"\t"<<transversions.computeDST() / transversionsDEV.second <<"\n";

    return toreturn.str();
}

DstatResult &  DstatResult::operator+=(const DstatResult & other){   

    this->all             += other.all;
    this->noCpg           += other.noCpg;
    this->onlyCpg         += other.onlyCpg;
    this->transversions   += other.transversions;
    this->transitions     += other.transitions;

    return *this;

}
