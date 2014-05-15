/*
 * DivergenceResult
 * Date: Jan-24-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include "DivergenceResult.h"

DivergenceResult::DivergenceResult(){
    all.reinitializedCounters();
    noCpg.reinitializedCounters();
    onlyCpg.reinitializedCounters();
    transitions.reinitializedCounters();
    transversions.reinitializedCounters();
    noDamage.reinitializedCounters();
}

DivergenceResult::~DivergenceResult(){

}



string DivergenceResult::getHeader(){
    //    return "noMut\tcommon\tind1Spec\tind2Spec\tdivInd1\tdivInd1Low\tdivInd1High";
    return 
	all.getHeader("all")+"\t"+
	onlyCpg.getHeader("justCpg")+"\t"+
	noCpg.getHeader("noCpg")+"\t"+
	transitions.getHeader("transi")+"\t"+
	transversions.getHeader("transv")+"\t"+
	noDamage.getHeader("noDam");
}
