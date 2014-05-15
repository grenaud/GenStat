/*
 * ComputeDivergence_core
 * Date: Jan-24-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef ComputeDivergence_core_h
#define ComputeDivergence_core_h

#include "DivergenceResult.h"

using namespace std;



inline bool isTransition(char allel_sample,char allel_reference,char allel_chimpHumanAncestor);
inline bool     isDamage(char allel_sample,char allel_reference,char allel_chimpHumanAncestor);
void computeDiv(char allel_chimpHumanAncestor,char allel_reference,char allel_sample, bool isCpg, DivergenceResult * divr);

/* class ComputeDivergence_core{ */
/* private: */
/* public: */
/* ComputeDivergence_core(); */
/* ~ComputeDivergence_core(); */
/* }; */

#endif
