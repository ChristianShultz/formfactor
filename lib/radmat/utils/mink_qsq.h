#ifndef MINK_QSQ_H
#define MINK_QSQ_H 

#include "io/adat_xmlio.h"

namespace radmat
{

    // send in some integer array of 3-momentum for initial and final guys
    //      as well as their REST mass and the size of a momentum unit 
    //
    //      spit out a double corresponding to the qsq 
   double Mink_qsq(const ADATXML::Array<int> &pf, const double mf, 
       const ADATXML::Array<int> &pi, const double mi,
       const double factor); 

  double mom_factor(const double xi, const int L_s) ;
}

#endif /* MINK_QSQ_H */
