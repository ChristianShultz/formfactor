#ifndef SPECTRUM_3F_743_H
#define SPECTRUM_3F_743_H 

#include "spectrum_state.h"

namespace radmat
{

  namespace L_3F_743_StateFactoryEnv
  {

    struct sPion   : public SpectrumState<J0mRep_t, PI,   L_3F_743, 11, 1483, enumStateHandler> {};
    struct sPionP  : public SpectrumState<J0mRep_t, PIp,  L_3F_743, 11, 3683, enumStateHandler> {};
    struct sPionPP : public SpectrumState<J0mRep_t, PIpp, L_3F_743, 11, 4563, enumStateHandler> {};


    struct sRho    : public SpectrumState<J1mRep_t, RHO,    L_3F_743, 11, 2162, enumStateHandler> {};
    struct sRhoP   : public SpectrumState<J1mRep_t, RHOp,   L_3F_743, 11, 3959, enumStateHandler> {};
    struct sRhoPP  : public SpectrumState<J1mRep_t, RHOpp,  L_3F_743, 11, 4207, enumStateHandler> {};
    struct sRhoPPP : public SpectrumState<J1mRep_t, RHOppp, L_3F_743, 11, 4867, enumStateHandler> {};

    struct sA1 : public SpectrumState<J1pRep_t, A1, L_3F_743, 11, 3222, enumStateHandler> {};

    bool registerAll(); 
  }

}


#endif /* SPECTRUM_3F_743_H */
