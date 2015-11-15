#ifndef LORENTZFF_CANONICAL_ROTATIONS_CHECKER_H
#define LORENTZFF_CANONICAL_ROTATIONS_CHECKER_H 


#include "radmat/construct_data/lattice_multi_data_object.h"
#include "radmat/utils/handle.h"


namespace radmat
{


  struct LatticeRotationRelationChecker
  {
    void check( const rHandle<LLSQLatticeMultiData> &) const; 
  }; 






} // radmat 



#endif /* LORENTZFF_CANONICAL_ROTATIONS_CHECKER_H */
