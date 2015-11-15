#ifndef REDSTAR_CANONICAL_LATTICE_ROTATIONS_H
#define REDSTAR_CANONICAL_LATTICE_ROTATIONS_H 


/*
 *
 *  Return lattice rotation matricies,
 *  the rotation that takes us from 
 *  p_ref to p, where both p_ref and p
 *  are allowed lattice rotations
 *
 */

#include "redstar_canonical_rotations.h"

namespace radmat
{

  namespace CanonicalLatticeRotationEnv
  {

    typedef Util::SingletonHolder<
      Util::ObjectFactory<
      RotationMatrix_t,
      std::string, 
      void,
      RotationMatrix_t* (*)(void),
      Util::StringFactoryError>
        >
        TheCanonicalLatticeRotationFactory;


    bool registerAll(void); 

    // returns the lattice rotation
    RotationMatrix_t* call_factory(const mom_t &rot); 


  } // CanonicalLatticeRotationEnv

} // radmat



#endif /* REDSTAR_CANONICAL_LATTICE_ROTATIONS_H */
