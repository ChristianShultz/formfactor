#ifndef REDSTAR_CANONICAL_ROTATIONS_H
#define REDSTAR_CANONICAL_ROTATIONS_H 


/*
 *
 *    Return rotation matrix to take us from 
 *    p_z to p where p is an allowed lattice 
 *    momentum.
 *
 *    These rotations are a combination of 
 *    the reference rotation followed by the 
 *    lattice rotation 
 *
 *
 */


#include "radmat/utils/tensor.h"
#include "radmat/utils/handle.h"
#include "adat/singleton.h"
#include "adat/objfactory.h"
#include "io/adat_xmlio.h"


namespace radmat
{
  typedef Tensor<double,2> RotationMatrix_t; 
  typedef ADATXML::Array<int> mom_t; 

  namespace CanonicalRotationEnv
  {


    typedef Util::SingletonHolder<
      Util::ObjectFactory<
      RotationMatrix_t, 
      std::string, 
      void,
      RotationMatrix_t* (*)(void),
      Util::StringFactoryError>
        >
        TheCanonicalRotationFactory;


    // do reg
    bool registerAll(void); 

    // call using a string, (LG returns ref rots)
    RotationMatrix_t* call_factory(const std::string &);

    // call using a momentum, (returns Rp_z = p)
    RotationMatrix_t* call_factory(const mom_t &);

  } // CanonicalRotationEnv

} // radmat



#endif /* REDSTAR_CANONICAL_ROTATIONS_H */
