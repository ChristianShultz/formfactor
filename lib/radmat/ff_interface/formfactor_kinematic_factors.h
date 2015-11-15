#ifndef FORMFACTOR_KINEMATIC_FACTORS_H
#define FORMFACTOR_KINEMATIC_FACTORS_H 


#include "formfactor.h"
#include "semble/semble_vector.h"
#include "semble/semble_matrix.h"
#include <complex>

namespace radmat
{


  //
  // hide the transformation from tags to decomps behind this class 
  //
  //    the tag can be anything, the cc file in here must know how
  //    to deal with the derived type 
  //
  struct FFKinematicFactors_t
  {
    typedef SEMBLE::SembleMatrix<std::complex<double> > KinematicFactorMatrix;
    typedef SEMBLE::SembleVector<std::complex<double> > KinematicFactorRow;   
  

    virtual ~FFKinematicFactors_t(void) {} // handle cleans itself up

    virtual KinematicFactorRow genFactors(const DataTagPrimitive *);
    virtual KinematicFactorMatrix genFactorsMat(const DataTagPrimitive *); 
  };

} // radmat

#endif /* FORMFACTOR_KINEMATIC_FACTORS_H */
