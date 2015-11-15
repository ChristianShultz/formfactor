#ifndef PHOTON_SUBDUCTION_H
#define PHOTON_SUBDUCTION_H 


#include "radmat/utils/obj_expr_t.h"
#include "io/adat_xmlio.h"
#include <complex>
  
// combine the polarization and subduction into a single step 

namespace radmat
{

  typedef ListObjExpr_t<std::complex<double> , int> 
    lorentz_to_cubic_t;  

  lorentz_to_cubic_t photon_subduction(const std::string &rep, 
      const std::string &spher_rep,  
      const int row, 
      const ADATXML::Array<int> &q); 

} // radmat 

#endif /* PHOTON_SUBDUCTION_H */
