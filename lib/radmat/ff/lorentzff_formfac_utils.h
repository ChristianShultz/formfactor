#ifndef LORENTZFF_FORMFAC_UTILS_H
#define LORENTZFF_FORMFAC_UTILS_H 

#include "io/adat_xmlio.h"
#include <complex>
#include "radmat/utils/tensor.h"
#include "radmat/rotation_interface/rotation_interface.h"

namespace radmat
{
  int 
    round_double_nearest_int(const double &d);

  double 
    round_to_zero(const double &d, const double thresh=1e-6); 

  std::complex<double>
    round_to_zero(const std::complex<double> &d, const double thresh=1e-6); 


  mom_t 
    get_space_mom(const Tensor<double, 1> &p, const double mom_kick);
}

#endif /* LORENTZFF_FORMFAC_UTILS_H */
