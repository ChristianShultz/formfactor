/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

* File Name : lorentzff_formfac_utils.cc

* Purpose :

* Creation Date : 12-12-2013

* Last Modified : Thu 09 Jan 2014 05:22:42 PM EST

* Created By : shultz

_._._._._._._._._._._._._._._._._._._._._.*/

#include "lorentzff_formfac_utils.h"
#include <cmath>

namespace radmat
{
  int 
    round_double_nearest_int(const double &d)
    {
     // return d < 0 ? -1*int(fabs(d) + 0.5 ) : int(fabs(d) +0.5); 
     return floor( d + 0.5 ) ; 
    }

  double 
    round_to_zero(const double &d, const double thresh)
    {
      if(fabs(d) < thresh )
        return 0; 
      return d; 
    } 

  std::complex<double>
    round_to_zero(const std::complex<double> &d, const double thresh)
    {
      return std::complex<double>( round_to_zero(d.real(),thresh),
          round_to_zero(d.imag(),thresh) ); 
    }


  mom_t 
    get_space_mom(const Tensor<double, 1> &p, const double mom_kick)
    {
      mom_t r = gen_mom<0,0,0>(); 
      r[0] = round_double_nearest_int(p[1]/mom_kick); 
      r[1] = round_double_nearest_int(p[2]/mom_kick); 
      r[2] = round_double_nearest_int(p[3]/mom_kick); 
      return r; 
    }

}
