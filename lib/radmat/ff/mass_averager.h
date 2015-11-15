#ifndef MASS_AVERAGER_H
#define MASS_AVERAGER_H 

#include "radmat/utils/tensor.h"

namespace radmat
{

  struct AVERAGE_MASSES
  {
    static void operate(Tensor<double,1> &l, Tensor<double,1> &r)
    {}
    //  {
    //    double ml = sqrt( l[0]*l[0] - l[1]*l[1] - l[2]*l[2] - l[3]*l[3]); 
    //    double mr = sqrt( r[0]*r[0] - r[1]*r[1] - r[2]*r[2] - r[3]*r[3]); 

    //    double m = 0.5*( ml + mr ); 

    //    l[0] = sqrt( m*m + l[1]*l[1] + l[2]*l[2] + l[3]*l[3]); 
    //    r[0] = sqrt( m*m + r[1]*r[1] + r[2]*r[2] + r[3]*r[3]); 
    //  }
  };

  struct DONT_AVERAGE_MASSES
  {
    static void operate(Tensor<double,1> &l, Tensor<double,1> &r)
    { }
  };


} // radmat


#endif /* MASS_AVERAGER_H */
