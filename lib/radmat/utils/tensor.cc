// tensor.cc -
//
// Thursday, May 10 2012
//

#include "tensor.h"
#include "hadron/irrep_util.h"

namespace radmat
{

  Tensor<double,2> g_dd(void)
  {
    Tensor<double,2> gmunu((TensorShape<2>())[4][4],0.);
    gmunu[0][0] = 1.;
    gmunu[1][1] = -1.;
    gmunu[2][2] = -1.;
    gmunu[3][3] = -1;
    gmunu.lower_index(0);
    gmunu.lower_index(1);
    return gmunu;
  }

  Tensor<double,2> g_uu(void)
  {
    Tensor<double,2> gmunu((TensorShape<2>())[4][4],0.);
    gmunu[0][0] = 1.;
    gmunu[1][1] = -1.;
    gmunu[2][2] = -1.;
    gmunu[3][3] = -1;
    return gmunu;
  }
}
