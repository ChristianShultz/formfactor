#ifndef POLARISATION_TENSORS_H_H_GUARD
#define POLARISATION_TENSORS_H_H_GUARD

#include "tensor.h"
#include "aux.h"
#include "radmat/utils/handle.h"
#include "ensem/ensem.h"
#include "hadron/irrep_util.h"
#include "hadron/clebsch.h"
#include "itpp/itbase.h"
#include "xml_array.h"
#include <utility>
#include <map>





/*
NB: These rotate like Helicity creation operators or equivalently 
like helicity states using the rotation conventions from ADAT.

They are 4-dimensional lorentz vectors 
 */





namespace radmat
{

  //! one remapping to rule them all
  int remapHelicity_1based(const int lambda, const int J); 


  // these routines use +/- helicities ( not 1 based )
  Tensor<std::complex<double> , 1> creation_op_J1_3(const ADATXML::Array<int> &mom, const int lambda); 
  Tensor<std::complex<double> , 1> creation_op_J1_3z(const int lambda); 

  Tensor<std::complex<double> , 1> creation_op_J1_4(const ADATXML::Array<int> &mom, 
    const int lambda, 
    const double Energy, 
    const double momentum_factor);

  Tensor<std::complex<double> , 1> creation_op_J1_4z(const double mod_p, 
    const int lambda, 
    const double Energy);

  typedef std::pair<idx_t, short> pKey_t; // J, hel

} // radmat

#endif
