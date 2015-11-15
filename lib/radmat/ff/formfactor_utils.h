#ifndef FORMFACTOR_UTILS_H
#define FORMFACTOR_UTILS_H 



#include "radmat/utils/pow2assert.h"
#include "radmat/utils/tensor.h"
#include "radmat/utils/handle.h"
#include "itpp/itbase.h"
#include "io/adat_xmlio.h"




using namespace ADATXML;
//using namespace ADATIO;


namespace radmat
{

  template<typename T> Tensor<T,1> 
    pPlus(const Tensor<T,1> &p_f, const Tensor<T,1> &p_i)   // p_f + p_i
    {
      return p_f + p_i; 
    }

  template<typename T> Tensor<T,1>
  pMinus(const Tensor<T,1> &p_f, const Tensor<T,1> &p_i)  // p_f - p_i
  {
    return p_f - p_i; 
  }


  /////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////

  //! copy a rank one tensor into an itpp vector for easy manipulation of lorentz 4-vectors
  template<typename T>
    itpp::Vec<T> toItpp(const Tensor<T,1> &v)
    {
      idx_t dim = v.getDim(0);
      itpp::Vec<T> foo(dim);
      foo.zeros(); 

      for(idx_t i = 0; i < dim; i++)
        foo[i] = v[i];

      return foo;
    }

  template<typename T>
    Tensor<T,1> toTensor(const itpp::Vec<T> &v)
    {
      idx_t dim = v.size();
      Tensor<T,1> foo((TensorShape<1>())[dim],0.);
      for(idx_t i = 0; i < dim; i++)
        foo[i] = v[i];

      return foo;
    }


  //! copy a rank 2 tensor into an itpp mat 
  //  -- don't know why I bothered writing this but it may prove to be useful
  template<typename T>
    itpp::Mat<T> toItpp(const Tensor<T,2> &m)
    {
      idx_t row = m.getDim(0);
      idx_t col = m.getDim(1);
      itpp::Mat<T> foo(row,col);

      for(idx_t i = 0; i < row; i++)
        for(idx_t j = 0; j < col; j++)
          foo(i,j) = m[i][j];

      return foo;
    }


} // radmat 


#endif /* FORMFACTOR_UTILS_H */
