

#include "redstar_cartesian_interface.h"
#include "redstar_particle_handler_utils.h"
#include "radmat/utils/pow2assert.h"
#include "radmat/utils/tensor.h"
#include "redstar_photon_polarization_tensor.h"
#include "itpp/itbase.h"
#include "semble/semble_meta.h"
#include <complex>


namespace radmat
{
  //////////////////////////////////////////////////////////////////
  EnsemRedstarBlock 
    operator*(const std::complex<double> &c,
        const EnsemRedstarBlock &l)
    {
      return SEMBLE::toScalar(c)*l; 
    }

  //////////////////////////////////////////////////////////////////
  itpp::Vec<EnsemRedstarBlock> 
    operator*(const std::complex<double> &d,
        const itpp::Vec<EnsemRedstarBlock> &o)
    {
      itpp::Vec<EnsemRedstarBlock> foo(o.size()); 
      for(int i = 0; i < o.size(); ++i)
        foo[i] = d * o[i];
      return foo; 
    }



    //////////////////////////////////////////////////////////////////
    itpp::Vec<EnsemRedstarBlock> 
    operator*(const itpp::Mat<std::complex<double> > &m,
        const itpp::Vec<EnsemRedstarBlock> &v)
    {
      POW2_ASSERT(v.size() == m.cols()); 
      itpp::Vec<EnsemRedstarBlock> ret(m.rows()); 

      for(int row = 0; row < m.rows(); ++row)
        for(int col = 0; col < m.cols(); ++col)
        {
          if(std::norm(m(row,col)) > 0.0001)
          {
            // init 
            if(ret[row].m_expr.size() == 0)
              ret[row] =  m(row,col)*v(col);
            else
              ret[row] = ret[row] + m(row,col)*v(col);
          }
        }

      return ret; 
    }


  //////////////////////////////////////////////////////////////////
  // multiply a matrix against a vector of bools
  itpp::Vec<bool>
    operator*(const itpp::Mat<std::complex<double> > &m,
        const itpp::Vec<bool> &v)
    {
      POW2_ASSERT(v.size() == m.cols()); 
      itpp::Vec<bool> ret(m.rows());

      for(int row = 0; row < ret.size(); ++row)
        ret(row) = true; 


      for(int row = 0; row < m.rows(); ++row)
        for(int col = 0; col < v.size(); ++col)
        {
          // phases may not cancel exactly 
          // but in this context 0.000001 is zero
          if(std::norm(m(row,col)) > 0.0001)
            continue; 

          ret(row) &= v(col);
        }

      return ret; 
    }


  //////////////////////////////////////////////////////////////////
  // rows correspond to a helicity (+ 0 -)
  // cols are the cartesian coordinate (x,y,z)
  //         -->  j^lambda = M * j^cartesian 
  itpp::Mat<std::complex<double> > 
    eps3d(const ADATXML::Array<int> &mom ,
        const bool create)
  {
    Tensor<std::complex<double>, 1 > tmp;
    genPolTens3D<1> eps(mom);
    itpp::Mat<std::complex<double> > eps3(3,3); 

    for(int h = 1; h > -2; --h)
    {
      tmp = eps(h);

      for(int i = 0; i < 3; ++i)
        if(create)
          eps3(1-h,i) = tmp[i];
        else
          eps3(1-h,i) = std::conj(tmp[i]); 

    }

    return eps3; 
  }


  //////////////////////////////////////////////////////////////////
  // invert the matrix eps -- return the dagger of the transformation 
  // to helicity since this has to be a unitary transformation 
  //    ie: it is only a basis change
  itpp::Mat<std::complex<double> > 
    invert2Cart(const ADATXML::Array<int> mom, 
        const bool create)
  { 
    return itpp::hermitian_transpose(eps3d(mom,create));
  }

} // radmat
