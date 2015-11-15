#ifndef REDSTAR_CARTESIAN_INTERFACE_H
#define REDSTAR_CARTESIAN_INTERFACE_H 


#include "itpp/itbase.h"
#include "redstar_particle_handler_utils.h"


namespace radmat
{
  // multiplication overload
  EnsemRedstarBlock operator*(const std::complex<double> &c,
      const EnsemRedstarBlock &l);

  // multiplication overload
  itpp::Vec<EnsemRedstarBlock>
    operator*(const std::complex<double> &, 
        const itpp::Vec<EnsemRedstarBlock> &);

  // mul out a matrix against a vector of ERBs
  itpp::Vec<EnsemRedstarBlock> 
    operator*(const itpp::Mat<std::complex<double> > &m,
        const itpp::Vec<EnsemRedstarBlock> &v);

  // multiply a matrix against a vector of bools
  itpp::Vec<bool>
    operator*(const itpp::Mat<std::complex<double> > &m, const itpp::Vec<bool> &v);


  // rows correspond to a helicity (+ 0 -)
  // cols are the cartesian coordinate (x,y,z)
  // then M(row,col) = epsilon_cart(lambda)  -->  j^lambda = M * j^cartesian 
  itpp::Mat<std::complex<double> >
    eps3d(const ADATXML::Array<int> &mom , const bool create);


  // invert the matrix eps
  itpp::Mat<std::complex<double> >
    invert2Cart(const ADATXML::Array<int> mom, const bool create);


} // radmat 

#endif /* REDSTAR_CARTESIAN_INTERFACE_H */
