#ifndef FORMFACTOR_ABS_BASE_CFG_H
#define FORMFACTOR_ABS_BASE_CFG_H 

#include "formfactor_utils.h"
#include "radmat/utils/stringify.h"
#include "radmat/utils/pow2assert.h"
#include "radmat/utils/tensor.h"
#include "radmat/utils/handle.h"
#include "ensem/ensem.h"
#include <complex>
#include <utility>
#include <sstream>
#include <string>
#include <list>


namespace radmat
{

  typedef std::pair< Tensor<double,1> , int > MomRowPair_t; 

  // this is a base class to set up list of handles 
  // to polymorphic functors to construct a generalization
  // of a matrix element
  //
  // actual classes will implement the operator() 
  //
  // a decomposition is a list of functors

  template<typename T> 
  struct FFAbsBlockBase_t
  {
    virtual ~FFAbsBlockBase_t() {}
    virtual std::string ff() const {return "FFAbsBlockBase_t";}
    virtual std::string id() const = 0; // the derived class name

    virtual Tensor<T, 1> 
      operator()( const MomRowPair_t &lefty, 
          const MomRowPair_t &righty,
          const double mom_fac) const = 0; 
  };

  /////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////


  // generate some linear least squares system w/o having to know very much about the 
  // matrix element invariants. These are still at the configuration level. Just 
  // shove this into a SembleMatrix bin by bin to deal with the ensemble stats.

  // this will be polymorphic with different named classes corresponding to 
  // different quantum numbers, ie: PiPi  <-->   < 0+ | j_mu | 0+ >

  struct FFAbsBase_t; 
  REGISTER_STRINGIFY_TYPE( FFAbsBase_t );


  // this is pure right now
  struct FFAbsBase_t
  {

    virtual ~FFAbsBase_t() {}
    virtual int nFacs(void) const = 0; 
    virtual std::string ff(void) const = 0; 
    virtual std::map<int,std::string> ff_ids(void) const = 0; 
    virtual std::string id() const { return Stringify<FFAbsBase_t>(); }

    virtual itpp::Mat<std::complex<double> > 
      operator()( const MomRowPair_t &lefty, 
          const MomRowPair_t &righty,
          const double mom_fac) const = 0; 


    // only actual methods
    virtual MomRowPair_t to_mom_row_pair( const double E, 
        const Array<double> &p, 
        const int row)
    {
      Tensor<double,1> mom( (TensorShape<1>())[4] , 0. );
      mom[0] = E;
      mom[1] = p[0];
      mom[2] = p[1];
      mom[3] = p[2];

      return MomRowPair_t(mom,row); 
    }

  };


} // radmat 

#endif /* FORMFACTOR_ABS_BASE_CFG_H */
