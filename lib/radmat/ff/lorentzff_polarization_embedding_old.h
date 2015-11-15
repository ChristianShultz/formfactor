#ifndef LORENTZFF_POLARIZATION_EMBEDDING_OLD_H
#define LORENTZFF_POLARIZATION_EMBEDDING_OLD_H

#include <complex>

#include "lorentzff_formfac_utils.h"
#include "ensem/ensem.h"
#include "semble/semble_meta.h"
#include "hadron/clebsch.h"
#include "hadron/irrep_util.h"


namespace radmat
{

  ////////////////////////////////////////////////////////////////////////
  //
  //  embed target helicity
  //
  ////////////////////////////////////////////////////////////////////////
  template<idx_t J, int hel>
    struct embedHelicityPolarizationTensor_old
    {
      virtual ~embedHelicityPolarizationTensor_old() {}

      // do work 
      Tensor<std::complex<double>, J> 
        z_axis_helicity_tensor(const Tensor<double,1> &t,   // target 
            const double mom_factor)
        {
          double mod_p = sqrt( t[1]*t[1] + t[2]*t[2] + t[3]*t[3] );
          ZAxisHelicityTensor<J> foo;

          Tensor<std::complex<double> , J> tens = foo(mod_p,hel,t[0]); 

          typename Tensor<std::complex<double> , J>::iterator it;
          for(it = tens.begin(); it != tens.end(); ++it)
            *it = round_to_zero(*it); 

          return tens; 
        }

      virtual Tensor<std::complex<double>, J>
        rotate(const Tensor<std::complex<double>,J> &z_axis,
            const rHandle<Tensor<double,2> > &R)
        {
          Tensor<std::complex<double>,J> rot(z_axis); 
          Tensor<std::complex<double>,2> Rloc;
          Rloc = convertTensorUnderlyingType<std::complex<double>,double,2>(*R); 

          // rotate all indicies 
          for(int i = 0; i < J; ++i)
            rot = contract(Rloc,rot,1,i); 

          // clean up 
          typename Tensor<std::complex<double> , J>::iterator it;
          for(it = rot.begin(); it != rot.end(); ++it)
            *it = round_to_zero(*it); 

          return rot; 
        }


      virtual Tensor<std::complex<double> , J> 
        conjugate(const Tensor<std::complex<double> , J> &inp) const
        {
          Tensor<std::complex<double> , J> foo = inp;
          typename Tensor<std::complex<double> , J>::iterator it;

          for(it = foo.begin(); it != foo.end(); ++it)
            *it = std::conj(*it); 

          return foo;
        }
    };



  //  Classes that need polarization tensors should derive from the 
  //  leftPTensor and rightPTensor class to enforce consistency 
  //
  //  take care of the left/right who is complex convention
  //
  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////


  ////////////////////////////////////////////////////////////////////////
  //  left polarization tensors (final state, annih ops)
  template<idx_t J_left, int hel_left>
    struct leftPTensor_old
    {

      virtual ~leftPTensor_old() {}

      virtual Tensor<std::complex<double> , J_left > 
        left_p_tensor(const Tensor<double,1> &lefty, 
            const Tensor<double,1> &righty, 
            const double mom_factor) const
        {
          embedHelicityPolarizationTensor_old<J_left,hel_left> foo; 
          Tensor<std::complex<double> , J_left> eps_z = foo.z_axis_helicity_tensor(lefty,mom_factor); 
          mom_t l = get_space_mom(lefty,mom_factor);
          mom_t r = get_space_mom(righty,mom_factor); 
          rHandle<Tensor<double,2> > Rcomposite(radmat::CanonicalRotationEnv::call_factory(l));
          Tensor<std::complex<double>, J_left> eps = foo.rotate(eps_z,Rcomposite); 

          return foo.conjugate( eps ); 
        }
    };


  ////////////////////////////////////////////////////////////////////////
  //  right polarization tensors (initial state, creation ops)
  template<idx_t J_right, int hel_right>
    struct rightPTensor_old
    {

      virtual ~rightPTensor_old() {}

      virtual Tensor<std::complex<double> , J_right > 
        right_p_tensor(const Tensor<double,1> &lefty, 
            const Tensor<double,1> &righty, 
            const double mom_factor) const
        {
          embedHelicityPolarizationTensor_old<J_right,hel_right> foo; 
          Tensor<std::complex<double> , J_right> eps_z = foo.z_axis_helicity_tensor(righty,mom_factor); 
          mom_t l = get_space_mom(lefty,mom_factor);
          mom_t r = get_space_mom(righty,mom_factor); 
          rHandle<Tensor<double,2> > Rcomposite(radmat::CanonicalRotationEnv::call_factory(r));

          return foo.rotate(eps_z,Rcomposite); 
        }
    };

}



#endif /* LORENTZFF_POLARIZATION_EMBEDDING_OLD_H */
