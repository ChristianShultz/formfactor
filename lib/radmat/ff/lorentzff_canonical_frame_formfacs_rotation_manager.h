#ifndef LORENTZFF_CANONICAL_FRAME_FORMFACS_ROTATION_MANAGER_H
#define LORENTZFF_CANONICAL_FRAME_FORMFACS_ROTATION_MANAGER_H 

#include "formfactor_abs_base_cfg.h"
#include "hadron/irrep_util.h"
#include "radmat/utils/tensor.h"
#include "radmat/rotation_interface/rotation_interface.h"
#include "lorentzff_formfac_utils.h"
#include "mass_averager.h"
#include <complex>
#include <sstream>


// this does not use the following headder but put it in so all 
// lorentzff_canonical classes see it
#include "lorentzff_formfactor_abs_base_cfg.h"

// print on helicity 
// #define LORENTZFF_CAN_PRINT_HELICITY

namespace radmat
{

  template<class DerivedFF, 
    typename Data_t, 
    int J_left , 
    int J_right,
    class MASS_AVERAGER >  // not being used -- could imagine averaging over irrep masses
    struct FormFacRotationManager
    : public FFAbsBlockBase_t<Data_t>,
    public DMatrixManager
  {
    typedef Tensor<double,1> p4_t; 

    virtual ~FormFacRotationManager() {}

    virtual std::pair<mom_t,mom_t> 
      pair_mom(const p4_t &l, const p4_t &r, const double kick) const
      {
        return std::pair<mom_t,mom_t>(get_space_mom(l,kick),get_space_mom(r,kick)); 
      }

    virtual std::pair<p4_t,p4_t>
      can_mom(const RotationMatrix_t *R,
          const p4_t &l, 
          const p4_t &r,
          const double kick) const
      {
        typedef radmat::LatticeRotationEnv::TheRotationGroupGenerator RG; 
        Tensor<double,1> ll( (TensorShape<1>())[4], 0.);
        Tensor<double,1> rr( (TensorShape<1>())[4], 0.);
        std::pair<mom_t,mom_t> fmom = pair_mom(l,r,kick); 
        mom_pair_key frame = this->get_frame(fmom.first,fmom.second); 
        std::pair<mom_t,mom_t> cmom = frame.moms(); 

        mom_t cl = cmom.first; 
        mom_t cr = cmom.second; 

#ifdef BREAK_ROTATIONS_REST
        mom_pair_key k = mom_pair_key(fmom.first,fmom.second); 
        fmom = k.frame_orientation(); 
        cmom = frame.frame_orientation(); 
#endif 
        // NB: do the check with the frame orientation
        // are these the momentum i think that they are? 
        if ( !!! check_total_frame_transformation( R, fmom.first,fmom.second,cmom.first,cmom.second,true ) )
        {
          std::cout << __func__ << ": frame transformation error" << std::endl;

          std::cout << "int moms l" << string_mom(fmom.first) << " r " << string_mom(fmom.second)
            << " ll " << string_mom(cl) << " rr " << string_mom(cr) << std::endl;

          std::cout << "pl " << l << " pll " << ll << " pr " << r << " prr "
            << rr << " R " << *R << std::endl;
          throw std::string("frame remap error"); 
        }

        ll[0] = l[0];
        ll[1] = cl[0]*kick;
        ll[2] = cl[1]*kick;
        ll[3] = cl[2]*kick;

        rr[0] = r[0];
        rr[1] = cr[0]*kick;
        rr[2] = cr[1]*kick;
        rr[3] = cr[2]*kick;


        // MASS_AVERAGER::operate(ll,rr); 

        return std::pair<p4_t,p4_t>(ll,rr); 
      }

    virtual Tensor<Data_t , 1> 
      rotate(const RotationMatrix_t *R, 
          const Tensor<Data_t,1> &in) const
      {
        Tensor<Data_t,1> out( (TensorShape<1>())[4] , Data_t() ); 

        for(int i = 0; i < 4; ++i)
          for(int j = 0; j < 4; ++j)
            out[i] = out[i] +  (*R)[i][j] * in[j];

        return out; 
      }


    virtual Tensor<Data_t,1> 
      operator()(const MomRowPair_t &lefty, 
          const MomRowPair_t &righty, 
          const double mom_fac ) const 
      {
        return canonicalize(lefty,righty,mom_fac); 
      }


    virtual Tensor<Data_t,1>
      canonicalize(const MomRowPair_t &lefty,
          const MomRowPair_t &righty, 
          const double kick) const
      {
        std::pair<mom_t,mom_t> moms = pair_mom(lefty.first,righty.first,kick);  
        rCompEulAngles eul = rotation_matrix(moms.first,moms.second); 
        RotationMatrix_t *R = new RotationMatrix_t( genRotationMatrix(eul) );  
        std::pair<p4_t,p4_t> can_moms = can_mom(R,lefty.first,righty.first,kick); 

        WignerMatrix_t * Wl = left_wigner_matrix(eul,moms.first,moms.second,J_left);
        WignerMatrix_t * Wr = right_wigner_matrix(eul,moms.first,moms.second,J_right);

        //  WignerMatrix_t * Wl = left_wigner_matrix(eul,moms.first,moms.second,J_left,false,1);
        //  WignerMatrix_t * Wr = right_wigner_matrix(eul,moms.first,moms.second,J_right,false,1);

        this->conjugate(Wl); 
        this->conjugate(Wr); 

        Tensor<Data_t,1> sum( (TensorShape<1>())[4] , Data_t() ); 

        int left_h = J_left - lefty.second; 
        int right_h = J_right - righty.second; 
        int left_bound = 2*J_left +1; 
        int right_bound = 2*J_right +1; 

        for(int lh = 0; lh < left_bound; ++lh)
          for(int rh = 0; rh < right_bound; ++rh)
          {
            std::complex<double> weight = ( (*Wl)[lh][left_h] * (*Wr)[right_h][rh] );

            if( std::norm(weight) > 1e-6) 
            {

#ifdef LORENTZFF_CAN_PRINT_HELICITY
              std::cout << "**************************** " << __func__ 
                << ": constructing hl = " << J_left - lh 
                << " hr = " << J_right - rh << std::endl; 
              std::cout << "left_weight = " << (*Wl)[lh][left_h] 
                << "\nright_weight = " << (*Wr)[right_h][rh] << std::endl;
#endif 

              Tensor<Data_t, 1> tmp = this->impl(can_moms.first,can_moms.second,kick,J_left-lh,J_right-rh); 
              for(int i = 0; i < 4; ++i)
                sum[i] = sum[i] + weight * tmp[i] ; // do by hand for overloads
            }
          }

        Tensor<Data_t,1> ret = rotate(R,sum); 

        delete Wl;
        delete Wr; 
        delete R; 

        return ret; 
      }

    // derived classes implement impl, ff_impl, id_impl 
    virtual Tensor<Data_t, 1> 
      impl(const p4_t &l, const p4_t &r, const double kick, const int hl, const int hr) const = 0; 
    //      {
    //        return static_cast<const DerivedFF*>(this)->impl(l,r,kick,hl,hr); 
    //      } 

    virtual std::string 
      ff_impl(void) const 
      {
        return static_cast<const DerivedFF*>(this)->ff_impl(); 
      }

    virtual std::string id() const 
    {
      return Stringify<DerivedFF>(); 
    }

  }; 


} // radmat



#endif /* LORENTZFF_CANONICAL_FRAME_FORMFACS_ROTATION_MANAGER_H */
