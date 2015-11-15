#ifndef LORENTZFF_PIPISTAR_H_H_GUARD
#define LORENTZFF_PIPISTAR_H_H_GUARD


#include "lorentzff_polarization_embedding.h"
#include "radmat/utils/pow2assert.h"
#include "radmat/utils/stringify.h"
#include "lorentzff_canonical_frame_formfacs_rotation_manager.h"
#include <complex>

namespace radmat
{

  struct PiPiStarF1;
  REGISTER_STRINGIFY_TYPE( PiPiStarF1 ); 

  // only one ff
  struct PiPiStarF1
    : public FormFacRotationManager<PiPiStarF1, std::complex<double> , 0 , 0 , DONT_AVERAGE_MASSES >
  {
    virtual ~PiPiStarF1() {}

    virtual  std::string ff_impl() const
    {
      return std::string("F_1(Q^2) \\left( p_+^{\\mu} + \\frac{m_l^2 - m_r^2}{Q^2}p_i^{\\mu} \\right)"); 
//      return std::string("F_1(Q^2)\\left(- p_+^{\\mu}\frac{Q^2}{m_{\\pi*}^2 -m_{\\pi}^2} + p_-\\right)");
    }

    // return a complex version of p_+
    virtual Tensor<std::complex<double> , 1> 
      impl(const Tensor<double,1> &p_f,
          const Tensor<double,1> &p_i,
          const double mom_fac,
          int h_f,
          int h_i) const
      {
        Tensor<std::complex<double>,2> metric; 
        Tensor<std::complex<double>,1> pp,pm;
        Tensor<std::complex<double>,0> num,denom; 

        metric = convertTensorUnderlyingType<std::complex<double>,double,2>(g_dd()); 
        pm = convertTensorUnderlyingType<std::complex<double>,double,1>(pMinus(p_f,p_i));
        pp = convertTensorUnderlyingType<std::complex<double>,double,1>(pPlus(p_f,p_i));

        num = contract(pm,applyMetric(pp,metric,0),0,0);
        denom = contract(pm,applyMetric(pm,metric,0),0,0);

        // Q2 = 0 is a problem 
        // if( std::norm(denom.value()) < 1e-14 )
        //  denom.value() = std::complex<double>(1e-14,0);  

        POW2_ASSERT_DEBUG(std::norm(denom.value()) > 1e-14);

        std::complex<double> coeff = -num.value() / denom.value(); 

        return coeff * pm + pp;  
      }
  };




  // generate a list for the PiPiStar constructor
  template<int embedl, int embedr>
    LorentzFFAbsBase_t::LorentzFFAbs_list PiPiStarGenList()
    {
      LorentzFFAbsBase_t::LorentzFFAbs_list retPiPiStar;
      LorentzFFAbsBase_t::BBType *blockPtr;
      blockPtr = new PiPiStarF1();
      POW2_ASSERT(blockPtr);  // blow up if something went wrong
      retPiPiStar.push_back(LorentzFFAbsBase_t::BBHandle_t(blockPtr));
      return retPiPiStar;
    }


  ////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////

  template<int embedl, int embedr> struct PiPiStar; 
  REGISTER_STRINGIFY_TYPE2( PiPiStar<0,0> ); 

  // only need to derive the constructor, everything else is in the 
  // base class, do this polymorphically, make some function that 
  // returns the appropriate handle based on the requested matrix element type
  template<int embedl, int embedr>
    struct PiPiStar : public LorentzFFAbsBase_t
  {
    PiPiStar()
      : LorentzFFAbsBase_t(radmat::PiPiStarGenList<embedl,embedr>())  
    {  } 

    PiPiStar& operator=(const PiPiStar &o)
    {

      if(this != &o)
        LorentzFFAbsBase_t::operator=(o);
      return *this;
    }

    // no slicing
    PiPiStar(const PiPiStar &o)
      : LorentzFFAbsBase_t(o)
    {  }

    virtual std::string reg_id() const { return Stringify< PiPiStar<embedl,embedr> >(); }
    virtual int left_spin() const {return embedl;}
    virtual int right_spin() const {return embedr;}
    virtual LorentzFFAbsBase_t * clone() const { return new PiPiStar(); }

    private:
    // I'm not sure if these could inherit so we will hide them as well
    PiPiStar(const LorentzFFAbsBase_t::LorentzFFAbs_list &);
    PiPiStar(const LorentzFFAbsBase_t::LorentzFFAbs_list);
  };

} // close radmat

#endif /* PIPISTAR GUARD */
