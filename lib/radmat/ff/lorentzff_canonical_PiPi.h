#ifndef LORENTZFF_PIPI_H_H_GUARD
#define LORENTZFF_PIPI_H_H_GUARD


#include "lorentzff_polarization_embedding.h"
#include "radmat/utils/pow2assert.h"
#include "lorentzff_canonical_frame_formfacs_rotation_manager.h"
#include "radmat/utils/stringify.h"
#include <complex>

namespace radmat
{

  struct PiPiF1; 
  REGISTER_STRINGIFY_TYPE( PiPiF1 ); 

  // only one ff
  struct PiPiF1
    : public FormFacRotationManager<PiPiF1, std::complex<double>, 0 , 0 , AVERAGE_MASSES>
  {
    virtual ~PiPiF1() {}

    virtual  std::string ff_impl() const
    {
      return std::string(" F_1(Q^2) p_+^{\\mu} ");
    }

    // return a complex version of p_+
    virtual Tensor<std::complex<double> , 1> 
      impl(const Tensor<double,1> &p_f,
          const Tensor<double,1> &p_i,
          const double mom_fac,
          int h_f,
          int h_i) const
      {
        return convertTensorUnderlyingType<std::complex<double>,double,1>( pPlus(p_f,p_i) );
      }
  };



  ////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////


  // generate a list for the PiPi constructor
  template<int embedl, int embedr>
    LorentzFFAbsBase_t::LorentzFFAbs_list PiPiGenList()
    {
      LorentzFFAbsBase_t::LorentzFFAbs_list retPiPi;
      LorentzFFAbsBase_t::BBType *blockPtr;
      blockPtr = new PiPiF1();
      POW2_ASSERT(blockPtr);  // blow up if something went wrong
      retPiPi.push_back(LorentzFFAbsBase_t::BBHandle_t(blockPtr));
      return retPiPi;
    }


  ////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////

  // can only instantiate the <0,0> type now
  template<int embedl, int embedr> struct PiPi; 
  REGISTER_STRINGIFY_TYPE2( PiPi<0,0> ); 


  // only need to derive the constructor, everything else is in the 
  // base class, do this polymorphically, make some function that 
  // returns the appropriate handle based on the requested matrix element type
  template<int embedl, int embedr>
    struct PiPi : public LorentzFFAbsBase_t
  {
    PiPi()
      : LorentzFFAbsBase_t(radmat::PiPiGenList<embedl,embedr>())  
    {  } 

    PiPi& operator=(const PiPi &o)
    {

      if(this != &o)
        LorentzFFAbsBase_t::operator=(o);
      return *this;
    }

    // no slicing
    PiPi(const PiPi &o)
      : LorentzFFAbsBase_t(o)
    {  }

    virtual ~PiPi() {} 

    virtual std::string reg_id() const { return Stringify< PiPi<embedl,embedr> >(); }
    virtual int left_spin() const {return embedl;}
    virtual int right_spin() const {return embedr;}

    virtual LorentzFFAbsBase_t* clone() const { return new PiPi(); }; 


    private:
    // I'm not sure if these could inherit so we will hide them as well
    PiPi(const LorentzFFAbsBase_t::LorentzFFAbs_list &);
    PiPi(const LorentzFFAbsBase_t::LorentzFFAbs_list);
  };


} // close radmat

#endif
