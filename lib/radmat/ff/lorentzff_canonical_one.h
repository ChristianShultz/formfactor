#ifndef LORENTZFF_CANONICAL_ONE_H
#define LORENTZFF_CANONICAL_ONE_H 



#include "lorentzff_polarization_embedding.h"
#include "radmat/utils/pow2assert.h"
#include "lorentzff_canonical_frame_formfacs_rotation_manager.h"
#include "radmat/utils/stringify.h"
#include <complex>

namespace radmat
{

  template<int Jl,int Jr>
  struct ONEF1; 

  REGISTER_STRINGIFY_TYPE2( ONEF1<0,0> ); 
  REGISTER_STRINGIFY_TYPE2( ONEF1<1,0> ); 
  REGISTER_STRINGIFY_TYPE2( ONEF1<2,0> ); 
  REGISTER_STRINGIFY_TYPE2( ONEF1<3,0> ); 
  REGISTER_STRINGIFY_TYPE2( ONEF1<0,1> ); 
  REGISTER_STRINGIFY_TYPE2( ONEF1<1,1> ); 
  REGISTER_STRINGIFY_TYPE2( ONEF1<2,1> ); 
  REGISTER_STRINGIFY_TYPE2( ONEF1<3,1> ); 
  REGISTER_STRINGIFY_TYPE2( ONEF1<0,2> ); 
  REGISTER_STRINGIFY_TYPE2( ONEF1<1,2> ); 
  REGISTER_STRINGIFY_TYPE2( ONEF1<2,2> ); 
  REGISTER_STRINGIFY_TYPE2( ONEF1<3,2> ); 
  REGISTER_STRINGIFY_TYPE2( ONEF1<0,3> ); 
  REGISTER_STRINGIFY_TYPE2( ONEF1<1,3> ); 
  REGISTER_STRINGIFY_TYPE2( ONEF1<2,3> ); 
  REGISTER_STRINGIFY_TYPE2( ONEF1<3,3> ); 

  template<int Jl , int Jr >
  struct ONEF1
    : public FormFacRotationManager<ONEF1<Jl,Jr>, std::complex<double>, Jl , Jr , AVERAGE_MASSES>
  {
    virtual ~ONEF1() {}

    virtual  std::string ff_impl() const
    {
      return std::string(" <1,1,1,1>");
    }

    // return a complex version of p_+
    virtual Tensor<std::complex<double> , 1> 
      impl(const Tensor<double,1> &p_f,
          const Tensor<double,1> &p_i,
          const double mom_fac,
          int h_f,
          int h_i) const
      {
       return Tensor<std::complex<double>,1>((TensorShape<1>())[4],std::complex<double>(1.,0.)); 
      }
  };



  ////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////


  // generate a list for the ONE constructor
  template<int embedl, int embedr>
    LorentzFFAbsBase_t::LorentzFFAbs_list ONEGenList()
    {
      LorentzFFAbsBase_t::LorentzFFAbs_list retONE;
      LorentzFFAbsBase_t::BBType *blockPtr;
      blockPtr = new ONEF1<embedl,embedr>();
      POW2_ASSERT(blockPtr);  // blow up if something went wrong
      retONE.push_back(LorentzFFAbsBase_t::BBHandle_t(blockPtr));
      return retONE;
    }


  ////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////

  // can only instantiate the <i,j> types now
  template<int embedl, int embedr> struct ONE; 
  REGISTER_STRINGIFY_TYPE2( ONE<0,0> ); 
  REGISTER_STRINGIFY_TYPE2( ONE<1,0> ); 
  REGISTER_STRINGIFY_TYPE2( ONE<2,0> ); 
  REGISTER_STRINGIFY_TYPE2( ONE<3,0> ); 
  REGISTER_STRINGIFY_TYPE2( ONE<0,1> ); 
  REGISTER_STRINGIFY_TYPE2( ONE<1,1> ); 
  REGISTER_STRINGIFY_TYPE2( ONE<2,1> ); 
  REGISTER_STRINGIFY_TYPE2( ONE<3,1> ); 
  REGISTER_STRINGIFY_TYPE2( ONE<0,2> ); 
  REGISTER_STRINGIFY_TYPE2( ONE<1,2> ); 
  REGISTER_STRINGIFY_TYPE2( ONE<2,2> ); 
  REGISTER_STRINGIFY_TYPE2( ONE<3,2> ); 
  REGISTER_STRINGIFY_TYPE2( ONE<0,3> ); 
  REGISTER_STRINGIFY_TYPE2( ONE<1,3> ); 
  REGISTER_STRINGIFY_TYPE2( ONE<2,3> ); 
  REGISTER_STRINGIFY_TYPE2( ONE<3,3> ); 



  // only need to derive the constructor, everything else is in the 
  // base class, do this polymorphically, make some function that 
  // returns the appropriate handle based on the requested matrix element type
  template<int embedl, int embedr>
    struct ONE : public LorentzFFAbsBase_t
  {
    ONE()
      : LorentzFFAbsBase_t(radmat::ONEGenList<embedl,embedr>())  
    {  } 

    ONE& operator=(const ONE &o)
    {

      if(this != &o)
        LorentzFFAbsBase_t::operator=(o);
      return *this;
    }

    // no slicing
    ONE(const ONE &o)
      : LorentzFFAbsBase_t(o)
    {  }

    virtual ~ONE() {} 

    virtual std::string reg_id() const { return Stringify< ONE<embedl,embedr> >(); }
    virtual int left_spin() const {return embedl;}
    virtual int right_spin() const {return embedr;}

    virtual LorentzFFAbsBase_t* clone() const { return new ONE(); }; 


    private:
    // I'm not sure if these could inherit so we will hide them as well
    ONE(const LorentzFFAbsBase_t::LorentzFFAbs_list &);
    ONE(const LorentzFFAbsBase_t::LorentzFFAbs_list);
  };


} // close radmat


#endif /* LORENTZFF_CANONICAL_ONE_H */
