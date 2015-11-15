#ifndef LORENTZFF_CANONICAL_PIRHO_H
#define LORENTZFF_CANONICAL_PIRHO_H 


#include "lorentzff_polarization_embedding.h"
#include "radmat/utils/levi_civita.h"
#include "radmat/utils/pow2assert.h"
#include "lorentzff_canonical_frame_formfacs_rotation_manager.h"


namespace radmat
{

  struct PiRhoF1; 
  REGISTER_STRINGIFY_TYPE( PiRhoF1 ); 


  // actual implementation of the kinematic factor
  struct PiRhoF1
    : public FormFacRotationManager<PiRhoF1, std::complex<double> , 0 , 1 , DONT_AVERAGE_MASSES> , 
    public rightSpinPTensor<1>
  {
    typedef std::complex<double> Data_t; 
    virtual ~PiRhoF1() {}

    virtual std::string
      ff_impl() const
      {
        std::string s; 
        s = "F_1(Q^2) \\epsilon^{\\mu,\\nu,\\rho,\\sigma}\\epsilon_{\\nu}";
        s += "(p,\\lambda)p_{\\rho}^{+}p_{\\sigma}^{-}";
        return s; 
      }

    virtual Tensor<std::complex<double> , 1> 
      impl(const Tensor<double,1> &p_f, 
          const Tensor<double,1> &p_i, 
          const double mom_fac,
          int pihel,
          int rhohel) const
      {
        // come up with the ingredient list
        Tensor<std::complex<double>, 1> epsilon = this->right_p_tensor(p_f,p_i,mom_fac,rhohel); 
        Tensor<std::complex<double>, 1> pleft , pright;
        pleft = convertTensorUnderlyingType<std::complex<double>,double,1>( p_f );
        pright = convertTensorUnderlyingType<std::complex<double>,double,1>( p_i );
        Tensor<std::complex<double>,4> levi = levi_civita<std::complex<double>,4>(); ; 
        Tensor<std::complex<double>, 2> gdd;
        gdd = convertTensorUnderlyingType<std::complex<double>,double,2>(g_dd());

        pleft = applyMetric(pleft,gdd,0); 
        pright = applyMetric(pright,gdd,0); 
        epsilon = applyMetric(epsilon,gdd,0); 


        Tensor<double, 0> m_left, m_right;
        std::complex<double> norm; 
        m_left = contract(p_f,applyMetric(p_f,g_dd(),0),0,0);
        m_right = contract(p_i,applyMetric(p_i,g_dd(),0),0,0);
        norm = std::complex<double>( 2./( sqrt(m_left.value()) + sqrt(m_right.value()) ), 0.); 

#if 0
        Tensor<std::complex<double> , 0> inner_prod = contract( epsilon, p_i , 0 , 0 ) ; 
        if ( std::norm ( inner_prod.value() ) >  1e-6 ) 
          std::cout << "mom dotted into polarization was " << inner_prod.value() << std::endl; 
#endif 

        return  norm * contract(
            contract(
              contract(levi,
                pleft , 3 , 0),
              pright , 2 , 0 ),
            epsilon , 1 , 0 );
      }
  };


  // generate a list for the PiPi constructor
  //
  //  use an embedding so we can play with subduction later
  //
  template< int embedl , int embedr  > 
    LorentzFFAbsBase_t::LorentzFFAbs_list PiRhoGenList()
    {
      LorentzFFAbsBase_t::LorentzFFAbs_list retCanonicalPiRho;
      LorentzFFAbsBase_t::BBType *blockPtr;
      blockPtr = new PiRhoF1();
      POW2_ASSERT(blockPtr);
      retCanonicalPiRho.push_back(LorentzFFAbsBase_t::BBHandle_t(blockPtr));
      return retCanonicalPiRho;
    }


  ////////////////////////////
  ////////////////////////////

  template<int embedl, int embedr> struct PiRho; 
  REGISTER_STRINGIFY_TYPE2( PiRho<0,1> ); 


  //  use an embedding so we can play with subduction later
  template<int embedl, int embedr>
    struct PiRho : public LorentzFFAbsBase_t
  {
    PiRho()
      : LorentzFFAbsBase_t(radmat::PiRhoGenList<embedl,embedr>())
    {  }

    PiRho& operator=(const PiRho &o)
    {
      if(this != &o)
        LorentzFFAbsBase_t::operator=(o);
      return *this;
    }

    PiRho(const PiRho &o)
      : LorentzFFAbsBase_t(o)
    { }

    virtual ~PiRho() {}

    virtual std::string reg_id() const { return Stringify<PiRho<embedl,embedr> >(); }
    virtual int left_spin() const { return embedl; }
    virtual int right_spin() const { return embedr; }
    virtual LorentzFFAbsBase_t* clone() const { return new PiRho(); }

    private:
    PiRho(const LorentzFFAbsBase_t::LorentzFFAbs_list &);
    PiRho(const LorentzFFAbsBase_t::LorentzFFAbs_list); 

  };

}








#endif /* LORENTZFF_CANONICAL_PIRHO_H */
