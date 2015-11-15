#ifndef LORENTZFF_CANONICAL_RHORHO_H
#define LORENTZFF_CANONICAL_RHORHO_H 


#include "lorentzff_canonical_frame_formfacs_rotation_manager.h"
#include "lorentzff_polarization_embedding.h"
#include "radmat/utils/levi_civita.h"
#include "radmat/utils/pow2assert.h"
#include "radmat/utils/stringify.h"
#include <exception>
#include <sstream>
#include <complex>

// #define PRINT_RR_DECOMP_G1
// #define PRINT_RR_DECOMP_G2
// #define PRINT_RR_DECOMP_G3
//

#define ZERO_EXPLICIT_RHO_RHO

namespace radmat
{


  // PRD 73, 074507 (2006) 
  //

  //////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////

  struct RhoRhoG1; 
  REGISTER_STRINGIFY_TYPE( RhoRhoG1 ); 


  struct RhoRhoG1
    : public FormFacRotationManager<RhoRhoG1, std::complex<double> , 1, 1, AVERAGE_MASSES >,
    public leftSpinPTensor<1> , 
    public rightSpinPTensor<1>
  {
    typedef std::complex<double> Data_t;
    virtual ~RhoRhoG1() {}

    virtual std::string
      ff_impl() const
      {
        return "-G_1(Q^2)(p_f + p_i)^\\mu \\epsilon^*_\\nu(p_f,\\lambda_f)\\epsilon^\\nu(p_i,\\lambda_i) \\\\";
      }

    virtual Tensor<std::complex<double> , 1> 
      impl(const Tensor<double,1> &p_f, 
          const Tensor<double,1> &p_i, 
          const double mom_fac,
          const int lh,
          const int rh)  const
      {
        Tensor<std::complex<double> , 1> ret( (TensorShape<1>())[4], std::complex<double>(0.,0.)); 
        Tensor<std::complex<double> , 1> eps_left, eps_right; 
        eps_left = left_p_tensor(p_f,p_i,mom_fac,lh);
        eps_right = right_p_tensor(p_f,p_i,mom_fac,rh);
        Tensor<std::complex<double>,0> val; 
        Tensor<std::complex<double>,2> metric; 
        metric = convertTensorUnderlyingType<std::complex<double>,double,2>(g_dd()); 
        val = contract(eps_left,applyMetric(eps_right,metric,0),0,0); 

        ret = convertTensorUnderlyingType<std::complex<double>, double,1>(p_f + p_i); 

#ifdef PRINT_RR_DECOMP_G1
        std::cout << "RRG1:" << " val " << val.value() << " pp " << ret << std::endl; 
#endif

        return -( val.value() * ret ); 
      }
  };

  //////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////

  struct RhoRhoG2; 
  REGISTER_STRINGIFY_TYPE( RhoRhoG2 ); 

  struct RhoRhoG2 
    : public FormFacRotationManager<RhoRhoG2, std::complex<double> , 1 , 1,AVERAGE_MASSES>,
    public leftSpinPTensor<1> , 
    public rightSpinPTensor<1>
  {
    typedef std::complex<double> Data_t;
    virtual ~RhoRhoG2() {}
    virtual std::string 
      ff_impl() const
      {
        std::string s = "G_2(Q^2)\\left[ \\epsilon^\\mu(p_i,\\lambda_i)\\epsilon^*_\\nu(p_f,\\lambda_f)p_i^\\nu";
        s += "+ \\epsilon^{*\\mu}(p_f,\\lambda_f)\\epsilon_\\nu(p_i,\\lambda_i)p_f^\\nu \\right] \\\\";
        return s;
      }

    bool check_space(const Tensor<double,1> &lefty , const Tensor<double,1> &righty) const
    {
      return ( (lefty[1] == righty[1]) && (lefty[2] == righty[2]) && (lefty[3] == righty[3])); 
    }

    virtual Tensor<std::complex<double> , 1> 
      impl(const Tensor<double,1> &p_f, 
          const Tensor<double,1> &p_i, 
          const double mom_fac,
          const int lh, 
          const int rh)  const
      {
        Tensor<std::complex<double> , 1> ret( (TensorShape<1>())[4], std::complex<double>(0.,0.)); 
        Tensor<std::complex<double> , 1> p_left,p_right; 
        Tensor<std::complex<double> , 1> eps_left, eps_right; 
        Tensor<std::complex<double> , 0> val_a, val_b; 
        Tensor<std::complex<double> , 2> metric; 
        metric = convertTensorUnderlyingType<std::complex<double>,double,2>(g_dd()); 
        eps_left = left_p_tensor(p_f,p_i,mom_fac,lh);
        eps_right = right_p_tensor(p_f,p_i,mom_fac,rh);
        p_left = convertTensorUnderlyingType< Data_t , double , 1>(p_f); 
        p_right = convertTensorUnderlyingType< Data_t , double , 1>(p_i); 
        val_a = contract(eps_left, applyMetric(p_right,metric,0),0,0);  
        val_b = contract(eps_right, applyMetric(p_left,metric,0),0,0);  

#ifdef ZERO_EXPLICIT_RHO_RHO
        if (check_space(p_f,p_i))
        {
          val_a.value() = std::complex<double>(0.,0.); 
          val_b.value() = std::complex<double>(0.,0.); 
        }
#endif 


        ret = val_a.value() * eps_right + val_b.value() * eps_left; 

#ifdef PRINT_RR_DECOMP_G2
          std::cout << "RRG2:" << ": pleft = " << p_left 
            << " pright = " << p_right 
            << " epsl = " << eps_left 
            << " epsr = " << eps_right 
            << " el.pr = " << val_a.value() << "     er.pl = " << val_b.value() 
            << std::endl;  
          std::cout << "returning ******* \n" << ret << std::endl;
#endif 

        return ret; 
      }
  };


  //////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////

  struct RhoRhoG3; 
  REGISTER_STRINGIFY_TYPE(RhoRhoG3); 

  struct RhoRhoG3 
    : public FormFacRotationManager<RhoRhoG3, std::complex<double> , 1, 1,AVERAGE_MASSES>,
    public leftSpinPTensor<1> , 
    public rightSpinPTensor<1>
  {
    typedef std::complex<double> Data_t;
    virtual ~RhoRhoG3() {}
    virtual std::string
      ff_impl() const
      {
        std::string s =  "-\\frac{G_3(Q^2)}{2m_v^2}(p_f + p_i)^\\mu";
        s += " \\epsilon^*_\\nu(p_f,\\lambda_f)p_i^\\nu";
        s += " \\epsilon_\\alpha(p_i,\\lambda_f)p_f^\\alpha \\\\ ";
        return s; 
      }


    bool check_space(const Tensor<double,1> &lefty , const Tensor<double,1> &righty) const 
    {
      return ( (lefty[1] == righty[1]) && (lefty[2] == righty[2]) && (lefty[3] == righty[3])); 
    }

    virtual Tensor<std::complex<double> , 1> 
      impl(const Tensor<double,1> &p_f, 
          const Tensor<double,1> &p_i, 
          const double mom_fac,
          const int lh,
          const int rh)  const
      {
        Tensor<std::complex<double> , 1> ret( (TensorShape<1>())[4], std::complex<double>(0.,0.)); 
        Tensor<std::complex<double> , 1> p_left,p_right; 
        Tensor<std::complex<double> , 1> eps_left, eps_right; 
        Tensor<std::complex<double> , 0> val_a, val_b, mass; 
        Tensor<std::complex<double> , 2> metric; 
        metric = convertTensorUnderlyingType<std::complex<double>,double,2>(g_dd()); 
        eps_left = left_p_tensor(p_f,p_i,mom_fac,lh);
        eps_right = right_p_tensor(p_f,p_i,mom_fac,rh);
        p_left = convertTensorUnderlyingType< Data_t , double , 1>(p_f); 
        p_right = convertTensorUnderlyingType< Data_t , double , 1>(p_i); 
        val_a = contract(eps_left, applyMetric(p_right,metric,0),0,0);  
        val_b = contract(eps_right, applyMetric(p_left,metric,0),0,0);  
        mass = contract(p_right, applyMetric(p_right,metric,0),0,0); 

        ret = convertTensorUnderlyingType<std::complex<double>, double,1>(p_f + p_i); 

#ifdef ZERO_EXPLICIT_RHO_RHO
        if (check_space(p_f,p_i))
        {
          val_a.value() = std::complex<double>(0.,0.); 
          val_b.value() = std::complex<double>(0.,0.); 
        }
#endif 

#ifdef PRINT_RR_DECOMP_G3
        std::cout << "RRG3:" << " lh " << lh << " rh " << rh << std::endl;
        std::cout << "pl " << p_f 
          << " pr " << p_i 
          << " eps_left " << eps_left 
          << " eps_r " << eps_right 
          << " pp " << ret << std::endl; 
#endif 

        // return - val_a.value() * val_b.value() * ret;
        return  - ( (val_a.value() * val_b.value() ) / (2.*mass.value()) ) * ret; 
      }
  };




  //////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////

  template<int embedl, int embedr>
    LorentzFFAbsBase_t::LorentzFFAbs_list RhoRhoGenList()
    {
      LorentzFFAbsBase_t::LorentzFFAbs_list retRhoRho; 
      LorentzFFAbsBase_t::BBType *g1 , *g2, *g3; 

      try
      {
        g1 = new radmat::RhoRhoG1();
        g2 = new radmat::RhoRhoG2();
        g3 = new radmat::RhoRhoG3();

        // POW2_ASSERT(g1 && g2 && g3);
        POW2_ASSERT( g1 );
        POW2_ASSERT( g2 );
        POW2_ASSERT( g3 );

        retRhoRho.push_back(LorentzFFAbsBase_t::BBHandle_t(g1)); 
        retRhoRho.push_back(LorentzFFAbsBase_t::BBHandle_t(g2)); 
        retRhoRho.push_back(LorentzFFAbsBase_t::BBHandle_t(g3)); 
      }
      catch(...)
      {
        POW2_ASSERT(false); 
      } 

      return retRhoRho;
    }



  template<int embedl,int embedr> struct RhoRho; 
  REGISTER_STRINGIFY_TYPE2( RhoRho<1,1> ); 


  //////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////
  template<int embedl, int embedr>
    struct RhoRho : public LorentzFFAbsBase_t
  {
    RhoRho()
      : LorentzFFAbsBase_t(radmat::RhoRhoGenList<embedl,embedr>())
    { }

    RhoRho& operator=(const RhoRho &o)
    {
      if (this != &o)
        LorentzFFAbsBase_t::operator=(o);

      return *this; 
    }

    RhoRho(const RhoRho &o)
      : LorentzFFAbsBase_t(o)
    {  }

    virtual ~RhoRho() {}

    virtual std::string reg_id() const { return Stringify< RhoRho<embedl,embedr> >(); }
    virtual int left_spin() const { return embedl; }
    virtual int right_spin() const { return embedr; }
    virtual LorentzFFAbsBase_t * clone() const { return new RhoRho(); }

    private: 
    RhoRho(const LorentzFFAbsBase_t::LorentzFFAbs_list &); 
    RhoRho(const LorentzFFAbsBase_t::LorentzFFAbs_list ); 
  };

  //
  // MULTIPOLE BASIS
  //

  //////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////

  struct RhoRhoEtaGen
  {
    double operator()(const Tensor<double,1> &l, const Tensor<double,1> &r)
    {
      Tensor<double,2 > metric = g_dd(); 
      Tensor<double,0 > Q2 = contract( l-r, applyMetric(l-r,metric,0) , 0, 0); 
      Tensor<double,0 > mass = contract(r, applyMetric(l,metric,0) ,0,0); 

      //  Tensor<double,0 > ml = contract(l, applyMetric(l,metric,0) ,0,0); 
      //  Tensor<double,0 > mr = contract(r, applyMetric(r,metric,0) ,0,0); 
      //  // average over irrep splitting ( Q^2 / 4 m^2 ) 
      //  return -Q2.value() / ( 2. *( ml.value() + mr.value() ) );  

      return -Q2.value() / ( 4. * mass.value() ) ;  
    }
  }; 

  //////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////


  struct RhoRhoGC; 
  REGISTER_STRINGIFY_TYPE( RhoRhoGC ); 

  struct RhoRhoGC
    : public FFAbsBlockBase_t<std::complex<double> >
  {
    typedef std::complex<double> Data_t;
    virtual ~RhoRhoGC() {}

    virtual std::string ff() const {return "";}

    virtual std::string id() const {return Stringify<RhoRhoGC>(); }

    virtual Tensor<Data_t,1> 
      operator()(const MomRowPair_t &lefty, 
          const MomRowPair_t &righty, 
          const double mom_fac) const
      {
        RhoRhoEtaGen eg; 
        double eta = eg(lefty.first,righty.first); 
        double two_thirds_eta = eta*2./3.;

        std::complex<double> g1coeff = std::complex<double>( 1. + two_thirds_eta, 0.); 
        std::complex<double> g2coeff = std::complex<double>( -two_thirds_eta, 0.); 
        std::complex<double> g3coeff = std::complex<double>( two_thirds_eta*( 1. + eta), 0.); 

        RhoRhoG1 G1; 
        RhoRhoG2 G2; 
        RhoRhoG3 G3; 

        return g1coeff * G1(lefty,righty,mom_fac) 
          + g2coeff * G2(lefty,righty,mom_fac)
          + g3coeff * G3(lefty,righty,mom_fac); 
      }

  };



  //////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////


  struct RhoRhoGM; 
  REGISTER_STRINGIFY_TYPE( RhoRhoGM ); 

  struct RhoRhoGM
    : public FFAbsBlockBase_t<std::complex<double> >
  {
    typedef std::complex<double> Data_t;
    virtual ~RhoRhoGM() {}

    virtual std::string ff() const {return "";}

    virtual std::string id() const {return Stringify<RhoRhoGM>(); }

    virtual Tensor<Data_t,1> 
      operator()(const MomRowPair_t &lefty, 
          const MomRowPair_t &righty, 
          const double mom_fac) const
      {
        std::complex<double> g2coeff = std::complex<double>( 1., 0.); 

        RhoRhoG2 G2; 

        return g2coeff * G2(lefty,righty,mom_fac);
      }

  };


  //////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////


  struct RhoRhoGQ; 
  REGISTER_STRINGIFY_TYPE( RhoRhoGQ ); 

  struct RhoRhoGQ
    : public FFAbsBlockBase_t<std::complex<double> >
  {
    typedef std::complex<double> Data_t;
    virtual ~RhoRhoGQ() {}

    virtual std::string ff() const {return "";}

    virtual std::string id() const {return Stringify<RhoRhoGQ>(); }

    virtual Tensor<Data_t,1> 
      operator()(const MomRowPair_t &lefty, 
          const MomRowPair_t &righty, 
          const double mom_fac) const
      {
        RhoRhoEtaGen eg; 
        double eta = eg(lefty.first,righty.first); 

        std::complex<double> g1coeff = std::complex<double>( 1. , 0.); 
        std::complex<double> g2coeff = std::complex<double>( -1., 0.); 
        std::complex<double> g3coeff = std::complex<double>(  1. + eta, 0.); 

        RhoRhoG1 G1; 
        RhoRhoG2 G2; 
        RhoRhoG3 G3; 

        return g1coeff * G1(lefty,righty,mom_fac) 
          + g2coeff * G2(lefty,righty,mom_fac)
          + g3coeff * G3(lefty,righty,mom_fac); 
      }

  };

  //////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////

  template<int embedl, int embedr>
    LorentzFFAbsBase_t::LorentzFFAbs_list RhoRhoGenListMultipole()
    {
      LorentzFFAbsBase_t::LorentzFFAbs_list retRhoRho; 
      LorentzFFAbsBase_t::BBType *gc , *gm, *gq; 

      try
      {
        gc = new radmat::RhoRhoGC();
        gm = new radmat::RhoRhoGM();
        gq = new radmat::RhoRhoGQ();

        POW2_ASSERT( gc );
        POW2_ASSERT( gm );
        POW2_ASSERT( gq );

        retRhoRho.push_back(LorentzFFAbsBase_t::BBHandle_t(gc)); 
        retRhoRho.push_back(LorentzFFAbsBase_t::BBHandle_t(gm)); 
        retRhoRho.push_back(LorentzFFAbsBase_t::BBHandle_t(gq)); 
      }
      catch(...)
      {
        POW2_ASSERT(false); 
      } 

      return retRhoRho;
    }



  template<int embedl,int embedr> struct RhoRhoMultipole; 
  REGISTER_STRINGIFY_TYPE2( RhoRhoMultipole<1,1> ); 


  //////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////
  template<int embedl, int embedr>
    struct RhoRhoMultipole : public LorentzFFAbsBase_t
  {
    RhoRhoMultipole()
      : LorentzFFAbsBase_t(radmat::RhoRhoGenListMultipole<embedl,embedr>())
    { }

    RhoRhoMultipole& operator=(const RhoRhoMultipole &o)
    {
      if (this != &o)
        LorentzFFAbsBase_t::operator=(o);

      return *this; 
    }

    RhoRhoMultipole(const RhoRhoMultipole &o)
      : LorentzFFAbsBase_t(o)
    {  }

    virtual ~RhoRhoMultipole() {}

    virtual std::string reg_id() const { return Stringify< RhoRhoMultipole<embedl,embedr> >(); }
    virtual int left_spin() const { return embedl; }
    virtual int right_spin() const { return embedr; }
    virtual LorentzFFAbsBase_t * clone() const { return new RhoRhoMultipole(); }

    private: 
    RhoRhoMultipole(const LorentzFFAbsBase_t::LorentzFFAbs_list &); 
    RhoRhoMultipole(const LorentzFFAbsBase_t::LorentzFFAbs_list ); 
  };


} // namespace radmat



#endif /* LORENTZFF_CANONICAL_RHORHO_H */
