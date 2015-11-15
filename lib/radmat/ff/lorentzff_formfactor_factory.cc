// lorentzff_formfactor_factory.cc -
//
// Saturday, June  2 2012
//

#include "lorentzff_formfactor_factory.h"
#include <string>
#include <complex>
#include <exception>
#include "adat/singleton.h"
#include "adat/objfactory.h"
#include "radmat/utils/printer.h"

// ffs
//#include "lorentzff_canonical_PiPi.h"
//#include "lorentzff_canonical_PiPiStar.h"
//#include "lorentzff_canonical_PiRho.h"
//#include "lorentzff_canonical_RhoPi.h"
//#include "lorentzff_canonical_RhoRho.h"


#include "lorentzff_canonical_cont_spin_formfactors.h"


#include <omp.h>

namespace FacEnv = radmat::LorentzffFormFactorDecompositionFactoryEnv;
typedef radmat::TheLorentzffFormFactorDecompositionFactory Factory;


namespace radmat
{

  namespace LorentzffFormFactorDecompositionFactoryEnv
  {

    namespace
    {
      // helper function
      template<class T, class U> 
        T* upCast(void)
        {
          T *t = new U();
          POW2_ASSERT(t);
          return t;
        }

      template<typename Derived> 
        bool 
        do_reg(void)
        {
          Derived d; 
          bool reg = Factory::Instance().registerObject(
              d.reg_id(), upCast<FFAbsBase_t,Derived> ); 

          if ( !!! reg ) 
          {
            std::cout << __PRETTY_FUNCTION__ << ": reg error for " << d.reg_id() << std::endl;
          } 
          return reg; 
        }

      bool registered = false;
    } // anonomyous 

    // make it blow up if anything goes wrong by wrapping another 
    // call around the factory.createObj method
    rHandle<FFAbsBase_t > callFactory(const std::string &matElemID)
    {
      FFAbsBase_t *foo;
      foo = NULL;
      try
      {
        foo = TheLorentzffFormFactorDecompositionFactory::Instance().createObject(matElemID);
      }
      catch(std::exception &e)
      {
        std::cout << "elem - " << matElemID << std::endl;
        std::cout << __PRETTY_FUNCTION__ << e.what(); 
        throw e; 
      }
      catch(std::string &s)
      {
        std::cout << "elem - " << matElemID << std::endl;
        std::cout << __PRETTY_FUNCTION__ << s << std::endl;
        throw s;
      }
      catch(...)
      {
        std::cout << "elem - " << matElemID << std::endl;
        std::cout << __PRETTY_FUNCTION__ << ": some error" << std::endl;
        POW2_ASSERT(false); 
      }

      // not a null pointer
      POW2_ASSERT(foo);
      return rHandle<FFAbsBase_t >(foo);
    }

    // dump all keys in the factory
    std::vector<std::string> 
      all_keys(void)
      {
        return TheLorentzffFormFactorDecompositionFactory::Instance().keys(); 
      }


    // register the factory "inventory"
    bool registerAll( void )
    {
      printer_function<debug_print_reg_all>(__PRETTY_FUNCTION__); 
      bool success = true;

      if(!!!registered)
      {

        //    success &= do_reg<radmat::PiPi<0,0> >();
        //    success &= do_reg<radmat::PiPiStar<0,0> >();
        //    success &= do_reg<radmat::PiRho<0,1> >();
        //    success &= do_reg<radmat::RhoPi<1,0> >();
        //    success &= do_reg<radmat::RhoRho<1,1> >();


        success &= do_reg<J0pJ0p_diag>();
        success &= do_reg<J0mJ0m_diag>();
        success &= do_reg<J0pJ0p_tran>();
        success &= do_reg<J0mJ0m_tran>();
        success &= do_reg<J1mJ0m_tran>();
        success &= do_reg<J1pJ0p_tran>();
        success &= do_reg<J0mJ1m_tran>();
        success &= do_reg<J0pJ1p_tran>();
        success &= do_reg<J1pJ1p_diag>();
        success &= do_reg<J1mJ1m_diag>();
        
        // test classes
        success &= do_reg<J1mJ1m_test>(); 


        registered = true;
      }

      if( !!! success )
      {
        throw std::string("failed to reg in LorentzffFormFactorDecompositionFactoryEnv"); 
      }

      return success;
    }

  } // close LorentzffFormFactorDecompositionFactoryEnv

} // close radmat
