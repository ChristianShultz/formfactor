/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : formfactor_factory.cc

 * Purpose :

 * Creation Date : 22-02-2014

 * Last Modified : Mon 17 Mar 2014 03:55:10 PM EDT

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/



#include "formfactor_factory.h"
#include "formfactor_helicity_formfactors.h"
#include "formfactor_subduced_formfactors.h"
#include "formfactor_cubic_formfactors.h"
#include <string>
#include <complex>
#include <exception>
#include "radmat/utils/printer.h"

namespace FacEnv = radmat::FormFactorDecompositionFactoryEnv;
typedef radmat::TheFormFactorInjectionRecipeFactory Factory;


namespace radmat
{

  namespace FormFactorDecompositionFactoryEnv
  {


    namespace
    {
      void splash_keys(void) 
      {
        std::vector<std::string> keys = FacEnv::all_keys(); 
        std::vector<std::string>::const_iterator k; 
        std::cout << __PRETTY_FUNCTION__ << ": available keys" << std::endl;
        for(k = keys.begin(); k != keys.end(); ++k)
          std::cout << *k << std::endl;

      }

      bool registered = false;
    } // anonomyous 

    // make it blow up if anything goes wrong by wrapping another 
    // call around the factory.createObj method
    rHandle<FormFactorBase_t > 
      callFactory(const std::string &matElemID)
      {
        FormFactorBase_t *foo;
        foo = NULL;
        try
        {
          foo = Factory::Instance().createObject(matElemID,matElemID);
        }
        catch(std::exception &e)
        {
          std::cout << "tried to find - " << matElemID << std::endl;
          std::cout << __PRETTY_FUNCTION__ << e.what(); 
          splash_keys();
          throw e; 
        }
        catch(std::string &s)
        {
          std::cout << "tried to find - " << matElemID << std::endl;
          std::cout << __PRETTY_FUNCTION__ << s << std::endl;
          splash_keys();
          throw s;
        }
        catch(...)
        {
          std::cout << "tried to find  - " << matElemID << std::endl;
          std::cout << __PRETTY_FUNCTION__ << ": some error" << std::endl;
          splash_keys();
          POW2_ASSERT(false); 
        }

        // not a null pointer
        POW2_ASSERT(foo);
        return rHandle<FormFactorBase_t>(foo); 
      }

    // dump all keys in the factory
    std::vector<std::string> 
      all_keys(void)
      {
        return Factory::Instance().keys(); 
      }


    // register the factory "inventory"
    bool registerAll( void )
    {
      printer_function<debug_print_reg_all>(__PRETTY_FUNCTION__); 
      bool success = true;

      if(!!!registered)
      {
        success &= ::radmat::HelicityFormFactorDecompositionFactoryEnv::registerAll(); 
        success &= ::radmat::SubducedFormFactorDecompositionFactoryEnv::registerAll(); 
        success &= ::radmat::CubicFormFactorDecompositionFactoryEnv::registerAll(); 
        registered = true;
      }

      if( !!! success )
      {
        throw std::string("failed to reg in FormFactorDecompositionFactoryEnv"); 
      }

      return success;
    }

  } // close FormFactorDecompositionFactoryEnv

} // close radmat
