/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

* File Name : data_representation_factory.cc

* Purpose :

* Creation Date : 19-03-2014

* Last Modified : Thu 20 Mar 2014 09:16:24 AM EDT

* Created By : shultz

_._._._._._._._._._._._._._._._._._._._._.*/


#include "data_representation_factory.h"
#include "data_representation_cubic_groups.h"
#include "data_representation_lorentz_groups.h"
#include "radmat/utils/pow2assert.h"
#include "radmat/utils/printer.h"



namespace FacEnv = radmat::DataRepresentationFactoryEnv;
typedef radmat::TheDataRepresentationFactory Factory;


namespace radmat
{

  namespace DataRepresentationFactoryEnv
  {
    namespace 
    {
      bool local_registration = false; 
    }

    bool registerAll(void)
    {
      printer_function<debug_print_reg_all>(__PRETTY_FUNCTION__); 
      bool success = true; 
      if( !!! local_registration )
      {
        success &= LorentzRepresentationFactoryEnv::registerAll(); 
        success &= CubicRepresentationFactoryEnv::registerAll(); 
        local_registration = true; 
      }

      return success; 
    }

    // make it blow up if anything goes wrong by wrapping another 
    // call around the factory.createObj method
    rHandle<Rep_p > callFactory(const std::string &id)
    {
      Rep_p *foo;
      foo = NULL;
      try
      {
        foo = Factory::Instance().createObject(id);
      }
      catch(std::exception &e)
      {
        std::cout << "id - " << id << std::endl;
        std::cout << __PRETTY_FUNCTION__ << e.what(); 
        throw e; 
      }
      catch(std::string &s)
      {
        std::cout << "id - " << id << std::endl;
        std::cout << __PRETTY_FUNCTION__ << s << std::endl;
        throw s;
      }
      catch(...)
      {
        std::cout << "id - " << id << std::endl;
        std::cout << __PRETTY_FUNCTION__ << ": some error" << std::endl;
        POW2_ASSERT(false); 
      }

      // not a null pointer
      POW2_ASSERT(foo);
      return rHandle<Rep_p >(foo);
    }


  }// Fac Env

}

