/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

* File Name : data_representation_cubic_groups.cc

* Purpose :

* Creation Date : 19-03-2014

* Last Modified : Tue 08 Apr 2014 10:12:43 AM EDT

* Created By : shultz

_._._._._._._._._._._._._._._._._._._._._.*/

#include "data_representation_cubic_groups.h"
#include "data_representation_factory.h"
#include "radmat/utils/pow2assert.h"
#include "radmat/utils/printer.h"


namespace FacEnv = radmat::DataRepresentationFactoryEnv;
typedef radmat::TheDataRepresentationFactory Factory; 




namespace radmat
{

  namespace CubicRepresentationFactoryEnv
  {

    namespace
    {

      struct reg_printer
      {
        static void print(const std::string &msg)
        {}
      //  { std::cout << "cubic reps, regged " << msg << std::endl; }
      };

      struct key_printer
      {
        static void print(const std::string &msg)
        {}
//        { std::cout << "cubic reps, pulled " << msg << std::endl; }
      };

      template<class T, class U> 
        T* upCast()
        {
          T *t = new U(); 
          POW2_ASSERT(t);
          return t; 
        }

      template<typename Derived>
        bool do_reg()
        {
          Derived d; 
          bool reg = Factory::Instance().registerObject(
              d.rep_id(), upCast<Rep_p,Derived> );

          printer_function<reg_printer>(d.rep_id());

          if(!!!reg)
            printer_function<console_print>( d.rep_id() + " failed to register "); 

          return reg; 
        }


      bool do_reg_work(void)
      {
        bool success = true; 

        success &= do_reg<A1Rep_t>(); 
        success &= do_reg<A2Rep_t>(); 
        success &= do_reg<T1Rep_t>(); 
        success &= do_reg<T2Rep_t>(); 
        success &= do_reg<ERep_t >(); 

        success &= do_reg<H0D4A1Rep_t>(); 
        success &= do_reg<H0D4A2Rep_t>(); 
        success &= do_reg<H1D4E2Rep_t>(); 
        success &= do_reg<H2D4B1Rep_t>(); 
        success &= do_reg<H2D4B2Rep_t>(); 
        success &= do_reg<H3D4E2Rep_t>(); 
        success &= do_reg<H4D4A1Rep_t>(); 
        success &= do_reg<H4D4A2Rep_t>(); 

        success &= do_reg<H0D2A1Rep_t>(); 
        success &= do_reg<H0D2A2Rep_t>(); 
        success &= do_reg<H1D2B1Rep_t>(); 
        success &= do_reg<H1D2B2Rep_t>(); 
        success &= do_reg<H2D2A1Rep_t>(); 
        success &= do_reg<H2D2A2Rep_t>(); 
        success &= do_reg<H3D2B1Rep_t>(); 
        success &= do_reg<H3D2B2Rep_t>(); 
        success &= do_reg<H4D2A1Rep_t>(); 
        success &= do_reg<H4D2A2Rep_t>(); 

        success &= do_reg<H0D3A1Rep_t>(); 
        success &= do_reg<H0D3A2Rep_t>(); 
        success &= do_reg<H1D3E2Rep_t>(); 
        success &= do_reg<H2D3E2Rep_t>(); 
        success &= do_reg<H3D3A1Rep_t>(); 
        success &= do_reg<H3D3A2Rep_t>(); 
        success &= do_reg<H4D3E2Rep_t>(); 

        return success ;
      }


      bool local_registration = false; 

    } // anonomyous 


    bool registerAll(void)
    {
      printer_function<debug_print_reg_all>(__PRETTY_FUNCTION__); 

      bool success = true; 
      if( !!! local_registration)
      {
        success &= do_reg_work(); 
        local_registration = true; 
      }
      return success; 
    }


    // make it blow up if anything goes wrong by wrapping another 
    // call around the factory.createObj method
    rHandle<CubicRep> callFactory(const std::string &id)
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
      POW2_ASSERT( foo->rep_type() == ::radmat::Stringify<CubicRep_t>() ); 
      return rHandle<CubicRep>( dynamic_cast<CubicRep*>(foo) );
    }

    // dump all CubicRep keys in the factory
    std::vector<std::string> 
      cubic_keys(void)
      {
        std::vector<std::string> cubic_keys; 
        std::vector<std::string> all_keys = Factory::Instance().keys(); 
        std::vector<std::string>::const_iterator it; 

        for(it = all_keys.begin(); it != all_keys.end(); ++it )
        {
          printer_function<key_printer>(*it);
          rHandle<Rep_p> r = FacEnv::callFactory(*it); 
          if( r->rep_type() == ::radmat::Stringify<CubicRep_t>() )
            cubic_keys.push_back(*it); 
        }

        return cubic_keys; 
      }

  } // CubicRepresentationFactoryEnv

} // radmat





