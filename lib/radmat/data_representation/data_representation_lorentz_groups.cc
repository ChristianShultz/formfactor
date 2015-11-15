/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

* File Name : data_representation_lorentz_groups.cc

* Purpose :

* Creation Date : 19-03-2014

* Last Modified : Tue 08 Apr 2014 10:12:15 AM EDT

* Created By : shultz

_._._._._._._._._._._._._._._._._._._._._.*/

#include "data_representation_factory.h"
#include "data_representation_lorentz_groups.h"
#include "radmat/utils/pow2assert.h"
#include "radmat/utils/printer.h"

namespace FacEnv = radmat::DataRepresentationFactoryEnv;
typedef radmat::TheDataRepresentationFactory Factory;


namespace radmat
{

  namespace LorentzRepresentationFactoryEnv
  {
    

    namespace
    {

      struct reg_printer
      {
        static void print(const std::string &msg)
        {}
        // { std::cout << "spher rep, regged " << msg << std::endl; }
      };

      struct key_printer
      {
        static void print(const std::string &msg)
        {}
        // {std::cout << __PRETTY_FUNCTION__<<  "  " << msg << std::endl;}
      };

      bool registered = false; 

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
              d.rep_id() , upCast<Rep_p,Derived> ); 

          printer_function<reg_printer>(d.rep_id()); 

          if( !!! reg )
            std::cout << __PRETTY_FUNCTION__ 
              << ": reg error for " << d.rep_id() << std::endl;

          return reg; 
        }

    } // anonomyous 

    bool registerAll(void)
    {
      printer_function<debug_print_reg_all>(__PRETTY_FUNCTION__); 
      bool success = true; 

      if( !!! registered )
      {
        success &= do_reg<J0pRep_t>(); 
        success &= do_reg<J0mRep_t>();

        success &= do_reg<J1pRep_t>(); 
        success &= do_reg<J1mRep_t>();

        success &= do_reg<J2pRep_t>(); 
        success &= do_reg<J2mRep_t>();

        success &= do_reg<J3pRep_t>(); 
        success &= do_reg<J3mRep_t>();

        success &= do_reg<J4pRep_t>(); 
        success &= do_reg<J4mRep_t>();

        registered = true; 
      }

      return success; 
    }


    // make it blow up if anything goes wrong by wrapping another 
    // call around the factory.createObj method
    rHandle<LorentzRep> callFactory(const std::string &id)
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
      POW2_ASSERT( foo->rep_type() == ::radmat::Stringify<LorentzRep_t>() ); 
      return rHandle<LorentzRep>( dynamic_cast<LorentzRep*>(foo) );
    }

    // dump all LorentzRep_p keys in the factory
    std::vector<std::string> 
      spher_keys(void)
      {
        std::vector<std::string> sph_keys; 
        std::vector<std::string> all_keys = Factory::Instance().keys(); 
        std::vector<std::string>::const_iterator it; 

        for(it = all_keys.begin(); it != all_keys.end(); ++it )
        {
          printer_function<key_printer>(*it);
          rHandle<Rep_p> r = DataRepresentationFactoryEnv::callFactory(*it); 
          if( r->rep_type() == ::radmat::Stringify<LorentzRep_t>() )
            sph_keys.push_back(*it); 
        }

        return sph_keys; 
      }



  } // LorentzRepresentationFactoryEnv

} // radmat


