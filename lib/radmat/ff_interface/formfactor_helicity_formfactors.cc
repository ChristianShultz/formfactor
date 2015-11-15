/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : formfactor_helicity_formfactors.cc

 * Purpose :

 * Creation Date : 22-02-2014

 * Last Modified : Fri 03 Oct 2014 11:06:52 AM EDT

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/


#include "formfactor.h"
#include "formfactor_helicity_formfactors.h"
#include "formfactor_factory.h"
#include "radmat/data_representation/data_representation.h"
#include "radmat/ff/lorentzff_canonical_cont_spin_formfactors.h"
#include "radmat/ff/lorentzff_formfactor_factory.h"
#include "radmat/utils/stringify.h"
#include "radmat/utils/pow2assert.h"
#include "radmat/utils/printer.h"
#include <sstream>
#include <algorithm>


typedef radmat::TheFormFactorInjectionRecipeCookbook Cookbook; 
typedef radmat::TheFormFactorInjectionRecipeFactory Factory; 

namespace radmat
{

  namespace
  {

    struct print_elem_reg
    {
      static void print(const std::string &msg)
      {}
      // {std::cout << "helicity form factors " << msg << std::endl;}
    };


    struct momrp_printer
    {
      static void print(const std::string &msg)
      {}
      // { std::cout << "momrp_printer " << msg << std::endl; }
    };

    std::string to_string(const MomRowPair_t &p)
    {
      std::stringstream ss; 
      ss << "p" << p.first[0] << " " << p.first[1]
        << " " << p.first[2] << " " << p.first[3];
      ss << "r" << p.second; 
      return ss.str(); 
    }

  }

  // use whatever is in the lorentzff_canonical_cont_spin decomp
  itpp::Mat<std::complex<double> >
    HelicityFormFactor::generate_ffs( const MomRowPair_t &lefty, 
        const MomRowPair_t &righty, 
        const double mom_fac) const
    {
      FormFactorBase_t::recipe_h recipe = FormFactorBase_t::get_recipe(); 
      POW2_ASSERT( recipe->id() == Stringify<HelicityFormFactorRecipe_t>() ); 
      const HelicityFormFactorRecipe_t * helicity_recipe;
      helicity_recipe = dynamic_cast< const HelicityFormFactorRecipe_t* >( recipe.get_ptr() );

      printer_function<momrp_printer>( "l->" + to_string(lefty) + "  r->" + to_string(righty) ); 

      return helicity_recipe->mat->operator()(lefty,righty,mom_fac); 
    }


  // now worry about dumping recipes into the factory
  namespace HelicityFormFactorDecompositionFactoryEnv
  {

    // local utility
    namespace
    { 

      // snarf (or slurp?) all spherical reps
      std::vector<rHandle<LorentzRep> > 
        get_all_cont_reps()
        {
          std::vector<std::string> keys;
          std::vector<std::string>::const_iterator it; 
          std::vector<rHandle<LorentzRep> > vals; 

          keys = ::radmat::LorentzRepresentationFactoryEnv::spher_keys(); 
          for(it = keys.begin(); it != keys.end(); ++it)
          {
            printer_function<print_elem_reg>("grabbed a " + *it ); 
            vals.push_back( 
                ::radmat::LorentzRepresentationFactoryEnv::callFactory(*it) ); 
          }

          return vals;
        }

      // reg function
      typedef std::pair<rHandle<LorentzRep> , rHandle<LorentzRep> > rep_pair;

      // build the id from the lorentz spin factory, 
      // this is a class name like J0pJ3m
      std::string build_lorentz_spin_id(const rHandle<LorentzRep> &lefty, 
          const rHandle<LorentzRep> &righty)
      {
        std::stringstream ss; 
        ss << lefty->rep_id() << righty->rep_id();
        return ss.str(); 
      }

      std::string build_lorentz_spin_id(const rep_pair &reps)
      {
        return build_lorentz_spin_id(reps.first,reps.second); 
      }

      // generate a recipe
      FormFactorRecipe_t* gen_recipe( const std::string &s, const rep_pair &reps )
      {
        rHandle<LorentzFFAbsBase_t> mat; 
        mat = radmat::LorentzffFormFactorDecompositionFactoryEnv::callFactory(s); 
        POW2_ASSERT( &*mat ); 
        return new HelicityFormFactorRecipe_t(mat,reps.first,reps.second); 
      }

      FormFactorBase_t* callback( const std::string &recipe_id )
      {
        PtrRecipeHolder::map_t::iterator r;  

        r = Cookbook::Instance().mappy.find(recipe_id); 
        bool success = true; 

        if( r == Cookbook::Instance().mappy.end() )
        {
          printer_function<console_print>( "missed " + recipe_id); 
          throw std::string("recipe_error"); 
        }

        if( r->second->id() != radmat::Stringify<HelicityFormFactorRecipe_t>())
        {
          printer_function<console_print>( "expected a " + Stringify<HelicityFormFactorRecipe_t>());
          printer_function<console_print>( "got a " + r->second->id());
          throw std::string("recipe_error"); 
        }

        HelicityFormFactorRecipe_t *ptr = dynamic_cast<HelicityFormFactorRecipe_t*>( r->second ); 
        rHandle<FormFactorRecipe_t> recipe( new HelicityFormFactorRecipe_t(*ptr) ); 

        return new HelicityFormFactor( recipe ); 
      }

      bool do_reg( const std::string &s, const rep_pair &r)
      {
        printer_function<print_elem_reg>(std::string("regged <") + s + std::string(">"));
        Cookbook::Instance().mappy.insert(std::make_pair(s, gen_recipe(s,r)) ); 
        bool b = true; 
        b = Factory::Instance().registerObject(s,callback);  
        if( !!! b ) 
          printer_function<console_print>( s + " failed to register" ); 
        return b; 
      }

    } // anonomyous



    // the factory reg id for helicity form factors -- up to the diag/tran bit
    std::string build_id( const rHandle<LorentzRep> &lefty, 
        const rHandle<LorentzRep> &righty)
    {
      std::stringstream ss; 
      ss << lefty->rep_id()  << righty->rep_id();
      return ss.str(); 
    }

    // local registration flag
    namespace 
    {
      bool registered = false; 
    }

    // pump the factory 
    bool registerAll()
    {
      printer_function<debug_print_reg_all>(__PRETTY_FUNCTION__); 
      bool success = true; 

      if ( !!! registered )
      {
        typedef std::pair<rHandle<LorentzRep> , rHandle<LorentzRep> > rep_pair;
        typedef std::map<std::string, rep_pair > map_t; 
        typedef map_t::value_type value_type; 

        map_t ff_spin_keys;
        map_t::const_iterator it; 
        std::vector<std::string> ff_allowed_spin_keys; 
        std::vector<rHandle<LorentzRep> > cont_rep_keys; 
        std::vector<rHandle<LorentzRep> >::const_iterator i,j; 

        cont_rep_keys = get_all_cont_reps();

        for( i = cont_rep_keys.begin(); i != cont_rep_keys.end(); ++i)
          for( j = cont_rep_keys.begin(); j != cont_rep_keys.end(); ++j)
          {
            std::string s = build_id(*i,*j); 
            ff_spin_keys.insert( 
                value_type( s , rep_pair(*i,*j) )); 
          }

        // snarf some more keys
        ff_allowed_spin_keys 
          = ::radmat::LorentzffFormFactorDecompositionFactoryEnv::all_keys(); 

        // if it is allowed stick it into the map that radmat uses 
        for( it = ff_spin_keys.begin(); it != ff_spin_keys.end(); ++it)
        {
          std::string tran_id = build_lorentz_spin_id( it->second ) + std::string("_tran");
          std::string diag_id = build_lorentz_spin_id( it->second ) + std::string("_diag");
          std::string test_id = build_lorentz_spin_id( it->second ) + std::string("_test");

          if( std::find( ff_allowed_spin_keys.begin(), 
                ff_allowed_spin_keys.end(), 
                tran_id) != ff_allowed_spin_keys.end() )
            success &= do_reg(tran_id, it->second); 

          if( std::find( ff_allowed_spin_keys.begin(), 
                ff_allowed_spin_keys.end(), 
                diag_id) != ff_allowed_spin_keys.end() )
            success &= do_reg(diag_id, it->second); 

          if( std::find( ff_allowed_spin_keys.begin(), 
                ff_allowed_spin_keys.end(), 
                test_id) != ff_allowed_spin_keys.end() )
            success &= do_reg(test_id, it->second); 
        }

        registered = true; 
      }

      return success; 
    }

  } // HelicityFormFactorFactoryDecompositionFactoryEnv



} // radmat


#undef CHECK_HELICITY_EXPLICIT
