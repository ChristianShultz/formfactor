/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : formfactor_subduced_formfactors.cc

 * Purpose :

 * Creation Date : 13-03-2014

 * Last Modified : Fri 03 Oct 2014 04:36:09 PM EDT

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/

#include "formfactor_subduced_formfactors.h"
#include "formfactor.h"
#include "formfactor_helicity_formfactors.h"
#include "formfactor_factory.h"
#include "radmat/utils/stringify.h"
#include "radmat/utils/pow2assert.h"
#include "radmat/utils/printer.h"
#include "radmat/data_representation/data_representation.h"
#include "semble/semble_meta.h"
#include <sstream>
#include <algorithm>


typedef radmat::TheFormFactorInjectionRecipeCookbook Cookbook; 
typedef radmat::TheFormFactorInjectionRecipeFactory Factory; 



namespace radmat
{

  namespace 
  {

    struct subduce_table_printer
    {
      static void print(const std::string &msg)
      {}
      // { std::cout << "subduce_table_printer " + msg << std::endl; }
    };

    struct subduce_names_printer
    {
      static void print(const std::string &msg)
      {}
      // { std::cout << __PRETTY_FUNCTION__ <<  msg << std::endl; }
    };


    const SubduceTableMap::irrep_sub_table* 
      query_subduce_table(const std::string &table_id)
      {
        SubduceTableMap::map_t::const_iterator it; 
        it = TheSmarterSubduceTableMap::Instance().mappy.find(table_id);
        if( it == TheSmarterSubduceTableMap::Instance().mappy.end() )
        {
          printer_function<console_print>( "missed subduce table " + table_id ); 
          throw std::string("missed subduce table"); 
        }

        return it->second; 
      }

  } // anonomyous 

  // all the nasty comes together here 
  itpp::Mat<std::complex<double> >
    SubducedFormFactor::generate_ffs(const MomRowPair_t &lefty, 
        const MomRowPair_t &righty, 
        const double mom_fac) const 
    {
      FormFactorBase_t::recipe_h recipe_h = FormFactorBase_t::get_recipe(); 
      POW2_ASSERT( recipe_h->id() == Stringify<SubducedFormFactorRecipe_t>() ); 
      const SubducedFormFactorRecipe_t *recipe; 
      recipe = dynamic_cast<const SubducedFormFactorRecipe_t *>(recipe_h.get_ptr()); 

      // this is the return object, it is a 4 x nff matrix, zero it
      // out explicitly and then just sum on the heicities with 
      // a weight to generate a subduced decomposition 
      itpp::Mat<std::complex<double> > ret(4,recipe->nFacs());
      ret.zeros(); 

      int left_row = lefty.second; 
      int right_row = righty.second; 
      std::complex<double> cmplx_zero(0.,0.); 

      // left and right representations to make the strings to 
      // call the factory 
      const CubicRep *l_cub_rep_ptr , *r_cub_rep_ptr;
      const LorentzRep *l_sph_rep_ptr , *r_sph_rep_ptr; 

      // get some pointers, cant use handles b/c of threading 
      l_cub_rep_ptr = recipe->lefty.get_ptr();
      r_cub_rep_ptr = recipe->righty.get_ptr();
      l_sph_rep_ptr = recipe->hel.lefty.get_ptr();
      r_sph_rep_ptr = recipe->hel.righty.get_ptr(); 

      // use clone method so our handles are function local 
      // and equivalently thread local 
      rHandle<CubicRep> l_cub_rep( l_cub_rep_ptr->clone() ); 
      rHandle<CubicRep> r_cub_rep( r_cub_rep_ptr->clone() ); 
      rHandle<LorentzRep> l_sph_rep( l_sph_rep_ptr->clone() );
      rHandle<LorentzRep> r_sph_rep( r_sph_rep_ptr->clone() );

      printer_function<subduce_names_printer>( "l_cub " + l_cub_rep->rep_id() );
      printer_function<subduce_names_printer>( "r_cub " + r_cub_rep->rep_id() );
      printer_function<subduce_names_printer>( "l_sph " + l_sph_rep->rep_id() );
      printer_function<subduce_names_printer>( "r_sph " + r_sph_rep->rep_id() );


      // generate the id for the subduce table map 
      std::string l_map_id , r_map_id; 
      l_map_id = make_subduce_table_map_id( l_sph_rep , l_cub_rep );
      r_map_id = make_subduce_table_map_id( r_sph_rep , r_cub_rep ); 

      // the subduction tables -- these were filled out 
      // earlier in the registration stage 
      const SubduceTableMap::irrep_sub_table* l_table;
      const SubduceTableMap::irrep_sub_table* r_table;

      // find the table we want 
      l_table = query_subduce_table(l_map_id);
      r_table = query_subduce_table(r_map_id); 

      // got a const ptr, use a const iterator 
      SubduceTableMap::sub_list::const_iterator l,r; 

      // make the linear combinations 
      //     -- this call eventually descends all the way down to the 
      //     lorentzff_canonical classes, these classes know about how
      //     to transform under rotations so all the hard work is hidden 
      
      itpp::Mat<std::complex<double> > tmp(4,recipe->nFacs());
      std::complex<double> weight; 

      for( l = l_table->begin(left_row); l != l_table->end(left_row); ++l)
        for( r = r_table->begin(right_row); r != r_table->end(right_row); ++r)
        {
          weight =  ( std::conj(l->first) * r->first ); 
          tmp = recipe->hel.mat->operator()( std::make_pair( lefty.first , l->second),
                std::make_pair( righty.first , r->second ),
                mom_fac); 
          
          // skip baddies ( the subduce tables only have non-zero entries so this
          // line is redundant but guards against future stupidity 
          if( std::norm(weight) < 1e-6) 
            continue; 

          //   std::cout << __func__ << " left " << l->second << " " << l->first 
          //     << "\n right " << r->second << " " << r->first  
          //     <<  "\n prod weight " << weight << " res\n" << tmp << std::endl;

          ret += weight * tmp;
        }

      return ret; 
    }




  namespace SubducedFormFactorDecompositionFactoryEnv
  {

    namespace
    {
      struct helicity_recipe_printer
      {
        static void print(const std::string &s)
        {}
        //  { std::cout << "subduced form factors, found a " << s << std::endl;}
      };

      struct subduce_printer
      {
        static void print(const std::string &s)
        {}
        //  { std::cout << "found subduction table " << s << std::endl;}
      };

      struct subduce_reg_printer
      {
        static void print(const std::string &s)
        {}
      //   { std::cout << "subduced form factors, regging  " << s << std::endl;}
      };


      std::vector<std::string> all_subduce_tables()
      {
        typedef TheSmarterSubduceTableMap SMAP; 
        std::vector<std::string> ret;
        SubduceTableMap::map_t::const_iterator it; 
        for(it = SMAP::Instance().mappy.begin(); it != SMAP::Instance().mappy.end(); ++it)
          ret.push_back(it->first); 

        return ret; 
      }

      struct possible_sub_printer
      {
        static void print(const std::string &msg)
        {}
        // {std::cout << "possible subduction " << msg << std::endl;}
      };


      // handles are not thread safe and the map may be accessed concurrently breaking 
      // the reference counting, return an independent handle to circumvent stupidity  
      std::vector<std::pair<std::string,rHandle<CubicRep> > > 
        possible_subductions(const rHandle<LorentzRep> &rep)
        {
          std::vector<std::pair<std::string,rHandle<CubicRep> > > ret; 
          std::string comp = rep->rep_id();  
          std::vector<std::string> all_tables = all_subduce_tables(); 
          std::vector<std::string>::const_iterator it;

          for(it = all_tables.begin(); it != all_tables.end(); ++it)
            if( it->find(comp) != std::string::npos )
            {
              printer_function<possible_sub_printer>(comp + " matched to " + *it); 

              const SubduceTableMap::irrep_sub_table * foo; 
              foo = TheSmarterSubduceTableMap::Instance().mappy[*it]; 
              CubicRep * rep; 
              rep = foo->get_l_rep()->clone(); 

              ret.push_back(std::make_pair(comp,rHandle<CubicRep>(rep) ) ); 
            }

          return ret; 
        }

      // loop the cookbook and pull out helicity recipes for subduction 
      std::vector<rHandle<HelicityFormFactorRecipe_t> > 
        get_all_helicity_recipies()
        {
          std::vector<rHandle<HelicityFormFactorRecipe_t> > ret;
          PtrRecipeHolder::map_t::iterator r; 
          for( r = Cookbook::Instance().mappy.begin();
              r != Cookbook::Instance().mappy.end();
              ++r)
          {
            if( r->second->id() == ::radmat::Stringify<HelicityFormFactorRecipe_t>() )
            {
              printer_function<helicity_recipe_printer>(r->first); 
              HelicityFormFactorRecipe_t *ptr = dynamic_cast<HelicityFormFactorRecipe_t*>(r->second);
              ret.push_back( 
                  rHandle<HelicityFormFactorRecipe_t>( new HelicityFormFactorRecipe_t( *ptr ) ) ); 
            }
          }

          return ret; 
        }

      FormFactorRecipe_t* gen_recipe( const HelicityFormFactorRecipe_t &h_rep, 
          const rHandle<CubicRep> &lefty,
          const rHandle<CubicRep> &righty,
          const std::string &left_table, 
          const std::string &right_table)
      {
        return new SubducedFormFactorRecipe_t(h_rep,lefty,righty,left_table,right_table); 
      }

      FormFactorBase_t * callback(const std::string &recipe_id)
      {
        PtrRecipeHolder::map_t::iterator r; 
        r = Cookbook::Instance().mappy.find(recipe_id); 

        if( r == Cookbook::Instance().mappy.end())
        {
          printer_function<console_print>( "missed " + recipe_id ); 
          throw std::string("recipe error");
        }

        if( r->second->id() != Stringify<SubducedFormFactorRecipe_t>())
        {
          printer_function<console_print>( "expected a " + Stringify<SubducedFormFactorRecipe_t>());
          printer_function<console_print>( "got a " + r->second->id() ); 
          throw std::string("recipe error");  
        }

        SubducedFormFactorRecipe_t * ptr = dynamic_cast<SubducedFormFactorRecipe_t*>(r->second);
        rHandle<SubducedFormFactorRecipe_t> recipe( new SubducedFormFactorRecipe_t( *ptr ) ); 

        return new SubducedFormFactor( recipe ); 
      }

      std::string build_id(const rHandle<CubicRep> &l, const rHandle<CubicRep> &r)
      {
        return l->rep_id() + "," + r->rep_id(); 
      }

      bool do_reg( const rHandle<HelicityFormFactorRecipe_t> &h_rep, 
          const std::pair<std::string,rHandle<CubicRep> > &lefty,
          const std::pair<std::string,rHandle<CubicRep> > &righty)
      {
        std::string reg_id = h_rep->reg_id() + "__" + build_id(lefty.second,righty.second); 
        printer_function<subduce_reg_printer>(reg_id); 
        Cookbook::Instance().mappy.insert(
            std::make_pair(reg_id,
              gen_recipe( *(h_rep.get_ptr()),
                lefty.second,
                righty.second,
                lefty.first,
                righty.first))); 

        bool b = Factory::Instance().registerObject(reg_id,callback);
        if( !!! b )
          printer_function<console_print>( reg_id + " failed to register" ); 

        return b; 
      }

    } // anonomyous 


    std::string build_id(const rHandle<CubicRep> &lefty, const rHandle<CubicRep_p> &righty)
    {
      return lefty->rep_id() + "," + righty->rep_id(); 
    }

    namespace
    {
      bool local_registration = false; 
    } 


    // pump the factory 
    bool registerAll()
    {
      printer_function<debug_print_reg_all>(__PRETTY_FUNCTION__); 
      bool success = true; 

      if( !!! local_registration )
      {
        std::vector<rHandle<HelicityFormFactorRecipe_t> > h_recipies = get_all_helicity_recipies();
        std::vector<rHandle<HelicityFormFactorRecipe_t> >::const_iterator h_recipe_it; 

        for( h_recipe_it = h_recipies.begin(); h_recipe_it != h_recipies.end(); ++h_recipe_it)
        {
          std::vector<std::pair<std::string,rHandle<CubicRep> > > lefty, righty;
          std::vector<std::pair<std::string,rHandle<CubicRep> > >::const_iterator l, r;

          lefty = possible_subductions( (*h_recipe_it)->lefty );
          righty = possible_subductions( (*h_recipe_it)->righty );

          for( l = lefty.begin(); l != lefty.end(); ++l)
            for(r = righty.begin(); r != righty.end(); ++r)
              success &= do_reg(*h_recipe_it , *l, *r ); 

        } // h_recipe_it 
      } // local_registration 

      return success; 
    }


  } // SubducedFormFactorDecompositionFactoryEnv

} // radmat

