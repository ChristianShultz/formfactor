/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : data_representation_subduction_map.cc

 * Purpose :

 * Creation Date : 19-03-2014

 * Last Modified : Mon 29 Sep 2014 06:59:05 PM EDT

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/


#include "data_representation_subduction_map.h"
#include "data_representation_lorentz_groups.h"
#include "data_representation_cubic_groups.h"
#include "radmat/utils/stringify.h"
#include "radmat/utils/pow2assert.h"
#include "radmat/utils/printer.h"
#include "hadron/subduce_tables_factory.h"
#include "semble/semble_meta.h"
#include <sstream>
#include <algorithm>


namespace radmat
{

  namespace
  {
    struct subduce_table_printer
    {
      static void print(const std::string &msg)
      {}
      // { std::cout << "subduce_table_printer " << msg << std::endl;}
    };

    struct subduce_printer
    {
      static void print(const std::string &msg)
       {}
       //{ std::cout << "subduce_printer " << msg << std::endl;}
    };

    SubduceTableMap::irrep_sub_table*
      subduce_rest(const std::string &table_id, 
          const rHandle<LorentzRep> &cont, 
          const rHandle<CubicRep> &cub)
      {
        SubduceTableMap::irrep_sub_table::sub_table tab; 
        Hadron::SubduceTable* sub = Hadron::TheSubduceTableFactory::Instance().createObject(table_id); 
        std::complex<double> zero(0.,0.);

        tab.resize( cub->dim() ); 

        POW2_ASSERT( &*sub); 

        for(int row = 0; row < cub->dim(); ++row )
        {
          SubduceTableMap::sub_list entry; 
          for(int lambda = 0; lambda < 2*cont->rep_spin() + 1; ++lambda)
          {

         //  std::cout << __func__ << ": attempting " + table_id + " row " 
         //   << row << " hel " << lambda << std::endl; 


            // subduce tables are FORTRAN based 
            
            ENSEM::Complex foo = (*sub)(row +1, lambda +1); 
            std::complex<double> bar = SEMBLE::toScalar(foo); 

            std::complex<double> weight = SEMBLE::toScalar( (*sub)(row+1,lambda+1) ); 
            if(weight == zero)
              continue; 

            entry.push_back(std::make_pair(weight,cont->rep_spin() - lambda) ); 
          }
          tab[row] = entry; 
        }

        delete sub; 

        SubduceTableMap::irrep_sub_table* ret = new SubduceTableMap::irrep_sub_table(tab,cont,cub); 
        return ret; 
      }

    SubduceTableMap::irrep_sub_table* 
      subduce_flight(const std::string &table_id, 
          const rHandle<LorentzRep> &cont, 
          const rHandle<CubicRep> &cub)
      { 
        SubduceTableMap::irrep_sub_table::sub_table tab; 
        Hadron::SubduceTable* sub = Hadron::TheSubduceTableFactory::Instance().createObject(table_id); 
        const CubicRepFlight_p * cubb;
        cubb = dynamic_cast< const CubicRepFlight_p *> ( cub.get_ptr() ); 
        std::complex<double> zero(0.,0.);
        double eta_p = double(cont->rep_eta_p()); 
        tab.resize( cub->dim() ); 
        int helicity = cubb->helicity();

        POW2_ASSERT( &*sub); 

        if ( helicity == 0 )
        {
          SubduceTableMap::sub_list entry;
          std::complex<double> weight = SEMBLE::toScalar( (*sub)(1,1) );
          entry.push_back( std::make_pair( weight, helicity)); 
          tab[0] = entry;
        }
        else
        {
          for(int row = 0; row < cub->dim(); ++row)
          {
            SubduceTableMap::sub_list entry;
            std::complex<double> weightp = SEMBLE::toScalar( (*sub)(row+1,1));
            std::complex<double> weightm = SEMBLE::toScalar( (*sub)(row+1,2));

            entry.push_back(std::make_pair(weightp,helicity));
            entry.push_back(std::make_pair( eta_p * weightm , -helicity ) ); 

            tab[row] = entry; 
          }
        }

        delete sub;

        SubduceTableMap::irrep_sub_table* ret = new SubduceTableMap::irrep_sub_table(tab,cont,cub); 
        return ret; 
      }

    std::string print_table( const std::string &map_id)
    {
      const SubduceTableMap::irrep_sub_table* foo;
      foo = TheSmarterSubduceTableMap::Instance().get_table(map_id); 
      SubduceTableMap::sub_list::const_iterator list_it; 
      std::stringstream ss; 
      ss << map_id << std::endl;

      for(int i = 0; i < foo->size(); ++i)
      {
        ss << i << std::endl;
        for( list_it = foo->begin(i+1); list_it != foo->end(i+1); ++list_it) 
          ss << "[" << list_it->first << "]x<" << list_it->second << ">    ";  
        ss << std::endl;
      }

      return ss.str(); 
    }

    void add_subduce_table( const std::string &table_id, 
        const rHandle<LorentzRep> &cont, 
        const rHandle<CubicRep> &cub)
    {
      std::string map_id = make_subduce_table_map_id(cont,cub); 

      if( TheSmarterSubduceTableMap::Instance().mappy.find(map_id)
          == TheSmarterSubduceTableMap::Instance().mappy.end() )
      {
        printer_function<subduce_table_printer>( "attempting " + map_id ); 
        if( cub->rep_group() == Stringify<Oh>() )
          TheSmarterSubduceTableMap::Instance().mappy.insert(
              std::make_pair( map_id , subduce_rest(table_id, cont , cub ) ) ); 
        else
          TheSmarterSubduceTableMap::Instance().mappy.insert(
              std::make_pair( map_id , subduce_flight(table_id, cont , cub ) ) ); 

      //  printer_function<subduce_table_printer>( "made a " + map_id ); 
      //  printer_function<subduce_table_printer>( print_table( map_id) ); 
      }
    }


    std::vector<std::string> 
      all_subduce_table_keys()
      {
        return Hadron::TheSubduceTableFactory::Instance().keys(); 
      }

    std::vector<std::string> 
      spher_strings(const rHandle<LorentzRep> &s)
      {
        std::vector<std::string> ret; 
        std::stringstream ss; 
        ss << "J" << s->rep_spin(); 
        ret.push_back(ss.str());
        ss.str("");

        ss << "H0"; 
        s->rep_eta_p() == 1 ? ss << "+" : ss << "-";
        ret.push_back(ss.str());


        for(int i = 1; i <= s->rep_spin(); ++i)
        {
          ss.str("");
          ss << "H" << i; 
          ret.push_back(ss.str()); 
        }

        return ret; 
      };

    // string is the subduce table key, rep is what radmat calls the rep 
    //
    // do a double loop over possible left and right part of subduce table
    // keys to see if it exists 
    std::vector<std::pair<std::string,rHandle<CubicRep> > >
      possible_subductions( const rHandle<LorentzRep> &s )
      {
        std::vector<std::pair<std::string,rHandle<CubicRep> > > ret; 
        std::string delim = "->";
        std::string tail = ",1";
        std::vector<std::string> subduce_tables = all_subduce_table_keys();

        std::vector<std::string> sph_strs = spher_strings(s);
        std::vector<std::string> cub_strs = ::radmat::CubicRepresentationFactoryEnv::cubic_keys(); 
        std::vector<std::string>::const_iterator sph_it,cub_it; 

        for(sph_it = sph_strs.begin(); sph_it != sph_strs.end(); ++sph_it)
          for(cub_it = cub_strs.begin(); cub_it != cub_strs.end(); ++cub_it)
          {
            std::string id = *sph_it + delim + *cub_it + tail; 
            if( std::find( subduce_tables.begin(), subduce_tables.end(), id ) != subduce_tables.end() )
            {
              printer_function<subduce_printer>(id); 
              ret.push_back(std::make_pair(id,::radmat::CubicRepresentationFactoryEnv::callFactory(*cub_it)));
            }
          }

        return ret; 
      }


    bool work_handler()
    {
      std::vector<std::string> sph_strs = radmat::LorentzRepresentationFactoryEnv::spher_keys(); 
      std::vector<std::string>::const_iterator sph_it; 
      for(sph_it = sph_strs.begin(); sph_it != sph_strs.end(); ++sph_it)
      {
        rHandle<LorentzRep> rep = LorentzRepresentationFactoryEnv::callFactory(*sph_it); 
        std::vector<std::pair<std::string,rHandle<CubicRep> > > possible_sub_list; 
        possible_sub_list = possible_subductions(rep); 

        std::vector<std::pair<std::string,rHandle<CubicRep> > >::const_iterator it; 
        for(it = possible_sub_list.begin(); it != possible_sub_list.end(); ++it)
          add_subduce_table(it->first,rep,it->second); 

      } // sph_it

      return true; 
    }

  } // anonomyous 




  std::string make_subduce_table_map_id(const rHandle<LorentzRep> &cont, 
      const rHandle<CubicRep> &cub)
  {
    return cont->rep_id() + cub->rep_id(); 
  }

  namespace  TheSmarterSubduceTableMapFactoryEnv
  {
    namespace
    {
      bool local_registration = false; 
    }

    std::vector<std::string> all_keys()
    {
      std::vector<std::string> ret; 
      SubduceTableMap::map_t::const_iterator it; 
      typedef TheSmarterSubduceTableMap sstm; 
      ret.reserve(sstm::Instance().mappy.size()); 
      for(it = sstm::Instance().mappy.begin(); it != sstm::Instance().mappy.end(); ++it)
        ret.push_back(it->first); 

      return ret; 
    }

    bool registerAll()
    {
      bool success = true; 

      if( !!! local_registration )
      {
        success &= work_handler(); 
        local_registration = true; 
      }

      return success; 
    } 

  } // TheSmarterSubduceTableMap

}

