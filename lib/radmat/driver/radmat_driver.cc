/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : radmat_driver.cc

 * Purpose :

 * Creation Date : 25-02-2013

 * Last Modified : Mon 10 Nov 2014 10:24:37 AM EST

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/

#include "radmat/driver/radmat_driver.h"
#include "radmat/driver/radmat_driver_aux.h"
#include "radmat/utils/splash.h"
#include "radmat/driver/radmat_driver_props.h"
#include "radmat/data_representation/data_representation.h"
#include "radmat/rotation_interface/rotation_interface.h"
#include "jackFitter/plot.h"
#include "semble/semble_semble.h"
#include "io/adat_xmlio.h"
#include "ensem/ensem.h"
#include "adat/adat_stopwatch.h"
#include "hadron/ensem_filenames.h"
#include "formfac/formfac_qsq.h"

#include <complex>
#include <map>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <sstream>
#include <exception>
#include <algorithm>

using namespace radmat;
using namespace ENSEM;
using namespace ADAT;
using namespace ADATIO;



#define LOAD_LLSQ_PARALLEL 
#define FIT_LLSQ_PARALLEL
#define CHISQ_ANALYSIS_PARALLEL

#ifdef LOAD_LLSQ_PARALLEL 
#include <omp.h>
#endif 

#ifdef FIT_LLSQ_PARALLEL
#include <omp.h>
#endif 

#ifdef CHISQ_ANALYSIS_PARALLEL
#include <omp.h>
#endif 

namespace radmat
{

  namespace 
  {
    void 
      do_void_return_print(const std::string &s)
      {
        std::cout << s << std::endl; 
        return;
      }


    struct ff_data_store
    {
      ff_data_store() {}
      ff_data_store( const ENSEM::EnsemReal &QQ, 
          const ENSEM::EnsemReal &FFF, 
          const ThreePointDataTag &t)
        : Q2(QQ) , FF(FFF) , tag(t) 
      { }

      ENSEM::EnsemReal Q2; 
      ENSEM::EnsemReal FF;
      ThreePointDataTag tag; 
    };

    struct ff_data_store_sort_class
    {
      bool operator()(const ff_data_store &l, const ff_data_store &r)
      {
        return SEMBLE::toScalar(ENSEM::mean(l.Q2)) < SEMBLE::toScalar(ENSEM::mean(r.Q2)); 
      }
    } ff_data_sort;


    bool is_canon_mom_insertion(const Hadron::KeyHadronNPartNPtCorr_t &npt)
    {
      return npt.npoint[2].irrep.mom == FF::canonicalOrder(npt.npoint[2].irrep.mom); 
    }

  } // anonomyous 


  void RadmatDriver::run_program(const std::string &inifile)
  {
    init_false(); 
    if ( !!! read_xmlini(inifile))
      return do_void_return_print("failed to read xml ini");    

    if(m_ini.maxThread < 1)
    {
      std::stringstream err;
      err << "\n\n\n" << __PRETTY_FUNCTION__ << ": maxThread must be positive " 
        << "(maxthread = " << m_ini.maxThread << ")\n\n\n" << std::endl;
      check_exit(false,err.str().c_str());  
    }

    // If the scoping works the way omp says it does then 
    // this should control the max number of threads for the entire program 
    omp_set_num_threads(m_ini.maxThread);

    // build corrs or abort
    if ( !!! build_correlators() )
      return do_void_return_print("failed to build_corrs"); 

    // solve llsq or abort
    if( !!! solve_llsq() ) 
      return do_void_return_print("failed to solve llsq");

    // fit ffs or abort
    if( !!! fit_ffs() ) 
      return do_void_return_print("falied to fit ffs"); 

    // make plots or abort
    if( !!! make_FF_of_Q2_plots() ) 
      return do_void_return_print("failed to make FF of Q2 plots");

    // do chisq or abort
    if( !!! do_chisq_analysis() ) 
      return do_void_return_print("failed to do chisq analysis"); 
  }


  void RadmatDriver::xml_handler(const std::string &ini, const std::string &mode)
  {
    std::map<std::string, void (RadmatDriver::*)(const std::string &)> handler; 
    std::map<std::string, void (RadmatDriver::*)(const std::string &)>::iterator it; 
    handler["all"] = &RadmatDriver::build_xml; 
    handler["split_mom"] = &RadmatDriver::build_xml_split_p2;
    handler["canon_mom"] = &RadmatDriver::build_xml_canon_p2;
    handler["split_mom_N"] = &RadmatDriver::build_xml_split_p2_N;
    handler["two_point"] = &RadmatDriver::build_xml_twopoint;
    handler["split"] = &RadmatDriver::build_xml_split;


    it = handler.find(mode); 

    if (it != handler.end())
    {
      std::cout << __PRETTY_FUNCTION__ << ": FOUND " << mode << std::endl;
      (this->*(it->second))(ini);
    }
    else
    {
      std::cerr << __PRETTY_FUNCTION__ << ": error mode, " << mode
        << " not recognized, try one of the following" << std::endl;
      for (it = handler.begin(); it != handler.end(); ++it)
        std::cout << it->first << std::endl; 
    }
  }


  void RadmatDriver::build_xml(const std::string &inifile)
  {
    read_xmlini(inifile);

    std::vector<Hadron::KeyHadronNPartNPtCorr_t> keys;
    keys = m_correlators.construct_correlator_xml(m_ini.threePointIni); 

    ADATXML::XMLBufferWriter corrs;
    ADATXML::Array<Hadron::KeyHadronNPartNPtCorr_t> bc;

    bc.resize(keys.size()); 
    for(unsigned int i = 0; i < keys.size(); ++i)
      bc[i] = keys[i];

    write(corrs,"NPointList",bc);

    std::ofstream out("npt.list.xml");
    corrs.print(out);
    out.close();

    std::vector<Hadron::KeyHadronNPartNPtCorr_t>::const_iterator it; 

    out.open("npt.ensemFileNames.list"); 
    for(it = keys.begin(); it != keys.end(); ++it)
      out << Hadron::ensemFileName(*it) << "\n";
    out.close(); 
  }

  namespace
  {
    std::string mom_str(const ADATXML::Array<int> &p)
    {
      std::stringstream ss; 
      ss << "p" << p[0] << p[1] << p[2]; 
      return ss.str(); 
    }

    std::string can_mom_str(const Hadron::KeyHadronNPartNPtCorr_t &k)
    {
      return mom_str( FF::canonicalOrder(k.npoint[2].irrep.mom) );
    }




    template<int,int,int> 
      int split_mom_key(void)
      {
        __builtin_trap(); 
      }

    int modulus = 10; 

    // D3
    template<> int split_mom_key<1,1,1>()   { return  1 % modulus; }
    template<> int split_mom_key<-1,-1,-1>(){ return  2 % modulus; }
    template<> int split_mom_key<1,1,-1>()  { return  3 % modulus; }
    template<> int split_mom_key<1,-1,1>()  { return  4 % modulus; }
    template<> int split_mom_key<-1,1,1>()  { return  5 % modulus; }
    template<> int split_mom_key<1,-1,-1>() { return  6 % modulus; }
    template<> int split_mom_key<-1,-1,1>() { return  7 % modulus; }
    template<> int split_mom_key<-1,1,-1>() { return  8 % modulus; }

    // D2
    template<> int split_mom_key<1,1,0>()   { return  9 % modulus; }
    template<> int split_mom_key<1,0,1>()   { return 10 % modulus; }
    template<> int split_mom_key<0,1,1>()   { return 11 % modulus; }
    template<> int split_mom_key<1,-1,0>()  { return 12 % modulus; }
    template<> int split_mom_key<1,0,-1>()  { return 13 % modulus; }
    template<> int split_mom_key<0,1,-1>()  { return 14 % modulus; }
    template<> int split_mom_key<-1,1,0>()  { return 15 % modulus; }
    template<> int split_mom_key<-1,0,1>()  { return 16 % modulus; }
    template<> int split_mom_key<0,-1,1>()  { return 17 % modulus; }
    template<> int split_mom_key<-1,-1,0>() { return 18 % modulus; }
    template<> int split_mom_key<-1,0,-1>() { return 19 % modulus; }
    template<> int split_mom_key<0,-1,-1>() { return 20 % modulus; }

    // D4
    template<> int split_mom_key<0,0,1>()   { return 21 % modulus; }
    template<> int split_mom_key<0,0,-1>()  { return 22 % modulus; }
    template<> int split_mom_key<0,1,0>()   { return 23 % modulus; }
    template<> int split_mom_key<0,-1,0>()  { return 24 % modulus; }
    template<> int split_mom_key<1,0,0>()   { return 25 % modulus; }
    template<> int split_mom_key<-1,0,0>()  { return 26 % modulus; }

    // Oh
    template<> int split_mom_key<0,0,0>()   { return 27 % modulus; }


    template<int A,int B> 
      int split_mom_key(const ADATXML::Array<int> &i)
      {
        if( i[2] > 0 ) 
          return split_mom_key<A,B,1>(); 
        else if( i[2] == 0 )
          return split_mom_key<A,B,0>();
        else
          return split_mom_key<A,B,-1>(); 
      }


    template<int A> 
      int split_mom_key(const ADATXML::Array<int> &i)
      {
        if( i[1] > 0 ) 
          return split_mom_key<A,1>(i); 
        else if( i[1] == 0 )
          return split_mom_key<A,0>(i);
        else
          return split_mom_key<A,-1>(i); 
      }

    int split_mom_key(const ADATXML::Array<int> &i)
    {
      if( i[0] > 0 ) 
        return split_mom_key<1>(i); 
      else if( i[0] == 0 )
        return split_mom_key<0>(i);
      else
        return split_mom_key<-1>(i); 
    }


    Hadron::KeyHadronNPartNPtCorr_t
      twoPointCorr(const Hadron::KeyHadronNPartNPtCorr_t::NPoint_t &npt1,
          const Hadron::KeyHadronNPartNPtCorr_t::NPoint_t &npt2,
          const std::string &ensemble)
      {
        Hadron::KeyHadronNPartNPtCorr_t dest;
        dest.ensemble = ensemble; 
        dest.npoint.resize(2); 
        dest.npoint[1].t_slice = -2; 
        dest.npoint[2].t_slice = 0; 

        dest.npoint[1].irrep = npt1.irrep; 
        dest.npoint[2].irrep = npt2.irrep; 

        dest.npoint[1].irrep.creation_op = false; 
        dest.npoint[2].irrep.creation_op = true; 

        return dest; 
      }


    // make a twopoint list
    std::vector<Hadron::KeyHadronNPartNPtCorr_t> 
      twoPointList(const std::vector<Hadron::KeyHadronNPartNPtCorr_t::NPoint_t> &npts,
          const std::string &ensemble)
      {
        std::vector<Hadron::KeyHadronNPartNPtCorr_t::NPoint_t>::const_iterator a,b; 
        std::vector<Hadron::KeyHadronNPartNPtCorr_t> dest; 

        //  Square correlation matrix 
        //        for (a = npts.begin(); a != npts.end(); ++a)
        //          for(b = npts.begin(); b != npts.end(); ++b)
        //            dest.push_back( twoPointCorr(*a,*b,ensemble) );
        //

        // Diagonal elements
        for(a = npts.begin(); a != npts.end(); ++a)
          dest.push_back( twoPointCorr(*a,*a,ensemble) );

        return dest; 
      }

  } // anonomyous

  void RadmatDriver::build_xml_split_p2(const std::string &inifile)
  {
    read_xmlini(inifile);

    std::vector<Hadron::KeyHadronNPartNPtCorr_t> keys;
    std::vector<Hadron::KeyHadronNPartNPtCorr_t>::const_iterator unsorted_it;
    std::map<std::string,std::vector<Hadron::KeyHadronNPartNPtCorr_t> > sorted; 
    std::map<std::string,std::vector<Hadron::KeyHadronNPartNPtCorr_t> >::iterator sorted_it; 

    keys = m_correlators.construct_correlator_xml(m_ini.threePointIni); 

    std::cout << __PRETTY_FUNCTION__ << " writing xml for " 
      << keys.size() << " correlators" << std::endl; 

    // sort them based on momentum
    for(unsorted_it = keys.begin(); unsorted_it != keys.end(); ++unsorted_it)
    {
      sorted_it = sorted.find(can_mom_str(*unsorted_it)); 
      if (sorted_it != sorted.end())
        sorted_it->second.push_back(*unsorted_it); 
      else
        sorted.insert(
            std::pair<std::string,std::vector<Hadron::KeyHadronNPartNPtCorr_t> >(can_mom_str(*unsorted_it),
              std::vector<Hadron::KeyHadronNPartNPtCorr_t>(1,*unsorted_it) 
              )
            );
    }


    // now run them with some unique ids 
    for(sorted_it = sorted.begin(); sorted_it != sorted.end(); ++sorted_it)
    {
      ADATXML::XMLBufferWriter corrs;
      ADATXML::Array<Hadron::KeyHadronNPartNPtCorr_t> bc;
      keys = sorted_it->second; 

      bc.resize(keys.size()); 
      for(unsigned int i = 0; i < keys.size(); ++i)
        bc[i] = keys[i];

      write(corrs,"NPointList",bc);

      std::stringstream ss; 
      ss << "npt.list." << sorted_it->first << ".xml"; 

      std::ofstream out(ss.str().c_str());
      corrs.print(out);
      out.close();
    }
  }


  void RadmatDriver::build_xml_canon_p2(const std::string &inifile)
  {
    read_xmlini(inifile);

    std::vector<Hadron::KeyHadronNPartNPtCorr_t> all_keys,keys;
    std::vector<Hadron::KeyHadronNPartNPtCorr_t>::const_iterator it;

    all_keys = m_correlators.construct_correlator_xml(m_ini.threePointIni); 
    for(it = all_keys.begin(); it != all_keys.end(); ++it)
      if(is_canon_mom_insertion(*it))
        keys.push_back(*it); 


    ADATXML::XMLBufferWriter corrs;
    ADATXML::Array<Hadron::KeyHadronNPartNPtCorr_t> bc;

    bc.resize(keys.size()); 
    for(unsigned int i = 0; i < keys.size(); ++i)
      bc[i] = keys[i];

    write(corrs,"NPointList",bc);

    std::ofstream out("npt.list.canon.xml");
    corrs.print(out);
    out.close();

    out.open("npt.ensemFileNames.list"); 
    for(it = keys.begin(); it != keys.end(); ++it)
      out << Hadron::ensemFileName(*it) << "\n";
    out.close(); 

  }




  void RadmatDriver::build_xml_split_p2_N(const std::string &inifile)
  {
    const int N = 7; 
    read_xmlini(inifile);

    std::vector<Hadron::KeyHadronNPartNPtCorr_t> keys;
    std::vector<Hadron::KeyHadronNPartNPtCorr_t>::const_iterator unsorted_it;
    std::map<std::string,std::vector<Hadron::KeyHadronNPartNPtCorr_t> > sorted; 
    std::map<std::string,std::vector<Hadron::KeyHadronNPartNPtCorr_t> >::iterator sorted_it; 

    keys = m_correlators.construct_correlator_xml(m_ini.threePointIni); 

    std::cout << __PRETTY_FUNCTION__ << " writing xml for " 
      << keys.size() << " correlators" << std::endl; 

    std::map<std::string,std::vector<std::string> > seen; 
    std::map<std::string,std::vector<std::string> >::const_iterator seen_it; 
    std::map<std::string,char> letter_map; 
    char letter_start = 'a';

    //  now figure out how to split them up 
    std::map<std::string,std::string> mom_split_map; 
    std::map<std::string,std::string>::const_iterator mom_split_map_it; 

    // sort them based on momentum
    for(unsorted_it = keys.begin(); unsorted_it != keys.end(); ++unsorted_it)
    {
      std::string mom = mom_str(unsorted_it->npoint[2].irrep.mom); 
      std::string can_mom = can_mom_str(*unsorted_it); 

      if(letter_map.find(can_mom) == letter_map.end())
        letter_map.insert(std::make_pair(can_mom,letter_start));  

      mom_split_map_it = mom_split_map.find(mom); 

      // do nothing 
      if(mom_split_map_it != mom_split_map.end())
      {
        can_mom = mom_split_map_it->second; 
      }
      else
      {
        // check first time canonical mom 
        seen_it = seen.find(can_mom); 
        if(seen_it == seen.end())
          seen.insert(std::make_pair(can_mom, std::vector<std::string>())); 

        // add this momentum 
        std::vector<std::string> *see = &( seen.find(can_mom)->second ); 
        see->push_back(mom); 

        // increment 
        if(see->size() % N == 0)
          ++letter_map[can_mom]; 

        // redefine can mom 
        can_mom.append(1u,letter_map[can_mom]); 

        // now add it back to the mom_split map 
        mom_split_map.insert(std::make_pair(mom,can_mom)); 
      }

      // can_mom is now p110a or p110b etc  
      sorted_it = sorted.find(can_mom); 
      if (sorted_it != sorted.end())
        sorted_it->second.push_back(*unsorted_it); 
      else
        sorted.insert(
            std::pair<std::string,std::vector<Hadron::KeyHadronNPartNPtCorr_t> >(can_mom,
              std::vector<Hadron::KeyHadronNPartNPtCorr_t>(1,*unsorted_it) 
              )
            );
    }


    // now run them with some unique ids 
    for(sorted_it = sorted.begin(); sorted_it != sorted.end(); ++sorted_it)
    {
      ADATXML::XMLBufferWriter corrs;
      ADATXML::Array<Hadron::KeyHadronNPartNPtCorr_t> bc;
      keys = sorted_it->second; 

      bc.resize(keys.size()); 
      for(unsigned int i = 0; i < keys.size(); ++i)
        bc[i] = keys[i];

      write(corrs,"NPointList",bc);

      std::stringstream ss; 
      ss << "npt.list." << sorted_it->first << ".xml"; 

      std::ofstream out(ss.str().c_str());
      corrs.print(out);
      out.close();
    }
  }


  void RadmatDriver::build_xml_split(const std::string &inifile)
  {
    read_xmlini(inifile);

    std::vector<Hadron::KeyHadronNPartNPtCorr_t> keys;
    std::vector<Hadron::KeyHadronNPartNPtCorr_t>::const_iterator unsorted_it;
    std::map<int,std::vector<Hadron::KeyHadronNPartNPtCorr_t> > sorted; 
    std::map<int,std::vector<Hadron::KeyHadronNPartNPtCorr_t> >::iterator sorted_it; 

    keys = m_correlators.construct_correlator_xml(m_ini.threePointIni); 

    std::cout << __PRETTY_FUNCTION__ << " writing xml for " 
      << keys.size() << " correlators" << std::endl; 

    // sort them based on momentum -- some integer that comes from above (in a non biblical sense)
    for(unsorted_it = keys.begin(); unsorted_it != keys.end(); ++unsorted_it)
    {
      sorted_it = sorted.find(split_mom_key(unsorted_it->npoint[2].irrep.mom)); 
      if (sorted_it != sorted.end())
        sorted_it->second.push_back(*unsorted_it); 
      else
        sorted.insert(
            std::pair<int,std::vector<Hadron::KeyHadronNPartNPtCorr_t> >
            (
             split_mom_key(unsorted_it->npoint[2].irrep.mom),
             std::vector<Hadron::KeyHadronNPartNPtCorr_t>(1,*unsorted_it) 
            )
            );
    }


    // now run them with some unique ids 
    for(sorted_it = sorted.begin(); sorted_it != sorted.end(); ++sorted_it)
    {
      ADATXML::XMLBufferWriter corrs;
      ADATXML::Array<Hadron::KeyHadronNPartNPtCorr_t> bc;
      keys = sorted_it->second; 

      bc.resize(keys.size()); 
      for(unsigned int i = 0; i < keys.size(); ++i)
        bc[i] = keys[i];

      write(corrs,"NPointList",bc);

      std::stringstream ss; 
      ss << "npt.list.q" << sorted_it->first << ".xml"; 

      std::ofstream out(ss.str().c_str());
      corrs.print(out);
      out.close();
    }
  }

  void RadmatDriver::build_xml_twopoint(const std::string &inifile)
  {
    read_xmlini(inifile);

    std::vector<Hadron::KeyHadronNPartNPtCorr_t> keys;
    std::vector<Hadron::KeyHadronNPartNPtCorr_t>::const_iterator kit;
    keys = m_correlators.construct_correlator_xml(m_ini.threePointIni); 

    if ( keys.size() <= 0 ) 
      exit(12034); 

    std::vector<Hadron::KeyHadronNPartNPtCorr_t::NPoint_t> npts;

    for (kit = keys.begin(); kit != keys.end(); ++kit)
      npts.push_back(kit->npoint[1]);


    std::vector<Hadron::KeyHadronNPartNPtCorr_t> list = twoPointList( npts, keys[0].ensemble );

    ADATXML::XMLBufferWriter corrs;
    ADATXML::Array<Hadron::KeyHadronNPartNPtCorr_t> bc;

    bc.resize(list.size()); 
    for(unsigned int i = 0; i < keys.size(); ++i)
      bc[i] = list[i];

    write(corrs,"NPointList",bc);

    std::ofstream out("npt.list.xml");
    corrs.print(out);
    out.close();

  }

  //  void RadmatDriver::nuke_graph(const std::string &inifile, 
  //      const std::string &graph_db,
  //      const std::string &nuke_xml_out)
  //  {
  //
  //    read_xmlini(inifile);
  //
  //    std::vector<Hadron::KeyHadronNPartNPtCorr_t> keys;
  //    keys = m_correlators.build_correlator_xml(m_ini.threePointIni); 
  //
  //    DisconnectedGraphNuker n;
  //    n.find_nukes(keys,graph_db); 
  //    n.dump_nukes(nuke_xml_out); 
  //  }
  //
  //
  //  void RadmatDriver::build_stub_xml(const std::string &inifile)
  //  {
  //    read_xmlini(inifile);
  //
  //    std::vector<Hadron::KeyHadronNPartNPtCorr_t> keys;
  //    keys = m_correlators.build_correlator_xml(m_ini.threePointIni); 
  //
  //
  //    stubify(keys);
  //
  //    ADATXML::XMLBufferWriter corrs;
  //    ADATXML::Array<Hadron::KeyHadronNPartNPtCorr_t> bc;
  //
  //    bc.resize(keys.size()); 
  //    for(unsigned int i = 0; i < keys.size(); ++i)
  //      bc[i] = keys[i];
  //
  //    write(corrs,"NPointList",bc);
  //
  //    std::ofstream out("npt.list.xml");
  //    corrs.print(out);
  //    out.close();
  //
  //    std::vector<Hadron::KeyHadronNPartNPtCorr_t>::const_iterator it; 
  //
  //    out.open("npt.ensemFileNames.list"); 
  //    for(it = keys.begin(); it != keys.end(); ++it)
  //      out << Hadron::ensemFileName(*it) << "\n";
  //    out.close(); 
  //
  //  }



  void RadmatDriver::init_false(void)
  {
    read_ini = false;
    built_correlators = false;
    init_llsq = false; 
    solved_llsq = false;
    fit_formfacs = false;
    chisq_analysis = false; 
  }

  void RadmatDriver::check_exit(const bool &b, const char *c) const 
  {
    if(!!!b)
    {
      std::cerr << __func__ << ": error: called by " << c << ", exiting." << std::endl;
      exit(1); 
    }
  }


  bool RadmatDriver::read_xmlini(const std::string &xmlini)
  {


    try
    {
      XMLReader xml(xmlini);
      read(xml,"/DriverProps",m_ini);
    }
    catch(std::exception &e)
    {
      std::cout << "std exception: " << e.what();
    }
    catch(std::string &e)
    {
      std::cout << __func__ << ": ERROR: can't read xmlinifile ("
        << xmlini << "): " << e << std::endl;
      exit(1);
    }
    catch(...)
    {
      SPLASH("An error occured while loading the inifile");
      exit(1);
    }

    read_ini = true; 

    return read_ini;
  }



  bool RadmatDriver::build_correlators(void)
  {
    check_exit_ini(); 

    Util::StopWatch my_stopwatch; 
    my_stopwatch.start(); 

    std::cout << "Loading correlators.. " << std::endl; 

    multi_lattice_data = m_correlators.construct_multi_correlators(m_ini.threePointIni);

    my_stopwatch.stop(); 
    std::cout << "Loading correlators and constructing matrix elements took " 
      << my_stopwatch.getTimeInSeconds() << " seconds " << std::endl;

    built_correlators = true; 

    return built_correlators; 
  }

  bool RadmatDriver::solve_llsq(void)
  {
    check_exit_corrs(); 

    int idx, sz = multi_lattice_data.size(); 
    std::string soln_ID = std::string ("SVDNonSquare");

    if(sz == 0)
    {
      std::cerr << __func__ << ": error nothing to solve" << std::endl;
      return false;  
    }

    linear_systems_of_Q2.resize(multi_lattice_data.size());

    std::cout << "Solving LLSQ.. " << std::endl;

    Util::StopWatch my_stopwatch; 
    my_stopwatch.start(); 


    good_qs.resize(sz,false); 

    std::cout << __func__ << ": sz = " << sz << std::endl;

    std::map<std::string,bool> irrep_ops; 
    irrep_ops.insert(std::make_pair("subduce",false)); 
    irrep_ops.insert(std::make_pair("mix_irreps",true)); 
    bool mix_irreps = irrep_ops.find( m_ini.threePointIni.matElemMode )->second; 

#ifdef LOAD_LLSQ_PARALLEL 

#pragma omp parallel for shared(idx)  schedule(dynamic,1)

#endif 
    // POSSIBLE PARALLEL HERE
    for(idx =0; idx < sz; ++idx)
    {
      good_qs[idx] =  linear_systems_of_Q2[idx].load_llsq(multi_lattice_data[idx],
          m_ini.poleMass,
          m_ini.tolerance,
          mix_irreps);
    }
    // END PARALLEL

#pragma omp barrier

    int ngood(0);
    for(idx = 0; idx < sz; ++idx)
      if(good_qs[idx])
        ++ngood; 

    init_llsq = true;

    std::cout << __func__ << ": " << ngood << " good Q^2 points out of " << sz << std::endl;



    // print the list here in case the solver flakes we can easily determine where it went wrong
    print_Q2_list(); 


    // threading over q2
#ifdef LOAD_LLSQ_PARALLEL 
#pragma omp parallel for shared(idx)  schedule(dynamic,1)
#endif   
    for(idx = 0; idx < sz; ++idx)
      if(good_qs[idx])
        linear_systems_of_Q2[idx].solve_llsq(soln_ID); 
#pragma omp barrier

    my_stopwatch.stop();
    std::cout << "Solving LLSQ took "     
      << my_stopwatch.getTimeInSeconds() << " seconds " << std::endl;

    solved_llsq = true;


    // save the llsq state into a database
#ifdef LOAD_LLSQ_PARALLEL // thread this  
#pragma omp parallel for shared(idx)  schedule(dynamic,1)
#endif 
    for(idx = 0; idx < sz; ++idx)
      if(good_qs[idx])
        linear_systems_of_Q2[idx].save_llsq_state(); 

    // leave a barrier since to prevent any possibility of a jump out from below
#pragma omp barrier

    // save the ff_of_t
#ifdef LOAD_LLSQ_PARALLEL // thread this  
#pragma omp parallel for shared(idx)  schedule(dynamic,1)
#endif 
    for(idx = 0; idx < sz; ++idx)
      if(good_qs[idx])
        linear_systems_of_Q2[idx].save_ff_of_t(); 

    // leave a barrier since to prevent any possibility of a jump out from below
#pragma omp barrier

    // save the result 
#ifdef LOAD_LLSQ_PARALLEL // thread this  
#pragma omp parallel for shared(idx)  schedule(dynamic,1)
#endif 
    for(idx = 0; idx < sz; ++idx)
      if(good_qs[idx])                        // can only print the successful guys
        linear_systems_of_Q2[idx].dump_llsq(); 
#pragma omp barrier

    return solved_llsq;
  }

  bool RadmatDriver::fit_ffs(void)
  {
    check_exit_llsq(); 

    int idx, sz = multi_lattice_data.size(); 

    std::cout << "Fitting FF(t_ins) " << std::endl;

    Util::StopWatch my_stopwatch; 
    my_stopwatch.start(); 

    ADATXML::Array<int> timeslice_info;
    timeslice_info = m_ini.threePointIni.threePointCorrXMLIni.redstar->timeslice_info(); 
    POW2_ASSERT(timeslice_info.size() == 3); 
    int tsrc,tsnk;

    tsnk = timeslice_info[0]; 
    tsrc = timeslice_info[2]; 

    POW2_ASSERT(tsrc < tsnk); 

#ifdef FIT_LLSQ_PARALLEL
#pragma omp parallel for shared(idx,tsrc,tsnk)  schedule(dynamic,1)
#endif 
    // POSSIBLE PARALLEL HERE
    for(idx = 0; idx < sz; ++idx)
      if(good_qs[idx])
        linear_systems_of_Q2[idx].fit_data(m_ini.threePointComparatorProps,tsrc,tsnk);

#pragma omp barrier

    my_stopwatch.stop();
    std::cout << "Fitting took " 
      << my_stopwatch.getTimeInSeconds() << " seconds " << std::endl;

    fit_formfacs = true; 

#ifdef FIT_LLSQ_PARALLEL
#pragma omp parallel for shared(idx,tsrc,tsnk)  schedule(dynamic,1)
#endif 
    // POSSIBLE PARALLEL HERE
    for(idx = 0; idx < sz; ++idx)
      if(good_qs[idx])                        // can only print the successful guys
        linear_systems_of_Q2[idx].dump_fits(); 

#pragma omp barrier

    return fit_formfacs;
  }


  bool RadmatDriver::do_chisq_analysis(void)
  {
    if ( m_ini.chisq == "none") 
      return true; // hack job here

    check_exit_fit();

    int idx, sz = multi_lattice_data.size(); 

    std::cout << "chisq_analysis" << std::endl;

    int low, high ; 
    low = m_ini.threePointComparatorProps.tlow; 
    high = m_ini.threePointComparatorProps.thigh;


    Util::StopWatch my_stopwatch; 
    my_stopwatch.start(); 

#ifdef CHISQ_ANALYSIS_PARALLEL

#pragma omp parallel for shared(sz) schedule(dynamic,1)

#endif

    // POSSIBLE PARALLEL HERE
    for(idx = 0; idx < sz; ++idx)
      if(good_qs[idx])
        linear_systems_of_Q2[idx].chisq_analysis(low,high);

#pragma omp barrier

    my_stopwatch.stop();
    std::cout << "chisq_analysis took " 
      << my_stopwatch.getTimeInSeconds() << " seconds " << std::endl;

    chisq_analysis = true; 

    return chisq_analysis;
  }


  bool RadmatDriver::make_FF_of_Q2_plots(void)
  {
    std::string pth = SEMBLE::SEMBLEIO::getPath();
    std::stringstream path;
    path << pth << "FF_of_Q2/";
    SEMBLE::SEMBLEIO::makeDirectoryPath(path.str()); 

    // single Q2 data solutions
    std::vector<RadmatSingleQ2Solution> thingy; 
    std::vector<RadmatSingleQ2Solution>::const_iterator it; 
    thingy.reserve( good_qs.size() ); 

    // only pull down the good ones 
    for(int Q = 0; Q < linear_systems_of_Q2.size(); ++Q)
      if(good_qs[Q])
        thingy.push_back( linear_systems_of_Q2[Q].fetchSolution() ); 

    //key is name of form factor  
    typedef std::map<std::string, std::vector<ff_data_store> > ff_map_t; 
    ff_map_t ff_map;
    ff_map_t::iterator ff_map_it; 
    std::map<std::string,ENSEM::EnsemReal>::const_iterator little_it; 

    // reorganize so the form factor name is the key, data is pair of ensems 
    for(it = thingy.begin(); it != thingy.end(); ++it)
      for(little_it = it->ff_map.begin(); little_it != it->ff_map.end(); ++little_it)
      {
        std::string id = little_it->first; 
        ff_map_it = ff_map.find( id ); 

        // if it is not in the ff_map then insert 
        if( ff_map_it == ff_map.end() )
        {
          std::vector<ff_data_store> data; 
          data.push_back(ff_data_store(it->Q2,little_it->second,*(it->tags.begin()))); 
          ff_map.insert(std::make_pair(id,data)); 
        }
        else // grow the data 
        {
          ff_map_it->second.push_back( 
              ff_data_store(it->Q2,little_it->second,*(it->tags.begin()))); 
        }
      }

    // now loop and toss them into some plots/data files 
    ff_map_t::const_iterator ffit; 

    for( ffit = ff_map.begin(); ffit != ff_map.end(); ++ffit)
    {
      std::string id = ffit->first; 
      std::cout << __func__ << ": working on " << id << std::endl;
      std::vector<ff_data_store> data = ffit->second; 

      // run a sort on the data ( on the mean of Q2 ) 
      std::sort( data.begin() , data.end() , ff_data_sort ); 

      ENSEM::EnsemVectorReal FF; 
      FF.resize( data.begin()->Q2.size() ); 
      FF.resizeObs( data.size() );  

      std::vector<double> q(data.size()),qerr(data.size());  
      std::stringstream data_buffer; 

      data_buffer << "# Q2  Q2err FF FFerr | <left_rep> <left_mom> <right_rep> <right_mom> " << std::endl; 
      int count = 0; 


      // some cute little vector thingy 
      for(int i = 0; i < data.size(); ++i)
      {
        q[i] = SEMBLE::toScalar(ENSEM::mean( data[i].Q2 )); 
        qerr[i] = sqrt(SEMBLE::toScalar(ENSEM::variance( data[i].Q2 ))); 
        ENSEM::pokeObs(FF,data[i].FF,i); 

        ThreePointDataTag t = data[i].tag; 
        DataRep3pt data_rep = t.data_rep; 

        double ff = SEMBLE::toScalar(ENSEM::mean( data[i].FF )); 
        double fferr = sqrt(SEMBLE::toScalar(ENSEM::variance( data[i].FF ))); 

        // use the canonical momentum when reporting the solution 
        // if people pay attention then they will know that the llsq 
        // has been averaged over equivalent frames 
        std::pair<mom_t,mom_t> canonical_momentum; 
        canonical_momentum = radmat::LatticeRotationEnv::rotation_group_can_mom(t.left_mom,t.right_mom);  
        mom_t lp = canonical_momentum.first; 
        mom_t rp = canonical_momentum.second; 

        // bung in an overall phase 
        int bigPhase = m_ini.bigPhase; 
        double dphase = double(bigPhase); 

        // have all the data, put it into a useful log form 
        data_buffer << ++count; 
        data_buffer << " q2 " << q[i] << " " << qerr[i]; 
        data_buffer << " ff " << dphase*ff << " " << fferr; 
        data_buffer << " | " << t.full_irrep_id(data_rep,data_rep.l) 
          << " " << lp[0] << " " << lp[1] << " " << lp[2]; 
        data_buffer << " " << t.full_irrep_id(data_rep,data_rep.r) 
          << " " << rp[0] << " " << rp[1] << " " << rp[2]; 
        data_buffer << " g " << t.full_irrep_id(data_rep,data_rep.g); 
        data_buffer << " bigPhase " << dphase; 

        data_buffer << std::endl; 
      }

      // make a simple plot 
      AxisPlot plt; 
      plt.addEnsemData(q,FF,"\\sq",1);

      std::stringstream ss; 
      ss << path.str() << id << ".ax"; 
      plt.sendToFile(ss.str()); 

      // throw some junk to the screen so people think they 
      // actually did something 
      std::string cute_data = data_buffer.str(); 
      std::cout << cute_data << std::endl;

      // chunk out the useful info 
      std::stringstream s2; 
      s2 << path.str() << id << ".dat";
      std::ofstream out2( s2.str().c_str() ) ; 
      out2 << cute_data << std::endl;
      out2.close(); 

    } // end ffit 

    // success 
    return true; 
  }


  bool RadmatDriver::print_Q2_list(void) 
  {
    check_exit_init_llsq(); 
    std::stringstream ss;
    std::vector<RadmatSingleQ2Driver>::const_iterator it; 
    std::string pth = SEMBLE::SEMBLEIO::getPath() + std::string("Q2_to_mat_elems.txt"); 
    std::ofstream out(pth.c_str());
    std::string delim("--------------------------------\n"); 
    for(it = linear_systems_of_Q2.begin(); it != linear_systems_of_Q2.end(); ++it)
      if(it->check_linear_system())
        out << it->tags_at_this_Q2() << delim; 
    out.close();

    return true; 
  }











}


#undef LOAD_LLSQ_PARALLEL 
#undef FIT_LLSQ_PARALLEL
#undef CHISQ_ANALYSIS_PARALLEL



