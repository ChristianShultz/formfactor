/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : construct_correlators.cc

 * Purpose :

 * Creation Date : 04-12-2012

 * Last Modified : Mon 16 Jun 2014 09:42:15 AM EDT

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/



#include "construct_correlators.h"
#include "construct_correlators_utils.h"
#include "construct_correlators_bad_data_repository.h"
#include "radmat/utils/printer.h"
#include "adat/adat_stopwatch.h"
#include "adat/map_obj.h"
#include "hadron/ensem_filenames.h"
#include <string>
#include <vector>
#include <utility>

#define CONSTRUCT_CORRELATORS_PARALLEL
#ifdef CONSTRUCT_CORRELATORS_PARALLEL
#include <omp.h>
#endif 

#define TIME_CONSTRUCT_SINGLE_CORRS
#define TIME_CONSTRUCT_ALL_CORRS

// actually uses a subset since fitting fails with zero variance 
// #define USE_MEAN_CORRS_DEBUG_MODE 

namespace radmat
{

  namespace 
  {

    // garbage intermediary 
    template<typename T, typename U, typename V>
      struct triplet
      {
        triplet(void) {}
        triplet(const T &tt, const U &uu, const V &vv)
          : first(tt) , second(uu) , third(vv)
        {  }

        triplet(const std::pair<T,U> FandS, const V &vv)
          : first(FandS.first) , second(FandS.second) , third(vv)
        { }


        T first; 
        U second; 
        V third; 
      };


    ///////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////
    // decide if there is any data and move from tagged 
    // ensem redstar blocks to handles of LLSQL...Data type
    template<typename T> 
    std::pair<bool,rHandle<LLSQLatticeMultiData> >
      construct_lattice_data(
          const T &qsq,
          const std::vector<TaggedEnsemRedstarNPtBlock> &corrs,
          const std::string &sink_id, 
          const std::string &source_id, 
          const ThreePointCorrXMLIni_t::RenormalizationProp &Z_V,
          const DatabaseInterface_t &db )
      {

#ifdef TIME_CONSTRUCT_SINGLE_CORRS
        Util::StopWatch snoop;
        snoop.start();
#endif

        std::cout << __func__ << ": working on Q2 =" << qsq << std::endl;

        std::pair<bool , rHandle<LLSQLatticeMultiData> >  data;

#ifdef USE_MEAN_CORRS_DEBUG_MODE 
        std::cout << __func__ << ": WARNING COMPILED WITH MEAN CORRS, GARBAGE STATISTICS " << std::endl;
        data = build_correlators_no_copy_Ncfg_mean_fast(corrs,sink_id,source_id,Z_V,db,50); 
#else
        data = build_correlators_no_copy(corrs,sink_id,source_id,Z_V,db); 
#endif 

#ifdef TIME_CONSTRUCT_SINGLE_CORRS
        snoop.stop();
        if(data.first)
          std::cout << __func__ << ": q2 = " << qsq 
            << " NxM -> " << data.second->nrows() << "x" << data.second->ncols()
            << " took " << snoop.getTimeInSeconds() << " seconds" << std::endl;
        else
          std::cout << __func__ << ": q2 = " << qsq << " failure took " 
            << snoop.getTimeInSeconds() << " seconds " << std::endl;
#endif

        return data; 
      }



    ///////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////
    // wrap the above code and possible parallelization  
    // on the q2 sort values here 

    template<typename T> 
    std::vector< 
      triplet<bool, rHandle<LLSQLatticeMultiData>, T > 
      >
      build_llsq_corrs(
          const std::vector<std::pair<T,std::vector<TaggedEnsemRedstarNPtBlock> > > &loop_data,
          const std::string &snk_id, 
          const std::string &src_id,
          const ThreePointCorrXMLIni_t::RenormalizationProp &Z_V,
          const DatabaseInterface_t &db)
      {

        std::vector<triplet<bool,rHandle<LLSQLatticeMultiData>, T > > ret(loop_data.size()); 
        int idx;
        int sz ( loop_data.size() ) ; 

#ifdef CONSTRUCT_CORRELATORS_PARALLEL
#pragma omp parallel for shared(idx) schedule(dynamic,1)
#endif

        for(idx = 0; idx < sz; ++idx)
          ret[idx] =
            triplet<bool, rHandle<LLSQLatticeMultiData>, T > 
            (construct_lattice_data( loop_data[idx].first, loop_data[idx].second, 
                                     snk_id, src_id, Z_V, db), loop_data[idx].first);

#ifdef CONSTRUCT_CORRELATORS_PARALLEL
#pragma omp barrier
#endif 

        return ret; 
      }

    // template xml handler 
    template<typename T> 
      std::vector<TaggedEnsemRedstarNPtBlock> 
      pull_data_xml_function( const ThreePointCorrIni_t &, 
          rHandle<AbsRedstarXMLInterface_t> &)
      { __builtin_trap(); return std::vector<TaggedEnsemRedstarNPtBlock>(); }

    // handle a lorentz type 
    template<>
      std::vector<TaggedEnsemRedstarNPtBlock> 
      pull_data_xml_function<RedstarThreePointXMLLorentzHandler>( const ThreePointCorrIni_t &ini, 
          rHandle<AbsRedstarXMLInterface_t> &handle)
      {
        RedstarThreePointXMLLorentzHandler *red; 
        red = dynamic_cast<RedstarThreePointXMLLorentzHandler*>(handle.get_ptr()); 

        double p_factor = mom_factor(ini.xi,ini.L_s); 
        std::string elem_id = ini.matElemID; 
        const radmatDBProp_t *db_prop = &ini.radmatDBProp; 
        const ThreePointCorrXMLIni_t *three_pt = &ini.threePointCorrXMLIni; 

        RedstarThreePointXMLInput input; 
        input.db_props = ini.radmatDBProp; 
        input.pid_left = ini.threePointCorrXMLIni.sink_id; 
        input.pid_right = ini.threePointCorrXMLIni.source_id; 
        input.mom_fac = p_factor; 


        return tag_lattice_xml(
            red->handle_work(input),
            p_factor, 
            three_pt->maSink, 
            three_pt->maSource, 
            elem_id); 
      }


    // handle a lorentz type 
    template<>
      std::vector<TaggedEnsemRedstarNPtBlock> 
      pull_data_xml_function<RedstarThreePointXMLSubduceHandler>( const ThreePointCorrIni_t &ini, 
          rHandle<AbsRedstarXMLInterface_t> &handle)
      {
        RedstarThreePointXMLSubduceHandler *red; 
        red = dynamic_cast<RedstarThreePointXMLSubduceHandler*>(handle.get_ptr()); 

        double p_factor = mom_factor(ini.xi,ini.L_s); 
        std::string elem_id = ini.matElemID; 
        const radmatDBProp_t *db_prop = &ini.radmatDBProp; 
        const ThreePointCorrXMLIni_t *three_pt = &ini.threePointCorrXMLIni; 


        RedstarThreePointXMLInput input; 
        input.db_props = ini.radmatDBProp; 
        input.pid_left = ini.threePointCorrXMLIni.sink_id; 
        input.pid_right = ini.threePointCorrXMLIni.source_id; 
        input.mom_fac = p_factor; 

        // this must pick out the correct matrix element
        std::vector<TaggedEnsemRedstarNPtBlock> lorentz_tagged_data;
        lorentz_tagged_data =  tag_lattice_xml(
            red->handle_work(input),
            p_factor, 
            three_pt->maSink, 
            three_pt->maSource, 
            elem_id); 

        std::vector<TaggedEnsemRedstarNPtBlock>::iterator it; 
        for(it = lorentz_tagged_data.begin(); it != lorentz_tagged_data.end(); ++it)
          {
            // use a pointer to avoid a full copy 
            ThreePointDataTag * tag = &(it->data_tag); 
            rHandle<Rep_p> l_prim, r_prim; 
            l_prim = tag->data_rep.lefty(); 
            r_prim = tag->data_rep.righty(); 
            tag->mat_elem_id += "__" + l_prim->rep_id() + "," + r_prim->rep_id(); 
          }

        // actually has an updated cubic mat elem here despite the 
        // stupidity of the name 
        return lorentz_tagged_data; 
      }


    // handle a lorentz type 
    template<>
      std::vector<TaggedEnsemRedstarNPtBlock> 
      pull_data_xml_function<RedstarThreePointXMLMixedHandler>( const ThreePointCorrIni_t &ini, 
          rHandle<AbsRedstarXMLInterface_t> &handle)
      {
        RedstarThreePointXMLMixedHandler *red; 
        red = dynamic_cast<RedstarThreePointXMLMixedHandler*>(handle.get_ptr()); 

        double p_factor = mom_factor(ini.xi,ini.L_s); 
        std::string elem_id = ini.matElemID; 
        const radmatDBProp_t *db_prop = &ini.radmatDBProp; 
        const ThreePointCorrXMLIni_t *three_pt = &ini.threePointCorrXMLIni; 

        RedstarThreePointXMLInput input; 
        input.db_props = ini.radmatDBProp; 
        input.pid_left = ini.threePointCorrXMLIni.sink_id; 
        input.pid_right = ini.threePointCorrXMLIni.source_id; 
        input.mom_fac = p_factor; 

        // this must pick out the correct matrix element
        std::vector<TaggedEnsemRedstarNPtBlock> lorentz_tagged_data;
        lorentz_tagged_data =  tag_lattice_xml(
            red->handle_work(input),
            p_factor, 
            three_pt->maSink, 
            three_pt->maSource, 
            elem_id); 

        std::vector<TaggedEnsemRedstarNPtBlock>::iterator it; 
        for(it = lorentz_tagged_data.begin(); it != lorentz_tagged_data.end(); ++it)
          {
            // use a pointer to avoid a full copy 
            ThreePointDataTag * tag = &(it->data_tag); 
            rHandle<Rep_p> l_prim, r_prim; 
            l_prim = tag->data_rep.lefty(); 
            r_prim = tag->data_rep.righty(); 
            tag->mat_elem_id += "__" + l_prim->rep_id() + "," + r_prim->rep_id(); 
          }

        // actually has an updated cubic mat elem here despite the 
        // stupidity of the name 
        return lorentz_tagged_data; 
      }


    //////////////////////////////////////////////////////
    // pull data out of the xml class
    //      -- the name of the class defines the type of 
    //      representation we want out matrix element to 
    //      have but does not not not define the sorting 
    //      method that we apply to the data 
    std::vector<TaggedEnsemRedstarNPtBlock> 
      pull_data_xml(const ThreePointCorrIni_t &ini)
      {
        rHandle<AbsRedstarXMLInterface_t> red = ini.threePointCorrXMLIni.redstar; 

        if( red->type() == Stringify<RedstarThreePointXMLLorentzHandler>() )
          return pull_data_xml_function<RedstarThreePointXMLLorentzHandler>( ini, red ); 
        else if( red->type() == Stringify<RedstarThreePointXMLSubduceHandler>() )
          return pull_data_xml_function<RedstarThreePointXMLSubduceHandler>( ini, red ); 
        else if( red->type() == Stringify<RedstarThreePointXMLMixedHandler>() )
          return pull_data_xml_function<RedstarThreePointXMLMixedHandler>( ini, red ); 
        else
          printer_function<console_print>(" unknown type " + red->type()
              + " in pull_data_xml in construct_correlators.cc " ); 

        exit(1); 
      }



    ///////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////
    // the  brain 
    std::vector<rHandle<LLSQLatticeMultiData> >
      do_work_rotation_groups(const ThreePointCorrIni_t &ini)
      {

#ifdef TIME_CONSTRUCT_ALL_CORRS
        Util::StopWatch snoop;
        snoop.start(); 
#endif

        double p_factor = mom_factor(ini.xi,ini.L_s); 
        std::string elem_id = ini.matElemID; 
        const radmatDBProp_t *db_prop = &ini.radmatDBProp; 
        const ThreePointCorrXMLIni_t *three_pt = &ini.threePointCorrXMLIni; 

        std::vector<TaggedEnsemRedstarNPtBlock> unsorted_elems; 

        unsorted_elems = pull_data_xml(ini); 

        bool sort_mode; 
        std::map<std::string,bool> sort_mode_map; 
        sort_mode_map.insert( std::make_pair( "subduce" , false ) ); 
        sort_mode_map.insert( std::make_pair( "mix_irreps" , true ) ); 

        std::map<std::string,bool>::const_iterator sort_mode_it; 
        sort_mode_it = sort_mode_map.find( ini.matElemMode ); 
        if( sort_mode_it == sort_mode_map.end() ) 
        {
          std::cout << __PRETTY_FUNCTION__ << ": " << ini.matElemMode 
            << " is not supported, try one of the following" << std::endl;
          for(sort_mode_it = sort_mode_map.begin(); sort_mode_it != sort_mode_map.end(); ++sort_mode_it)
            std::cout << sort_mode_it->first << std::endl;
          exit(1); 
        }

        sort_mode = sort_mode_it->second; 


        std::map<std::string,std::vector<TaggedEnsemRedstarNPtBlock> > sorted_elems;
        std::map<std::string,std::vector<TaggedEnsemRedstarNPtBlock> >::const_iterator it;
        sorted_elems = sort_tagged_corrs_by_Q2_and_rotation_group(unsorted_elems, sort_mode); 

        std::vector<std::pair<std::string,std::vector<TaggedEnsemRedstarNPtBlock> > > loop_data; 

        for(it = sorted_elems.begin(); it != sorted_elems.end(); ++it)
          loop_data.push_back(std::pair<std::string,std::vector<TaggedEnsemRedstarNPtBlock> >(it->first,it->second)); 


        std::vector<triplet<bool, rHandle<LLSQLatticeMultiData>, std::string > > lattice_data; 
        DatabaseInterface_t db(*db_prop) ; 

        lattice_data = build_llsq_corrs(loop_data,three_pt->sink_id, three_pt->source_id, 
            three_pt->renormalization, db); 

        // stop and dump anything that we may be missing 
        ::radmat::BAD_DATA_REPO::dump_bad_data();

        std::vector<rHandle<LLSQLatticeMultiData> > ret; 
        std::vector<triplet<bool,rHandle<LLSQLatticeMultiData>,std::string > >::const_iterator dcheck; 

        for(dcheck = lattice_data.begin(); dcheck != lattice_data.end(); ++dcheck) 
        {
          // this is a have some data flag
          if( !!! dcheck->first )
          {
            std::cout << __func__ << ": dropping Q2 = " << dcheck->third 
              << " because it has 0 elems" << std::endl;  
            continue;    
          }
          ret.push_back(dcheck->second); 
        }

#ifdef TIME_CONSTRUCT_ALL_CORRS
        snoop.stop(); 
        std::cout << " ** time to build all correlators " << snoop.getTimeInSeconds() << " seconds" << std::endl;
        std::cout << " ** returning " << ret.size() << " possible Q2 points " << std::endl; 
#endif
        return ret;  
      }


  } // anonomyous 


  ///////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////


  ///////////////////////////////////////////////////////
  // construct lots of correlators
  std::vector<rHandle<LLSQLatticeMultiData> >
    ConstructCorrelators::construct_multi_correlators(void) const
    {
      return do_work_rotation_groups(m_ini); 
    }


  ///////////////////////////////////////////////////////
  // just make some xml 
  std::vector<Hadron::KeyHadronNPartNPtCorr_t>
    ConstructCorrelators::construct_correlator_xml(void) const
    {

#ifdef TIME_CONSTRUCT_ALL_CORRS
      Util::StopWatch snoop;
      snoop.start(); 
#endif

      double p_factor = mom_factor(m_ini.xi,m_ini.L_s); 
      std::string elem_id = m_ini.matElemID; 
      const ThreePointCorrXMLIni_t *three_pt = &m_ini.threePointCorrXMLIni; 

      std::vector<TaggedEnsemRedstarNPtBlock> unsorted_elems; 
      unsorted_elems = pull_data_xml(m_ini); 

      std::cout << __PRETTY_FUNCTION__ << ": $#unsorted = " << unsorted_elems.size() << std::endl; 

      std::vector<TaggedEnsemRedstarNPtBlock>::const_iterator block;
      EnsemRedstarNPtBlock::const_iterator npt; 
      std::map<std::string,EnsemRedstarNPtBlock::Obj_t> pull; 
      for(block = unsorted_elems.begin(); block != unsorted_elems.end(); ++block)
        for(npt = block->coeff_lattice_xml.begin(); npt != block->coeff_lattice_xml.end(); ++npt)
          pull[Hadron::ensemFileName(npt->m_obj)] = npt->m_obj; 


    std::vector<Hadron::KeyHadronNPartNPtCorr_t> keys; 
    keys.reserve(pull.size()); 
    std::map<std::string,EnsemRedstarNPtBlock::Obj_t>::const_iterator it; 
    for(it = pull.begin(); it != pull.end(); ++it)
      keys.push_back(it->second); 

#ifdef TIME_CONSTRUCT_ALL_CORRS
      snoop.stop(); 
      std::cout << " ** time to build xml for " <<  keys.size() 
        << " correlators " << snoop.getTimeInSeconds() << " seconds" << std::endl;
#endif

      return keys; 
    }

} // radmat


