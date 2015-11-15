/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : construct_correlators_utils.cc

 * Purpose :

 * Creation Date : 13-11-2013

 * Last Modified : Thu 24 Apr 2014 11:31:49 AM EDT

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/


#include "construct_correlators_utils.h"
#include "construct_correlators_bad_data_repository.h"
#include "radmat/data_representation/data_representation.h"
#include "semble/semble_meta.h"
#include "semble/semble_file_management.h"
#include "adat/adat_stopwatch.h"
#include "hadron/ensem_filenames.h"

#define DEBUG_MSG_OFF
#include "radmat/utils/debug_handler.h"

#define DO_TIMING_SORT_MAT_ELEMS_BY_Q2
#define DO_TIMING_CACHE_NORM_MAT_ELEMS
#define DO_TIMING_SUM_NORM_MAT_ELEMS

namespace radmat
{
  namespace BAD_DATA_REPO
  {
    BuildCorrsLocalBadDataRepo_t local_bad_data_repo; 
    void dump_bad_data(void)
    {
      local_bad_data_repo.dump_baddies(); 
    }
  }



  namespace
  {

    template<typename T> 
      struct ThreePtPropagationFactor
      {
        typename SEMBLE::PromoteEnsem<T>::Type operator()(const ENSEM::EnsemReal &E_sink,
            const typename SEMBLE::PromoteEnsem<T>::Type &Z_sink,
            const int t_sink,
            const int t_ins,
            const ENSEM::EnsemReal &E_source, 
            const typename SEMBLE::PromoteEnsem<T>::Type &Z_source, 
            const int t_source)
        {
          return ((ENSEM::conj(Z_source)*Z_sink/ (E_source * E_sink * SEMBLE::toScalar(4.)))
              * ENSEM::exp(-E_sink*(SEMBLE::toScalar(double(t_sink - t_ins))))
              * ENSEM::exp(-E_source*(SEMBLE::toScalar(double(t_ins - t_source))))
              );

        }

      };


    typedef ADAT::MapObject<Hadron::KeyHadronNPartNPtCorr_t,
            ENSEM::EnsemVectorComplex> singleThreadQ2NormalizedCorrCache;


    TaggedEnsemRedstarNPtBlock
      organizeCorrelatorSum(const TaggedEnsemRedstarNPtBlock &b)
      {
        typedef EnsemRedstarNPtBlock::ListObj_t ListObj_t; 
        typedef EnsemRedstarNPtBlock::Coeff_t Coeff_t;

        TaggedEnsemRedstarNPtBlock ret; 

        ret.data_tag = b.data_tag; 

        EnsemRedstarNPtBlock block = b.coeff_lattice_xml; 
        EnsemRedstarNPtBlock::const_iterator it; 
        std::map<std::string,ListObj_t> unique_elems; 
        std::map<std::string,ListObj_t>::iterator map_iterator; 

        for(it = block.begin(); it != block.end(); ++it)
        {
          std::string key = Hadron::ensemFileName(it->m_obj); 
          map_iterator = unique_elems.find(key); 

          if(map_iterator != unique_elems.end() )
          {
            ListObj_t foo = unique_elems[key]; 
            Coeff_t weight = foo.m_coeff + it->m_coeff; 
            foo.m_coeff = weight; 

            // update value
            if ( SEMBLE::toScalar(ENSEM::localNorm2(foo.m_coeff)) > 1e-6)
              map_iterator->second = foo;
            else
              unique_elems.erase(key);  
          }
          else
            unique_elems.insert(std::pair<std::string,ListObj_t>(key,*it)); 
        }

        std::map<std::string,ListObj_t>::const_iterator uniq; 
        std::stringstream coeff_obj_list ; 

        // construct the organized return xml list 
        // and record the list 
        for(uniq = unique_elems.begin(); uniq != unique_elems.end(); ++uniq)
        {
          ret.coeff_lattice_xml.push_back(uniq->second); 
          coeff_obj_list << SEMBLE::toScalar(uniq->second.m_coeff) << " * "; 
          coeff_obj_list << uniq->first << "\n";
        }

        std::string pth = SEMBLE::SEMBLEIO::getPath();
        std::stringstream path;

        // dump the ingredient list
        path << pth << "Q2_" << ret.qsq_tag(); 
        SEMBLE::SEMBLEIO::makeDirectoryPath(path.str()); 
        path << "/continuum_corr_ingredients";
        SEMBLE::SEMBLEIO::makeDirectoryPath(path.str());
        path << "/" << ret.data_tag.file_id;

        std::ofstream out(path.str().c_str()); 
        out << coeff_obj_list.str();
        out.close(); 

        return ret; 
      }


    ////////////////////////////////////////////////////////////////////

    // do the work and return the unique elems from the vector
    singleThreadQ2NormalizedCorrCache
      normalize_correlators( 
          const std::vector<TaggedEnsemRedstarNPtBlock> &corrs,
          std::vector<Hadron::KeyHadronNPartNPtCorr_t> &missed_xml, 
          std::vector<RadmatExtendedKeyHadronNPartIrrep_t> &missed_norm, 
          const std::string &sink_id, 
          const std::string &source_id, 
          const ThreePointCorrXMLIni_t::RenormalizationProp & Z_V,
          const DatabaseInterface_t &db,
          const int N=0
          )
      {
#ifdef DO_TIMING_CACHE_NORM_MAT_ELEMS
        Util::StopWatch snoop; 
        snoop.start();
#endif

        DEBUG_MSG(entering); 

        singleThreadQ2NormalizedCorrCache cache; 
        std::vector<TaggedEnsemRedstarNPtBlock>::const_iterator block; 
        EnsemRedstarNPtBlock::const_iterator npt; 

        // set up a zero
        ENSEM::EnsemVectorComplex zero; 
        zero.resize(1);
        zero.resizeObs(1);
        zero = SEMBLE::toScalar(0.); 

        for ( block = corrs.begin(); block != corrs.end(); ++block)
          for(npt = block->coeff_lattice_xml.begin(); npt != block->coeff_lattice_xml.end(); ++npt) 
          {
            if(cache.exist(npt->m_obj))
            {
              DEBUG_MSG(breaking early on success);
              continue; 
            }

            // can we make it?
            bool found = true; 
            DatabaseInterface_k k(npt->m_obj,sink_id,source_id);  
            if( !!!db.exists(npt->m_obj) ) 
            {
              missed_xml.push_back(npt->m_obj); 
              found = false;   
            }
            if( !!!db.exists(k.sink) ) 
            {
              missed_norm.push_back(k.sink); 
              found = false; 
            }
            if( !!!db.exists(k.source) ) 
            {
              missed_norm.push_back(k.source); 
              found = false; 
            }

            // break early if we're missing ingredients 
            if ( !!! found ) 
            {
              DEBUG_MSG(breaking out early b/c missing data);
              continue; 
            }

            // --- DO NORMALIZATION/RENORMALIZATION HERE
            // these are real but the * operator is only overloaded
            // for complex types so shove a zero down its throat
            // to get the desired result
            double Z = Z_V.RGE_prop;  
            ENSEM::Complex Z_c; 

            std::string g_rep = block->data_tag.origin_rep.g; 
            if ( g_rep == Stringify<J0p>() ) 
            {
              Z_c = SEMBLE::toScalar(std::complex<double>(Z*Z_V.Z_t)); 
            }
            else if ( g_rep == Stringify<J1m>() )
            {
              Z_c = SEMBLE::toScalar(std::complex<double>(Z*Z_V.Z_s)); 
            }
            else
            {
              std::cerr << __PRETTY_FUNCTION__ 
                << ": Error, I don't know what lorentz index " 
                << g_rep << " means " 
                << std::endl; 
              exit(1);  
            }

            // initialize some variables        
            ThreePtPropagationFactor<double> propagation_factor;

            ENSEM::EnsemVectorComplex corr_tmp = db.fetch(npt->m_obj);
            RadmatMassOverlapData_t source = db.fetch(k.source); 
            RadmatMassOverlapData_t sink = db.fetch(k.sink); 

            if( N != 0 )
            {
              ENSEM::EnsemVectorComplex corr_u; 
              ENSEM::EnsemReal src_E_u, src_Z_u;
              ENSEM::EnsemReal snk_E_u, snk_Z_u;
              int sz = N; 

              src_E_u.resize(sz); 
              src_Z_u.resize(sz); 
              snk_E_u.resize(sz); 
              snk_Z_u.resize(sz); 
              corr_u = corr_tmp;  
              corr_u.resize(sz); 

              for(int i = 0; i < sz; ++i)
              {
                ENSEM::pokeEnsem( src_E_u, ENSEM::peekEnsem(source.E(),i), i);
                ENSEM::pokeEnsem( src_Z_u, ENSEM::peekEnsem(source.Z(),i), i);

                ENSEM::pokeEnsem( snk_E_u, ENSEM::peekEnsem(sink.E(),i), i);
                ENSEM::pokeEnsem( snk_Z_u, ENSEM::peekEnsem(sink.Z(),i), i);

                ENSEM::pokeEnsem( corr_u, ENSEM::peekEnsem(corr_tmp,i), i); 
              }

              corr_tmp = corr_u; 
              source.E() = src_E_u; 
              source.Z() = src_Z_u; 
              sink.E() = snk_E_u; 
              sink.Z() = snk_Z_u; 
            }

            std::string outstem = Hadron::ensemFileName(npt->m_obj); 
            std::string pth = SEMBLE::SEMBLEIO::getPath();
            std::stringstream path;

            path << pth << "Q2_" << block->qsq_tag(); 
            SEMBLE::SEMBLEIO::makeDirectoryPath(path.str()); 

            path << "/correlator_normalization";
            SEMBLE::SEMBLEIO::makeDirectoryPath(path.str());

            path << "/" << outstem;
            ENSEM::write(path.str() + std::string("_corr_pre") , corr_tmp); 
            ENSEM::EnsemVectorComplex norm = corr_tmp * SEMBLE::toScalar(0.); 

            // the hadron key uses 1 based arrays
            // NB: assumption that npt is organized like <sink, ins , source>
            const int t_sink(npt->m_obj.npoint[1].t_slice); 
            const int t_source(npt->m_obj.npoint[3].t_slice);

            // sanity
            POW2_ASSERT(t_source < t_sink); 

            // NB: the indexing here assumes [tsource,tsink] ie: inclusive range
            for(int t_ins = t_source; t_ins <= t_sink; ++t_ins)
            {
              ENSEM::EnsemReal prop = propagation_factor(sink.E(),sink.Z(),t_sink,t_ins,
                  source.E(),source.Z(),t_source);

              ENSEM::pokeObs(corr_tmp,ENSEM::peekObs(corr_tmp,t_ins)/prop/Z_c,t_ins);

              ENSEM::pokeObs(norm,prop*Z_c,t_ins); 
            } // end loop over t_ins

            ENSEM::write(path.str() + std::string("_corr_post"), corr_tmp);
            ENSEM::write(path.str() + std::string("_norm") , norm);

            cache.insert(npt->m_obj,corr_tmp); 

          } // end loop over npt  

#ifdef DO_TIMING_CACHE_NORM_MAT_ELEMS
        snoop.stop(); 
        std::cout << __func__ << ": normalizing "
          << cache.size() << " unique npts took "
          << snoop.getTimeInSeconds() << " seconds" << std::endl; 
#endif

        return cache; 
      } 

  } // anonomyous 


  ////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////
  // brute force sort on qsq
  std::map<double,std::vector<TaggedEnsemRedstarNPtBlock> > 
    sort_tagged_corrs_by_Q2(const std::vector<TaggedEnsemRedstarNPtBlock> &unsorted)
    {

#ifdef DO_TIMING_SORT_MAT_ELEMS_BY_Q2
      Util::StopWatch snoop;
      snoop.start(); 
#endif

      std::map<double,std::vector<TaggedEnsemRedstarNPtBlock> > ret; 
      std::vector<TaggedEnsemRedstarNPtBlock>::const_iterator it; 

      // loop everything and toss it into a vector if it matches else make a new vector
      for(it = unsorted.begin(); it != unsorted.end(); ++it)
        if(ret.find(it->qsq_tag()) != ret.end())
          ret.find(it->qsq_tag())->second.push_back(*it);
        else
          ret.insert(std::map<double,std::vector<TaggedEnsemRedstarNPtBlock> >::value_type(
                it->qsq_tag(),std::vector<TaggedEnsemRedstarNPtBlock>(1,*it))); 


#ifdef DO_TIMING_SORT_MAT_ELEMS_BY_Q2
      snoop.stop(); 

      std::map<double,std::vector<TaggedEnsemRedstarNPtBlock> >::const_iterator it1; 
      int ct(0);
      for(it1 = ret.begin(); it1 != ret.end(); ++it1)
        ct += it1->second.size();

      std::cout << __func__ << ": sorting, " << ct << " elems took " 
        << snoop.getTimeInSeconds() << " seconds "  << std::endl; 
#endif 
      return ret;
    }

  ////////////////////////////////////////////////////////////////////
  // brute force sort on qsq
  std::map<std::string,std::vector<TaggedEnsemRedstarNPtBlock> > 
    sort_tagged_corrs_by_Q2_and_rotation_group(const std::vector<TaggedEnsemRedstarNPtBlock> &unsorted, 
        const bool mix_irreps)
    {

#ifdef DO_TIMING_SORT_MAT_ELEMS_BY_Q2
      Util::StopWatch snoop;
      snoop.start(); 
#endif

      std::cout << __func__ << ": sorting correlators"; 
      if( mix_irreps )
        std::cout << " using mix_irreps=true" << std::endl;
      else
        std::cout << " using mix_irreps=false" << std::endl;

        
      std::map<std::string,std::vector<TaggedEnsemRedstarNPtBlock> > ret; 
      std::vector<TaggedEnsemRedstarNPtBlock>::const_iterator it; 

      // loop everything and toss it into a vector if it matches else make a new vector
      for(it = unsorted.begin(); it != unsorted.end(); ++it)
        if(ret.find(it->rot_qsq_tag(mix_irreps)) != ret.end())
          ret.find(it->rot_qsq_tag(mix_irreps))->second.push_back(*it);
        else
          ret.insert(std::map<std::string,std::vector<TaggedEnsemRedstarNPtBlock> >::value_type(
                it->rot_qsq_tag(mix_irreps),std::vector<TaggedEnsemRedstarNPtBlock>(1,*it))); 


#ifdef DO_TIMING_SORT_MAT_ELEMS_BY_Q2
      snoop.stop(); 

      std::map<std::string,std::vector<TaggedEnsemRedstarNPtBlock> >::const_iterator it1; 
      int ct(0);
      for(it1 = ret.begin(); it1 != ret.end(); ++it1)
        ct += it1->second.size();

      std::cout << __func__ << ": sorting, " << ct << " elems " 
        << " into " << ret.size() << " rotation groups took " 
        << snoop.getTimeInSeconds() << " seconds "  << std::endl; 
#endif 
      return ret;
    }


  ////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////
  //build a correlator
  std::pair<bool, std::vector<ConstructCorrsMatrixElement> > 
    build_correlators(
        const std::vector<TaggedEnsemRedstarNPtBlock> &corrs,
        const std::string &sink_id, 
        const std::string &source_id, 
        const ThreePointCorrXMLIni_t::RenormalizationProp &Z_V,
        const DatabaseInterface_t &db)
    {
#ifdef DO_TIMING_SUM_NORM_MAT_ELEMS
      Util::StopWatch snoop;
      snoop.start(); 
#endif

      std::vector<Hadron::KeyHadronNPartNPtCorr_t> missed_xml; 
      std::vector<RadmatExtendedKeyHadronNPartIrrep_t> missed_norm; 

      std::vector<ConstructCorrsMatrixElement> ret; 
      double qsq; 

      // create a cache of all ingredients and 
      // strip all members of principle time 
      // dependence and overlap factors, 
      // ie: isolate the cubic matrix element 
      singleThreadQ2NormalizedCorrCache npoint_cache; 
      npoint_cache = normalize_correlators( corrs, missed_xml, missed_norm, 
          sink_id, source_id, Z_V, db);

      // push the bad data at the global repo 
      BAD_DATA_REPO::local_bad_data_repo.insert(omp_get_thread_num(),missed_xml);   
      BAD_DATA_REPO::local_bad_data_repo.insert(omp_get_thread_num(),missed_norm);   

      // abort if we missed everything
      if(npoint_cache.size() == 0) 
      {
        DEBUG_MSG(exiting early);
#ifdef DO_TIMING_SUM_NORM_MAT_ELEMS
        snoop.stop(); 
        std::cout << __func__ << "npoint_cache is empty"
          << " exiting early " << std::endl; 
#endif
        return std::pair<bool, std::vector<ConstructCorrsMatrixElement> >(false,ret); 
      }


      // intitialize some variables
      ENSEM::EnsemReal ScalarZeroR;
      ENSEM::EnsemVectorComplex VectorZeroC; 
      bool any_data = false; 
      std::vector<TaggedEnsemRedstarNPtBlock>::const_iterator block_it; 
      EnsemRedstarNPtBlock::const_iterator npt; 
      bool found = false; 

      // set the above variables (determine things like nbins, Lt etc)
      for ( block_it = corrs.begin(); block_it != corrs.end(); ++block_it)
      {
        if (found ) 
          break; 

        for(npt = block_it->coeff_lattice_xml.begin(); npt != block_it->coeff_lattice_xml.end(); ++npt) 
          if(npoint_cache.exist(npt->m_obj))
          {
            found = true; 
            DatabaseInterface_k k(npt->m_obj,sink_id,source_id);  
            RadmatMassOverlapData_t dummy = db.fetch(k.source); 
            ScalarZeroR = SEMBLE::toScalar(0.) * dummy.E();  
            VectorZeroC = SEMBLE::toScalar(std::complex<double>(0.,0.))* db.fetch(npt->m_obj); 
            break; 
          }
      }     

      // can't do anything -- no data 
      if ( !!! found ) 
        return std::pair<bool, std::vector<ConstructCorrsMatrixElement> >(false,ret);  

      for ( block_it = corrs.begin(); block_it != corrs.end(); ++block_it)
      {
        // A + B + C - B + A -> 2A +C
        TaggedEnsemRedstarNPtBlock block = organizeCorrelatorSum(*block_it); 

        // return corrs
        ENSEM::EnsemVectorComplex corr; 
        ENSEM::EnsemReal E_snk, E_src; 

        // zero out data
        corr = VectorZeroC; 
        E_snk = ScalarZeroR; 
        E_src = ScalarZeroR;
        int ct(0); 
        bool success = true;  

        for(npt = block.coeff_lattice_xml.begin(); npt != block.coeff_lattice_xml.end(); ++npt) 
        {
          success &= npoint_cache.exist(npt->m_obj);

          // abandon rest of loop if we miss one
          if ( !!! success ) 
            break; 

          // pull out data 
          DatabaseInterface_k k(npt->m_obj,sink_id,source_id);  
          RadmatMassOverlapData_t source = db.fetch(k.source); 
          RadmatMassOverlapData_t sink = db.fetch(k.sink); 

          // update energies
          E_snk = E_snk + sink.E(); 
          E_src = E_src + source.E(); 

          // update correlator sum 
          corr = corr + npt->m_coeff * npoint_cache[npt->m_obj];

          ++ct; 
        }

        // the return energy is a flat average 
        // of the things that went into the correlator
        if ( ct != 0)
        {
          E_snk = E_snk / SEMBLE::toScalar(double(ct)); 
          E_src = E_src / SEMBLE::toScalar(double(ct)); 
        }

        // thingy 
        ThreePointDataTag tag = block.data_tag;
        tag.left_E = E_snk; 
        tag.right_E = E_src; 

        any_data |= success; 

        // only push back good corrs
        if(success) 
          ret.push_back( ConstructCorrsMatrixElement(corr,tag,success) );
      }

#ifdef DO_TIMING_SUM_NORM_MAT_ELEMS
      snoop.stop(); 
      std::cout << __func__ << ": building " << ret.size() 
        << " matrix elements took " << snoop.getTimeInSeconds() 
        << " seconds" << std::endl; 
#endif

      return std::pair<bool, std::vector<ConstructCorrsMatrixElement> >(any_data,ret);  
    }


  ////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////
  //build a correlator
  std::pair<bool, rHandle<LLSQLatticeMultiData> > 
    build_correlators_no_copy_Ncfg_mean_fast(
        const std::vector<TaggedEnsemRedstarNPtBlock> &corrs,
        const std::string &sink_id, 
        const std::string &source_id, 
        const ThreePointCorrXMLIni_t::RenormalizationProp &Z_V,
        const DatabaseInterface_t &db,
        const int N)
    {
#ifdef DO_TIMING_SUM_NORM_MAT_ELEMS
      Util::StopWatch snoop;
      snoop.start(); 
#endif

      std::vector<Hadron::KeyHadronNPartNPtCorr_t> missed_xml; 
      std::vector<RadmatExtendedKeyHadronNPartIrrep_t> missed_norm; 

      rHandle<LLSQLatticeMultiData> ret(new LLSQLatticeMultiData) ; 
      double qsq; 

      // create a cache of all ingredients and 
      // strip all members of principle time 
      // dependence and overlap factors, 
      // ie: isolate the cubic matrix element 
      singleThreadQ2NormalizedCorrCache npoint_cache; 
      npoint_cache = normalize_correlators( corrs, missed_xml, missed_norm, 
          sink_id, source_id, Z_V, db, N);

      // push the bad data at the global repo 
      BAD_DATA_REPO::local_bad_data_repo.insert(omp_get_thread_num(),missed_xml);   
      BAD_DATA_REPO::local_bad_data_repo.insert(omp_get_thread_num(),missed_norm);   

      // abort if we missed everything
      if(npoint_cache.size() == 0) 
      {
        DEBUG_MSG(exiting early);
#ifdef DO_TIMING_SUM_NORM_MAT_ELEMS
        snoop.stop(); 
        std::cout << __func__ << "npoint_cache is empty"
          << " exiting early " << std::endl; 
#endif
        return std::pair<bool, rHandle<LLSQLatticeMultiData> >(false,ret); 
      }


      // intitialize some variables
      ENSEM::EnsemReal ScalarZeroR;
      ENSEM::EnsemVectorComplex VectorZeroC; 
      bool any_data = false; 
      std::vector<TaggedEnsemRedstarNPtBlock>::const_iterator block_it; 
      EnsemRedstarNPtBlock::const_iterator npt; 
      bool found = false; 

      // set the above variables (determine things like nbins, Lt etc)
      for ( block_it = corrs.begin(); block_it != corrs.end(); ++block_it)
      {
        if (found ) 
          break; 

        for(npt = block_it->coeff_lattice_xml.begin(); npt != block_it->coeff_lattice_xml.end(); ++npt) 
          if(npoint_cache.exist(npt->m_obj))
          {
            found = true; 
            DatabaseInterface_k k(npt->m_obj,sink_id,source_id);  
            RadmatMassOverlapData_t dummy = db.fetch(k.source); 
            ScalarZeroR = SEMBLE::toScalar(0.) * dummy.E();  
            VectorZeroC = SEMBLE::toScalar(std::complex<double>(0.,0.))* npoint_cache[npt->m_obj]; 
            break; 
          }
      }     

      // can't do anything -- no data 
      if ( !!! found ) 
        return std::pair<bool, rHandle<LLSQLatticeMultiData> >(false,ret);  

      for ( block_it = corrs.begin(); block_it != corrs.end(); ++block_it)
      {
        // A + B + C - B + A -> 2A +C
        TaggedEnsemRedstarNPtBlock block = organizeCorrelatorSum(*block_it); 

        // return corrs
        ENSEM::EnsemVectorComplex corr; 
        ENSEM::EnsemReal E_snk, E_src; 

        // zero out data
        corr = VectorZeroC; 
        E_snk = ScalarZeroR; 
        E_src = ScalarZeroR;
        int ct(0); 
        bool success = true;  

        for(npt = block.coeff_lattice_xml.begin(); npt != block.coeff_lattice_xml.end(); ++npt) 
        {
          success &= npoint_cache.exist(npt->m_obj);

          // abandon rest of loop if we miss one
          if ( !!! success ) 
            break; 

          // pull out data 
          DatabaseInterface_k k(npt->m_obj,sink_id,source_id);  
          RadmatMassOverlapData_t source = db.fetch(k.source); 
          RadmatMassOverlapData_t sink = db.fetch(k.sink); 

          // update energies
          E_snk = E_snk + sink.E(); 
          E_src = E_src + source.E(); 

          // update correlator sum 
          corr = corr + npt->m_coeff * npoint_cache[npt->m_obj];

          ++ct; 
        }

        // the return energy is a flat average 
        // of the things that went into the correlator
        if ( ct != 0)
        {
          E_snk = E_snk / SEMBLE::toScalar(double(ct)); 
          E_src = E_src / SEMBLE::toScalar(double(ct)); 
        }

        // hack time !! 
        ENSEM::EnsemReal l_sE,r_sE; 
        ENSEM::EnsemVectorComplex corr_s; 
        int sz = N; 
        l_sE.resize(sz); 
        r_sE.resize(sz); 
        corr_s = corr; 
        corr_s.resize(sz);

        for(int i =0; i < sz; ++i)
        {
          ENSEM::pokeEnsem(l_sE,ENSEM::peekEnsem(E_snk,i),i);
          ENSEM::pokeEnsem(r_sE,ENSEM::peekEnsem(E_src,i),i);
          ENSEM::pokeEnsem(corr_s,ENSEM::peekEnsem(corr,i),i); 
        }
      

        // thingy 
        ThreePointDataTag tag = block.data_tag;
        tag.left_E = l_sE; 
        tag.right_E = r_sE; 

        any_data |= success; 

        // only push back good corrs
        if(success) 
          ret->append_row_ensem(corr_s,tag); 
      }

#ifdef DO_TIMING_SUM_NORM_MAT_ELEMS
      snoop.stop(); 
      std::cout << __func__ << ": building " << ret->nrows() 
        << " matrix elements took " << snoop.getTimeInSeconds() 
        << " seconds" << std::endl; 
#endif

      return std::pair<bool, rHandle<LLSQLatticeMultiData> >(any_data,ret);  
    }


  ////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////
  //build a correlator
  std::pair<bool, rHandle<LLSQLatticeMultiData> > 
    build_correlators_no_copy(
        const std::vector<TaggedEnsemRedstarNPtBlock> &corrs,
        const std::string &sink_id, 
        const std::string &source_id, 
        const ThreePointCorrXMLIni_t::RenormalizationProp &Z_V,
        const DatabaseInterface_t &db)
    {
#ifdef DO_TIMING_SUM_NORM_MAT_ELEMS
      Util::StopWatch snoop;
      snoop.start(); 
#endif

      std::vector<Hadron::KeyHadronNPartNPtCorr_t> missed_xml; 
      std::vector<RadmatExtendedKeyHadronNPartIrrep_t> missed_norm; 

      rHandle<LLSQLatticeMultiData> ret(new LLSQLatticeMultiData) ; 
      double qsq; 

      // create a cache of all ingredients and 
      // strip all members of principle time 
      // dependence and overlap factors, 
      // ie: isolate the cubic matrix element 
      singleThreadQ2NormalizedCorrCache npoint_cache; 
      npoint_cache = normalize_correlators( corrs, missed_xml, missed_norm, 
          sink_id, source_id, Z_V, db);

      // push the bad data at the global repo 
      BAD_DATA_REPO::local_bad_data_repo.insert(omp_get_thread_num(),missed_xml);   
      BAD_DATA_REPO::local_bad_data_repo.insert(omp_get_thread_num(),missed_norm);   

      // abort if we missed everything
      if(npoint_cache.size() == 0) 
      {
        DEBUG_MSG(exiting early);
#ifdef DO_TIMING_SUM_NORM_MAT_ELEMS
        snoop.stop(); 
        std::cout << __func__ << "npoint_cache is empty"
          << " exiting early " << std::endl; 
#endif
        return std::pair<bool, rHandle<LLSQLatticeMultiData> >(false,ret); 
      }


      // intitialize some variables
      ENSEM::EnsemReal ScalarZeroR;
      ENSEM::EnsemVectorComplex VectorZeroC; 
      bool any_data = false; 
      std::vector<TaggedEnsemRedstarNPtBlock>::const_iterator block_it; 
      EnsemRedstarNPtBlock::const_iterator npt; 
      bool found = false; 

      // set the above variables (determine things like nbins, Lt etc)
      for ( block_it = corrs.begin(); block_it != corrs.end(); ++block_it)
      {
        if (found ) 
          break; 

        for(npt = block_it->coeff_lattice_xml.begin(); npt != block_it->coeff_lattice_xml.end(); ++npt) 
          if(npoint_cache.exist(npt->m_obj))
          {
            found = true; 
            DatabaseInterface_k k(npt->m_obj,sink_id,source_id);  
            RadmatMassOverlapData_t dummy = db.fetch(k.source); 
            ScalarZeroR = SEMBLE::toScalar(0.) * dummy.E();  
            VectorZeroC = SEMBLE::toScalar(std::complex<double>(0.,0.))* db.fetch(npt->m_obj); 
            break; 
          }
      }     

      // can't do anything -- no data 
      if ( !!! found ) 
        return std::pair<bool, rHandle<LLSQLatticeMultiData> >(false,ret);  

      for ( block_it = corrs.begin(); block_it != corrs.end(); ++block_it)
      {
        // A + B + C - B + A -> 2A +C
        TaggedEnsemRedstarNPtBlock block = organizeCorrelatorSum(*block_it); 

        // return corrs
        ENSEM::EnsemVectorComplex corr; 
        ENSEM::EnsemReal E_snk, E_src; 

        // zero out data
        corr = VectorZeroC; 
        E_snk = ScalarZeroR; 
        E_src = ScalarZeroR;
        int ct(0); 
        bool success = true;  

        for(npt = block.coeff_lattice_xml.begin(); npt != block.coeff_lattice_xml.end(); ++npt) 
        {
          success &= npoint_cache.exist(npt->m_obj);

          // abandon rest of loop if we miss one
          if ( !!! success ) 
            break; 

          // pull out data 
          DatabaseInterface_k k(npt->m_obj,sink_id,source_id);  
          RadmatMassOverlapData_t source = db.fetch(k.source); 
          RadmatMassOverlapData_t sink = db.fetch(k.sink); 

          // update energies
          E_snk = E_snk + sink.E(); 
          E_src = E_src + source.E(); 

          // update correlator sum 
          corr = corr + npt->m_coeff * npoint_cache[npt->m_obj];

          ++ct; 
        }

        // the return energy is a flat average 
        // of the things that went into the correlator
        if ( ct != 0)
        {
          E_snk = E_snk / SEMBLE::toScalar(double(ct)); 
          E_src = E_src / SEMBLE::toScalar(double(ct)); 
        }

        // thingy 
        ThreePointDataTag tag = block.data_tag;
        tag.left_E = E_snk; 
        tag.right_E = E_src; 

        any_data |= success; 

        // only push back good corrs
        if(success) 
          ret->append_row_ensem(corr,tag); 
      }

#ifdef DO_TIMING_SUM_NORM_MAT_ELEMS
      snoop.stop(); 
      std::cout << __func__ << ": building " << ret->nrows() 
        << " matrix elements took " << snoop.getTimeInSeconds() 
        << " seconds" << std::endl; 
#endif

      return std::pair<bool, rHandle<LLSQLatticeMultiData> >(any_data,ret);  
    }

} // radmat 

