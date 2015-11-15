#ifndef CONSTRUCT_CORRELATORS_UTILS_H
#define CONSTRUCT_CORRELATORS_UTILS_H 


#include "construct_correlators_xml.h"
#include "data_tag_redstar_interface.h"
#include "radmat/database/database.h"
#include "lattice_multi_data_object.h"
#include "hadron/hadron_npart_npt_corr.h"
#include "radmat/utils/mink_qsq.h"
#include "ensem/ensem.h"
#include "adat/map_obj.h"
#include <vector>
#include <string>
#include <map> 
#include <utility>

namespace radmat
{

  std::map<double,std::vector<TaggedEnsemRedstarNPtBlock> > 
    sort_tagged_corrs_by_Q2(const std::vector<TaggedEnsemRedstarNPtBlock> &); 

  std::map<std::string,std::vector<TaggedEnsemRedstarNPtBlock> > 
    sort_tagged_corrs_by_Q2_and_rotation_group(const std::vector<TaggedEnsemRedstarNPtBlock> &, const bool mix_irreps=false); 

  // the database type we will be using 
  typedef radmatAllConfDatabaseInterface< Hadron::KeyHadronNPartNPtCorr_t,
          ENSEM::EnsemVectorComplex,
          RadmatExtendedKeyHadronNPartIrrep_t,
          RadmatMassOverlapData_t> DatabaseInterface_t;

  // NB hardwire of types from above in here
  struct DatabaseInterface_k
  {    
    DatabaseInterface_k(const Hadron::KeyHadronNPartNPtCorr_t &k ,
        const std::string &id_sink, const std::string &id_source)
      : npt(k)
    {
      sink = RadmatExtendedKeyHadronNPartIrrep_t(id_sink,k.npoint[1].irrep);
      source = RadmatExtendedKeyHadronNPartIrrep_t(id_source,k.npoint[3].irrep); 
      sink.doLG_symmetry();
      source.doLG_symmetry();
    }

    Hadron::KeyHadronNPartNPtCorr_t npt;
    RadmatExtendedKeyHadronNPartIrrep_t sink,source;  
  };  


  struct ConstructCorrsMatrixElement
  {
    ConstructCorrsMatrixElement(const ENSEM::EnsemVectorComplex &d, 
        const ThreePointDataTag &t,
        const bool s
        )
      : data(d) , tag(t) , success(s)
    {  }

    ENSEM::EnsemVectorComplex data; 
    ThreePointDataTag tag; 
    bool success; 
  };

  std::pair<bool,std::vector<ConstructCorrsMatrixElement> >
    build_correlators(
        const std::vector<TaggedEnsemRedstarNPtBlock> &corrs,
        const std::string &sink_id, 
        const std::string &source_id, 
        const ThreePointCorrXMLIni_t::RenormalizationProp &,
        const DatabaseInterface_t & );

  // memory soft version of above -- faster 
  std::pair<bool, rHandle<LLSQLatticeMultiData> >
    build_correlators_no_copy(
        const std::vector<TaggedEnsemRedstarNPtBlock> &corrs,
        const std::string &sink_id, 
        const std::string &source_id, 
        const ThreePointCorrXMLIni_t::RenormalizationProp &,
        const DatabaseInterface_t & );

  // only take the first N cfgs 
  std::pair<bool, rHandle<LLSQLatticeMultiData> > 
    build_correlators_no_copy_Ncfg_mean_fast(
        const std::vector<TaggedEnsemRedstarNPtBlock> &corrs,
        const std::string &sink_id, 
        const std::string &source_id, 
        const ThreePointCorrXMLIni_t::RenormalizationProp &Z_V,
        const DatabaseInterface_t &db,
        const int N);

  namespace BAD_DATA_REPO
  {
    void dump_bad_data(void); 
  }

} // radmat


#endif /* CONSTRUCT_CORRELATORS_UTILS_H */
