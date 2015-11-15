#ifndef BUILD_CORRELATORS_BAD_DATA_REPOSITORY_H
#define BUILD_CORRELATORS_BAD_DATA_REPOSITORY_H 

#include "radmat/database/database.h"
#include "hadron/hadron_npart_npt_corr.h"
#include <omp.h>
#include <vector>
#include <map>


/*
 *  Roughly, when we go through and generate the xml for redstar we need 
 *  to record the missed corrs so that we can accumulate a list of the 
 *  lattice matrix elements we need to perform the calculation
 *
 *  NB: this class should also be thread safe so that we can build the corrs 
 *  in parallel in order to gain a speed boost -- building correlators  
 *  accounts for roughly one third of the calculation time on my mac 
 */



namespace radmat
{

  struct BuildCorrsLocalBadDataRepo_t
  {
  
    // pretty obvious
    void
      insert(const int tid, 
          const std::vector<RadmatExtendedKeyHadronNPartIrrep_t> &v); 

    void 
      insert(const int tid,
          const std::vector<Hadron::KeyHadronNPartNPtCorr_t> &v);

    // collect the result
    std::vector<Hadron::KeyHadronNPartNPtCorr_t> 
      merge_bad_data_npt(void) const;

    std::vector<RadmatExtendedKeyHadronNPartIrrep_t> 
      merge_bad_data_norm(void) const;

    // dump the result
    void 
      dump_baddies(void) const;

    // the integer shall be the thread number to avoid collisions 
    std::map<int,std::vector<RadmatExtendedKeyHadronNPartIrrep_t> > norms;
    std::map<int,std::vector<Hadron::KeyHadronNPartNPtCorr_t> > tpc; 
  };


} // radmat


#endif /* BUILD_CORRELATORS_BAD_DATA_REPOSITORY_H */
