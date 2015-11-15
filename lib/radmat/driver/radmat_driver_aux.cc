/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : radmat_driver_aux.cc

 * Purpose :

 * Creation Date : 25-06-2013

 * Last Modified : Thu Jun 27 08:13:27 2013

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/


#include "hadron/hadron_npart_npt_corr.h"
#include "hadron/hadron_npart_irrep.h"
#include "hadron/particle_op.h"
#include <vector>
#include <map>

namespace radmat
{
  namespace
  {
    void stubify(Hadron::KeyHadronNPartIrrep_t &rep)
    {
      rep.row = 1;

      /* -- redstar will fail if we mess with momentum
      rep.mom[0] = 0;
      rep.mom[1] = 0; 
      rep.mom[2] = 0; 
      */
    }

    void stubify(Hadron::KeyHadronNPartNPtCorr_t::NPoint_t &npt)
    {
      stubify(npt.irrep);
    }

    void stubify(Hadron::KeyHadronNPartNPtCorr_t &npt)
    {
      if(npt.npoint.size() != 3)
      {
        std::cerr << __func__ << ": error , don't know what to do, stop messing up" << std::endl;
        exit(1337); 
      }

      stubify(npt.npoint[1]);
      stubify(npt.npoint[3]);
    }

  } // anonomyous 


  void stubify(std::vector<Hadron::KeyHadronNPartNPtCorr_t> &keys)
  {
    std::vector<Hadron::KeyHadronNPartNPtCorr_t>::iterator it; 
    for(it = keys.begin(); it != keys.end(); ++it)
      stubify(*it); 
  }

}
