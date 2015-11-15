/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name :

 * Purpose :

 * Creation Date : 25-09-2012

 * Last Modified : Fri Dec  7 15:29:50 2012

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/


// #include "radmat/load_data/radmat_database_interface.h"
#include "hadron/hadron_npart_npt_corr.h"
#include "ensem/ensem.h"
#include "radmat/utils/splash.h"
#include "radmat/utils/pow2assert.h"
#include "semble/semble_meta.h"

#include <sstream>
#include <string>
#include <iostream>

// using namespace radmat;
using namespace Hadron;
using namespace ENSEM;
// using namespace FILEDB;

typedef KeyHadronNPartNPtCorr_t KEY;
typedef EnsemVectorComplex DATA;




int main(int argc, char *argv[])
{

  if(argc != 2)
  {
    SPLASH("usage: test_load_database : <path/to/database>");
    exit(1);
  }
#if 0

  std::istringstream arg1(argv[1]);
  std::string dbfile;
  arg1 >> dbfile;


  RadmatDBInterface<KEY,DATA>::DBParams dbpars;
  dbpars.dbFile = dbfile;
  dbpars.useKeyFile = false;



  RadmatDBInterface<KEY,DATA> foobar(dbpars);

  // check that we can get stuff
  std::vector<KEY> foobarKeys;
  foobarKeys = foobar.getKeys();
  KEY firstKey, lastKey; 
  firstKey =  *(foobarKeys.begin());
  lastKey =  *--(foobarKeys.end());
  std::vector<KEY> test_keys;
  std::vector<DATA> test_data;
  test_keys.push_back(firstKey);
  test_keys.push_back(lastKey);
  test_data = foobar.getEnsem(test_keys);

  // look at it to make sure it makes sense
  std::cout << "first key \n" << firstKey << std::endl;
  std::cout << "mean data first key \n" << ENSEM::mean(foobar.getEnsem(firstKey)) << std::endl;
  std::cout << "last key \n" << lastKey << std::endl;
  std::cout << "mean data last key \n" << ENSEM::mean(foobar.getEnsem(lastKey)) << std::endl;


  // check that we can get the stuff either way and that 
  // the data is equal by element in the ensemble as it should be
  bool success = true;
  DATA alpha = foobar.getEnsem(firstKey);
  DATA omega = foobar.getEnsem(lastKey);

  const unsigned int size = alpha.size();
  const unsigned int nelem = alpha.numElem();
  for(unsigned int obs = 0; obs < nelem; ++obs)
  {
    for(unsigned int idx = 0; idx < size; ++idx)
    {
      success &= SEMBLE::toScalar(ENSEM::peekObs(alpha,obs).elem(idx))
        ==  SEMBLE::toScalar(ENSEM::peekObs(test_data[0],obs).elem(idx));
      success &= SEMBLE::toScalar(ENSEM::peekObs(omega,obs).elem(idx))
        ==  SEMBLE::toScalar(ENSEM::peekObs(test_data[1],obs).elem(idx));
    }
  }

  if(!!!success)
  {
    SPLASH("Error pulling out data");
  }

  POW2_ASSERT(success);

#endif
  return 0;
}
