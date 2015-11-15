/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : mass_overlap_db_combine.cc

 * Purpose :

 * Creation Date : 08-01-2013

 * Last Modified : Wed 30 Apr 2014 02:01:57 PM EDT

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/


#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <exception>
#include <vector>
#include <utility>
#include <map>
#include "io/adat_io.h"
#include "io/adat_xmlio.h"
#include "io/key_val_db.h"
#include "AllConfStoreDB.h"
#include "ensem/ensem.h"
#include "radmat/database/database.h"
#include "hadron/hadron_npart_irrep.h"
#include "semble/semble_meta.h"


using namespace radmat;
using namespace Hadron;
using namespace ENSEM;
using namespace ADATXML;

typedef RadmatExtendedKeyHadronNPartIrrep_t K;
typedef RadmatMassOverlapData_t D;
typedef ADATIO::SerialDBKey<K> SK;
typedef ADATIO::SerialDBData<D> SD;

namespace
{
  int usage (int argc, char** argv)
  {
    std::cerr << "Usage: " << argv[0] << " <new dbasename> <first edb> [<second edb> ...]" << std::endl;
    return -1;
  }
}




struct dbInterface
{
  typedef RadmatExtendedKeyHadronNPartIrrep_t K;
  typedef RadmatMassOverlapData_t D;
  typedef ADATIO::SerialDBKey<K> SK;
  typedef ADATIO::SerialDBData<D> SD;

  dbInterface(void) {}
  ~dbInterface(void) {if (m_db); delete m_db; m_db = NULL;}

  bool alloc(const std::string &dbname)
  {
    try
    {
      m_db = new FILEDB::AllConfStoreDB<SK,SD>(std::vector<int>(1,1));
    }
    catch(...)
    {
      std::cerr << __func__ << ": ERROR: couldn't alloc a database" << std::endl;
      delete m_db;
      m_db = NULL;
      exit(1);
    }

    if(m_db->open(dbname,  O_RDWR | O_TRUNC | O_CREAT, 0664) != 0)
    {
      std::cerr << __func__ << ": error opening dbase= " <<dbname << std::endl;
      exit(1);
    }
    return true;
  }

  bool insert(const SK &k, const SD &d)
  {
    return (m_db->insert(k,std::vector<SD>(1,d)) == 0);
  }


  bool insert(const K &k, const D &d)
  {
    if(!!!m_db)
    {
      std::cerr << __func__ << ": ERROR: need to alloc a database first" << std::endl;
      exit(1);
    }

    SK key;
    SD data;
    key.key() = k; 
    data.data() = d; 

    return insert(key,data); 
  }

  FILEDB::AllConfStoreDB< SK , SD > *m_db;
};





  int
main (int argc, char** argv)
{
  // Parse the arguments
  if (argc < 3)
  {
    return usage(argc, argv);
  }

  dbInterface db;
  db.alloc(argv[1]); 


  for(int i = 2; i < argc; ++i)
  {

    FILEDB::AllConfStoreDB<SK,SD> tmp;
    if(tmp.open(argv[i], O_RDONLY, 0400) != 0)
    {
      std::cerr << __func__ << ": error opening db " << argv[i] << std::endl;
      exit(1); 
    }

    std::vector<SK> keys;
    std::vector<std::vector<SD> > data; 

    tmp.keysAndData(keys,data);
    bool success = true; 
    for(unsigned int idx = 0; idx < keys.size(); ++idx)
      success &= db.insert(keys[idx],data[idx][0]);

    if(!!!success)   
    {
      std::cerr << __func__ << ": error inserting data for db " << argv[i] << std::endl;
      exit(1); 
    }

  }

  return 0;
}
