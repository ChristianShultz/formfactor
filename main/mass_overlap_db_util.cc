/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : mass_overlap_db_combine.cc

 * Purpose :

 * Creation Date : 08-01-2013

 * Last Modified : Wed 30 Apr 2014 02:02:58 PM EDT

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/




/*
    !!NB: If the keys are automatically collapsing themselves the momenta printed will always be 0 0 0 and the row will always be 1
          as this allows us to only care about the momentum type (ie everything in the star of p has the same energy and overlap).
          that different rows of the little group have the same overlaps/energy  (controlled by a flag called LG_symmetry when you make the dbs
          and by a flag of the same(ish) name when running radmat. the key has a method to just rewrites some bit of itself on this flag)
  
   */






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
#include "radmat/utils/pow2assert.h"
#include "hadron/ensem_filenames.h"
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




struct dbInterface
{

  dbInterface(void) {}
  ~dbInterface(void) {if (m_db); delete m_db; m_db = NULL;}

  bool alloc(const std::string &dbname)
  {
    try
    {
      m_db = new FILEDB::AllConfStoreDB<SK,SD>();// (std::vector<int>(1,1));
    }
    catch(...)
    {
      std::cerr << __func__ << ": ERROR: couldn't alloc a database" << std::endl;
      delete m_db;
      m_db = NULL;
      exit(1);
    }

    if(m_db->open(dbname,  O_RDONLY, 400) != 0)
    {
      std::cerr << __func__ << ": error opening dbase= " <<dbname << std::endl;
      exit(1);
    }
    return true;
  }


  std::vector<K> getKeys(void)
  {
    POW2_ASSERT(&*m_db);
    std::vector<SK> keys;
    m_db->keys(keys); 
    std::vector<K> kk;
    std::vector<SK>::const_iterator it;
    for(it = keys.begin(); it != keys.end(); ++it)
      kk.push_back(it->key()); 

    return kk; 
  }

  D getData(const K &k)
  {
    SD dat;
    D eval;
    SK key;
    key.key() = k;
    std::vector<SD> val;
    int ret(0);

    if( (ret = m_db->get(key,val)) != 0) 
    {
      std::cerr << __func__ << ": key not found\n" << k;
      exit(1);
    }

    eval = val[0].data();

    return eval; 
  }


  ENSEM::EnsemReal getZ(const K &k)
  {
    D foobar = getData(k);
    return foobar.Z();
  }

  ENSEM::EnsemReal getE(const K &k)
  {
    D foobar = getData(k);
    return foobar.E(); 
  }


  void printEnsem(const K &k)
  {
    ENSEM::write(Hadron::ensemFileName(k.m_basic_key) + k.m_particle_id + std::string("_E"), getE(k)); 
    ENSEM::write(Hadron::ensemFileName(k.m_basic_key) + k.m_particle_id + std::string("_Z"), getZ(k)); 
  }


  void printKeysXML(const std::string & file)
  {
    std::vector<K> keys = getKeys(); 
    std::vector<K>::const_iterator it; 
    XMLFileWriter out(file); 
    push(out, "Keys"); 
    for(it = keys.begin(); it != keys.end(); ++it)
      write(out, "elem" , *it);

    pop(out);
    out.close(); 
  }

  void printEnsem(const std::string &keyfile)
  {

    try
    {
      XMLReader xml_in(keyfile);
      XMLArray::Array<K> keys;
      read(xml_in,"/Keys",keys);

      for(int i = 0; i < keys.size(); ++i)
        printEnsem(keys[i]);
    }
    catch(const std::string& e) 
    {
      std::cerr << __func__ << ": Caught Exception: " << e << std::endl;
      exit(1);
    }
    catch(std::exception& e) 
    {
      std::cerr << __func__ << ": Caught standard library exception: " << e.what() << std::endl;
      exit(1);
    }

  }

  void printKeys(void)
  {
    std::vector<K> keys = getKeys();
    std::vector<K>::const_iterator it;
    for(it = keys.begin(); it != keys.end(); ++it)
      std::cout << *it << std::endl; 
  }


  FILEDB::AllConfStoreDB< SK , SD > *m_db;
};


namespace
{
  int usage (int argc, char** argv)
  {
    std::cerr << "Usage: " << argv[0] << " <database> <operation> [<args> ...]" << std::endl;
    return -1;
  }
}




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

  std::string op = argv[2];

  if(op == "keys")
    db.printKeys();
  else if(op == "keysxml")
  {
    if(argc != 4)
      return usage(argc , argv);

    std::string file = argv[3];

    db.printKeysXML(file); 
  }
  else if(op == "get")
  {
    if(argc != 4)
      return usage(argc,argv);

    std::string file = argv[3];
    db.printEnsem(file);
  }
  else
  {
    std::cerr << __func__ << ": unrecognized operation" << std::endl;
    exit(1);
  }



  return 0;
}
