#ifndef RADMAT_DATABASE_INTERFACE_H
#define RADMAT_DATABASE_INTERFACE_H


#include <vector>
#include <exception>
#include "ensem/ensem.h"
#include "radmat/utils/handle.h"
#include "AllConfStoreDB.h"
#include "AllConfStoreMultipleDB.h"
#include "io/adat_io.h"
#include "io/adat_xmlio.h"
#include "io/key_val_db.h"


// #define DEBUGGING_RADMAT_DB_INTERFACE

#ifdef DEBUGGING_RADMAT_DB_INTERFACE
#include "semble/semble_meta.h" 
#endif

// #define DEBUGGING_ADATIO_KEYS

namespace radmat
{

  struct dbProp_t
  {
    std::vector<std::string> dbname;
    std::string badlist;
  };


  std::string toString(const dbProp_t &);
  std::ostream& operator<<(std::ostream&,const dbProp_t&);
  void read(ADATXML::XMLReader &xml, const std::string &path, dbProp_t &);
  void write(ADATXML::XMLWriter &xml, const std::string &path, const dbProp_t &);


  struct radmatDBProp_t
  {
    dbProp_t threePointDatabase;
    dbProp_t normalizationDatabase;
    bool allow_daggering;          // are the daggered normalizations and energies the same?  ie; is the data real
    bool LG_symmetry;              // does everything in the star of p have the same value AND row independent 
  };

  std::string toString(const radmatDBProp_t &);
  std::ostream& operator<<(std::ostream&, const radmatDBProp_t &);
  void read(ADATXML::XMLReader &xml, const std::string &path, radmatDBProp_t &);
  void write(ADATXML::XMLWriter &xml, const std::string &path, const radmatDBProp_t &);


  // some templated structure to deal with all of the database nastyness 
  template<typename CORRKEY, typename CORRDATA, typename NORMKEY, typename NORMDATA>
    struct radmatAllConfDatabaseInterface
    {
      typedef typename ADATIO::SerialDBKey<CORRKEY> S_C_KEY;
      typedef typename ENSEM::EnsemScalar<CORRDATA>::Type_t ESCALARDATA;
      typedef typename ADATIO::SerialDBData<ESCALARDATA> S_C_DATA;
      typedef typename ADATIO::SerialDBKey<NORMKEY> S_N_KEY;
      typedef typename ADATIO::SerialDBData<NORMDATA> S_N_DATA; 

      radmatAllConfDatabaseInterface(void); // hide ctor
      radmatAllConfDatabaseInterface(const radmatDBProp_t &);
      ~radmatAllConfDatabaseInterface(void) {outlog.close();}


      void alloc(void);
      bool exists(const CORRKEY &) const;
      bool exists(const NORMKEY &) const;
      CORRDATA fetch(const CORRKEY &) const;
      NORMDATA fetch(const NORMKEY &) const;

      private: 
      NORMDATA fetch_dagger(const NORMKEY &) const;  

      public:
      rHandle<FILEDB::AllConfStoreMultipleDB<S_C_KEY,S_C_DATA> > m_corr_db;
      rHandle<FILEDB::AllConfStoreMultipleDB<S_N_KEY,S_N_DATA> > m_norm_db;  
      radmatDBProp_t db_props;

      mutable std::ofstream outlog; // this is naughty but i wanted to move all of the 
      // output from the screen to a log w/o recoding any of the const methods.. 
    };

  template<typename CORRKEY, typename CORRDATA, typename NORMKEY, typename NORMDATA >
    radmatAllConfDatabaseInterface<CORRKEY,CORRDATA,NORMKEY,NORMDATA>::radmatAllConfDatabaseInterface(const radmatDBProp_t &prop)
    : db_props(prop)
    {
      alloc();
      outlog.open("database_search_log.txt", std::ios::out | std::ios::app);
    }


  // set up some handles and check that we can open it
  template<typename CORRKEY, typename CORRDATA, typename NORMKEY, typename NORMDATA>
    void radmatAllConfDatabaseInterface<CORRKEY,CORRDATA,NORMKEY,NORMDATA>::alloc(void)
    {

      m_corr_db = rHandle<FILEDB::AllConfStoreMultipleDB<S_C_KEY,S_C_DATA> > (new FILEDB::AllConfStoreMultipleDB<S_C_KEY,S_C_DATA>() );
      m_norm_db = rHandle<FILEDB::AllConfStoreMultipleDB<S_N_KEY,S_N_DATA> > (new FILEDB::AllConfStoreMultipleDB<S_N_KEY,S_N_DATA>() );

      m_corr_db->setCacheSize(1*1024*1024*1024); 
      m_norm_db->setCacheSize(1*1024*1024*1024); 

      // if(m_corr_db->open(db_props.threePointDatabase.dbname, O_RDONLY, 0400) != 0)
      if(m_corr_db->open(db_props.threePointDatabase.dbname) != 0)
      {
        std::vector<std::string>::const_iterator it; 
        outlog << __func__ << ": error opening one of the following" << std::endl;
        std::cerr << __func__ << ": error opening one of the following" << std::endl;
        for (it = db_props.threePointDatabase.dbname.begin(); it != db_props.threePointDatabase.dbname.end(); ++it)
        {
          outlog << *it << std::endl;
          std::cout << *it << std::endl;
        }
        exit(1);
      }


      // if(m_norm_db->open(db_props.normalizationDatabase.dbname, O_RDONLY, 0400) != 0)
      if(m_norm_db->open(db_props.normalizationDatabase.dbname) != 0)
      {
        std::vector<std::string>::const_iterator it; 
        outlog << __func__ << ": error opening one of the following" << std::endl;
        std::cerr << __func__ << ": error opening one of the following" << std::endl;
        for (it = db_props.normalizationDatabase.dbname.begin(); it != db_props.normalizationDatabase.dbname.end(); ++it)
        {
          outlog << *it << std::endl;
          std::cout << *it << std::endl;
        }
        exit(1);
      }
    }

  template<typename CORRKEY, typename CORRDATA, typename NORMKEY, typename NORMDATA>
    bool radmatAllConfDatabaseInterface<CORRKEY,CORRDATA,NORMKEY,NORMDATA>::exists(const CORRKEY &k) const
    {
      S_C_KEY key;
      key.key() = k;

      if( m_corr_db->exist(key) == 0)
      {
        outlog << __func__ << ": MISSED CORRELATOR KEY!!!" << std::endl;
        outlog << key.key() << std::endl;
      }

      return m_corr_db->exist(key);     
    }

  template<typename CORRKEY, typename CORRDATA, typename NORMKEY, typename NORMDATA>
    bool radmatAllConfDatabaseInterface<CORRKEY,CORRDATA,NORMKEY,NORMDATA>::exists(const NORMKEY &k) const
    {
      S_N_KEY key;
      key.key() = k;

      if(db_props.LG_symmetry)
        key.key().doLG_symmetry(); 

      if(!!! m_norm_db->exist(key))
      {
        if(!!!db_props.allow_daggering)
          return false;

        NORMKEY dag(k);
        dag.dagger();
        key.key() = dag;

        if(db_props.LG_symmetry)
          key.key().doLG_symmetry(); 


        if(!!! m_norm_db->exist(key))
        {
          outlog << __func__ << ": MISSED DAGGER TOO!!!!" << std::endl;
          outlog << key.key() << std::endl;
          return false; 
        }
      }

      return true; 
    }

  template<typename CORRKEY, typename CORRDATA, typename NORMKEY, typename NORMDATA>
    CORRDATA  radmatAllConfDatabaseInterface<CORRKEY,CORRDATA,NORMKEY,NORMDATA>::fetch(const CORRKEY &k) const
    {
      CORRDATA eval;


#ifdef DEBUGGING_ADATIO_KEYS
      std::cout << __func__ << ": trying to extract key " << k 
        << "\n" << __func__ 
        << ": trying to allocacte a key of type ADATIO::SerialDBKey<CORRKEY> " 
        << std::endl; 
#endif

      S_C_KEY key;

#ifdef DEBUGGING_ADATIO_KEYS
      std::cout << __func__ << ": trying to apply the equality operator " 
        << std::endl; 
#endif

      key.key() = k;
      std::vector<S_C_DATA> vals;
      int ret(0);

      try
      {
        if ((ret = m_corr_db->get(key, vals)) != 0)
        {
          outlog << __func__ << ": key not found\n" << k;
          std::cerr << __func__ << ": key not found\n" << k;
          exit(1);
        }

        eval.resize(vals.size());
        eval.resizeObs(vals[0].data().numElem());

        for(unsigned int i=0; i < vals.size(); ++i)
        {
          ESCALARDATA sval = vals[i].data();
          ENSEM::pokeEnsem(eval, sval, i);
        }

      }
      catch(const std::string& e) 
      {
        outlog << __func__ << ": Caught Exception: " << e << std::endl;
        std::cerr << __func__ << ": Caught Exception: " << e << std::endl; 
        exit(1);
      }
      catch(std::exception& e) 
      {
        outlog << __func__ << ": Caught standard library exception: " << e.what() << std::endl;
        std::cerr << __func__ << ": Caught standard library exception: " << e.what() << std::endl;
        std::cerr << __func__ << ": key -- " << key.key() << "\n vals.size = " << vals.size() 
          << "\n ret = " << ret << " eval.size = " << eval.size() << std::endl;
        exit(1);
      }

      return eval;
    }




  template<typename CORRKEY, typename CORRDATA, typename NORMKEY, typename NORMDATA>
    NORMDATA  radmatAllConfDatabaseInterface<CORRKEY,CORRDATA,NORMKEY,NORMDATA>::fetch(const NORMKEY &k) const
    {
      NORMDATA eval;
      try
      {
        S_N_KEY key;
        key.key() = k;

        if(db_props.LG_symmetry)
          key.key().doLG_symmetry(); 


        std::vector<S_N_DATA> val; // hacky thingy we did
        int ret(0);
        if ((ret = m_norm_db->get(key, val)) != 0)
        {
          outlog << __func__ << ": key not found\n" << k << std::endl;

          if(db_props.allow_daggering)
          {
            NORMKEY kdag(k);
            kdag.dagger(); 
            outlog << __func__ << ": looking for the daggered key\n" << kdag << std::endl;
            return fetch_dagger(kdag); 
          }
          else
          {
            std::cerr <<  __func__ << ": key not found\n" << k << std::endl;
            exit(1);
          }
        }

        eval = val[0].data(); // if this worked right we stored the entire thing in a single ensemble slot.. 

      }
      catch(const std::string& e) 
      {
        outlog << __func__ << ": Caught Exception: " << e << std::endl;
        std::cerr << __func__ << ": Caught Exception: " << e << std::endl; 
        exit(1);
      }
      catch(std::exception& e) 
      {
        outlog << __func__ << ": Caught standard library exception: " << e.what() << std::endl;
        std::cerr << __func__ << ": Caught standard library exception: " << e.what() << std::endl;
        exit(1);
      }

#ifdef DEBUGGING_RADMAT_DB_INTERFACE
      std::cout << __PRETTY_FUNCTION__ << std::endl;
      std::cout << "KEY ----- " << k << std::endl;
      std::cout << "ENERGY -- " << SEMBLE::toScalar(ENSEM::mean(eval.E())) << std::endl;
      std::cout << "       -- size " << eval.E().size() << std::endl;
      std::cout << "OVERLAP - " << SEMBLE::toScalar(ENSEM::mean(eval.Z())) << std::endl;
      std::cout << "       -- size " << eval.Z().size() << std::endl;
#endif

      return eval;
    } 

  // apply method dagger to key type and look again
  template<typename CORRKEY, typename CORRDATA, typename NORMKEY, typename NORMDATA>
    NORMDATA  radmatAllConfDatabaseInterface<CORRKEY,CORRDATA,NORMKEY,NORMDATA>::fetch_dagger(const NORMKEY &k) const
    {
      NORMDATA eval;
      try
      {
        S_N_KEY key;
        key.key() = k;

        if(db_props.LG_symmetry)
          key.key().doLG_symmetry(); 


        std::vector<S_N_DATA> val; // hacky thingy we did
        int ret(0);
        if ((ret = m_norm_db->get(key, val)) != 0)
        {
          outlog << __func__ << ": key not found\n" << k;
          exit(1);
        }

        eval = val[0].data(); // if this worked right we stored the entire thing in a single ensemble slot.. 

      }
      catch(const std::string& e) 
      {
        outlog << __func__ << ": Caught Exception: " << e << std::endl;
        std::cerr << __func__ << ": Caught Exception: " << e << std::endl;
        exit(1);
      }
      catch(std::exception& e) 
      {
        outlog << __func__ << ": Caught standard library exception: " << e.what() << std::endl;
        std::cerr << __func__ << ": Caught standard library exception: " << e.what() << std::endl;
        exit(1);
      }


#ifdef DEBUGGING_RADMAT_DB_INTERFACE
      std::cout << __PRETTY_FUNCTION__ << std::endl;
      std::cout << "KEY ----- " << k << std::endl;
      std::cout << "ENERGY -- " << SEMBLE::toScalar(ENSEM::mean(eval.E())) << std::endl;
      std::cout << "       -- size " << eval.E().size() << std::endl;
      std::cout << "OVERLAP - " << SEMBLE::toScalar(ENSEM::mean(eval.Z())) << std::endl;
      std::cout << "       -- size " << eval.Z().size() << std::endl;
#endif

      return eval;
    } 

} // namespace radmat



#undef DEBUGGING_RADMAT_DB_INTERFACE












#endif /* RADMAT_DATABASE_INTERFACE_H */
