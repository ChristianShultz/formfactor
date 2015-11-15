/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : corr_norm.cc

 * Purpose :

 * Creation Date : 07-05-2013

 * Last Modified : Tue 10 Dec 2013 11:43:49 AM EST

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/
#include "radmat/register_all/register_all.h"
#include "radmat/construct_data/radmat_database_interface.h"
#include "radmat/fake_data/fake_3pt_function_aux.h"
#include "radmat/utils/pow2assert.h"
#include "radmat/utils/perThreadStorage.h"
#include "radmat/llsq/llsq_gen_system.h"
#include "semble/semble_semble.h"
#include "radmat/construct_data/radmat_overlap_key_val_db.h"
#include "hadron/ensem_filenames.h"
#include "hadron/hadron_npart_npt_corr.h"
#include "io/adat_xmlio.h"
#include "ensem/ensem.h"
#include "formfac/formfac_qsq.h"
#include "jackFitter/plot.h"
#include "ensem/ensem.h"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <exception>


using namespace radmat; 


  template<typename T>
void do_XMLRead(ADATXML::XMLReader &ptop, const std::string &path, T &place, const char * f)
{
  if(ptop.count(path) > 0)
    read(ptop,path,place);
  else
  {
    std::cerr << __PRETTY_FUNCTION__ << ": Error, called by " << f << " trying to read path, " << path
      << ", path was empty, exiting" << std::endl;
    exit(1);
  }
}


struct input_prop
{
  radmatDBProp_t dbProps;
  Hadron::KeyHadronNPartNPtCorr_t npoint; 
  std::string source;
  std::string sink; 
};

void read(ADATXML::XMLReader &xml, const std::string &path, input_prop &prop)
{
  ADATXML::XMLReader ptop(xml,path); 
  do_XMLRead(ptop,"dbProps",prop.dbProps,__PRETTY_FUNCTION__); 
  do_XMLRead(ptop,"npoint",prop.npoint,__PRETTY_FUNCTION__); 
  do_XMLRead(ptop,"source",prop.source,__PRETTY_FUNCTION__); 
  do_XMLRead(ptop,"sink",prop.sink,__PRETTY_FUNCTION__); 
}


typedef radmatAllConfDatabaseInterface< Hadron::KeyHadronNPartNPtCorr_t,
        ENSEM::EnsemVectorComplex,
        RadmatExtendedKeyHadronNPartIrrep_t,
        RadmatMassOverlapData_t> DatabaseInterface_t;


// NB hardwire in here
struct dbKeyObject
{    
  typedef  Hadron::KeyHadronNPartNPtCorr_t nptkey;
  typedef RadmatExtendedKeyHadronNPartIrrep_t normkey;

  dbKeyObject(const nptkey &k , const std::string &id_sink, const std::string &id_source)
    : npt(k)
  {
    sink = normkey(id_sink,k.npoint[1].irrep);
    source = normkey(id_source,k.npoint[3].irrep); 
  }

  nptkey npt;
  normkey sink,source;  
};  



void do_work(const dbKeyObject &d, const DatabaseInterface_t &db)
{
  EnsemVectorComplex corr_tmp = db.fetch(d.npt);
  RadmatMassOverlapData_t source = db.fetch(d.source); 
  RadmatMassOverlapData_t sink = db.fetch(d.sink); 


  std::string outstem = Hadron::ensemFileName(d.npt); 
  std::string pth = SEMBLE::SEMBLEIO::getPath();
  std::stringstream path;
  path << pth << "some_data";
  ENSEM::write(path.str() + std::string("_corr_pre") , corr_tmp); 
  ENSEM::EnsemVectorComplex norm = corr_tmp * SEMBLE::toScalar(0.);


  int t_source = d.npt.npoint[3].t_slice;
  int t_sink = d.npt.npoint[1].t_slice; 

  ThreePtPropagationFactor<double> propagation_factor;

  // NB: the indexing here assumes [tsource,tsink] ie: inclusive range
  for(int t_ins = t_source; t_ins < t_sink; ++t_ins)
  {

    ENSEM::EnsemReal prop = propagation_factor(sink.E(),sink.Z(),t_sink,t_ins,
        source.E(),source.Z(),t_source);

    ENSEM::pokeObs(corr_tmp,ENSEM::peekObs(corr_tmp,t_ins)/prop,t_ins);

    ENSEM::pokeObs(norm,prop,t_ins); 
  } // end loop over t_ins

  ENSEM::write(path.str() + std::string("_corr_post"), corr_tmp);
  ENSEM::write(path.str() + std::string("_norm") , norm);

}




int main(int argc , char* argv[])
{

  AllFactoryEnv::registerAll(); 

  if(argc != 2)
  {
    std::cerr << "error: usage corr_norm <props.xml> " << std::endl;
    exit(1); 
  }

  input_prop ini; 
  std::string xmlini;

  try
  {
    std::istringstream val(argv[1]);
    val >> xmlini;
    ADATXML::XMLReader xml(xmlini); 
    read(xml,"Props",ini);
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
    SPLASH("An error occured while loading the fake data inifile");
    exit(1);
  }



  DatabaseInterface_t db(ini.dbProps);

  dbKeyObject foo(ini.npoint,ini.sink,ini.source);


  do_work(foo,db); 

  return 0; 
}
