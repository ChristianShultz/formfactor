/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : single_plot.cc

 * Purpose :

 * Creation Date : 10-01-2013

 * Last Modified : Wed 30 Apr 2014 02:05:04 PM EDT

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/


#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <exception>
#include <vector>
#include <utility>
#include <map>
#include <limits>
#include "io/adat_io.h"
#include "io/adat_xmlio.h"
#include "io/key_val_db.h"
#include "AllConfStoreDB.h"
#include "ensem/ensem.h"
#include "hadron/ensem_filenames.h"
#include "radmat/database/database.h"
#include "hadron/hadron_npart_npt_corr.h"
#include "semble/semble_meta.h"
#include "semble/semble_vector.h"


using namespace radmat;
using namespace Hadron;
using namespace ENSEM;
using namespace ADATXML;
using namespace std;
using namespace FILEDB;

typedef   Hadron::KeyHadronNPartNPtCorr_t KEY; 
typedef   ENSEM::EnsemVectorComplex DATA; 


struct XML_input_t
{
  std::string dbname;
  std::string gnu_name;
  std::string mode; // all will do all corrs in the db, single selects the one in the xml file
  Hadron::KeyHadronNPartNPtCorr_t corr; 
};


namespace
{
  int usage (int argc, char** argv)
  {
    cerr << "Usage: " << argv[0] << " <XML>" << endl;
    return -1;
  }


  template<typename T>
    void doXMLRead(ADATXML::XMLReader &ptop, const std::string &path, T &place, const char * f)
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
} // namespace anonomyous 


void read(ADATXML::XMLReader &xml, const std::string &path, XML_input_t &prop)
{
  ADATXML::XMLReader ptop(xml,path);
  doXMLRead(ptop,"dbname",prop.dbname,__PRETTY_FUNCTION__);
  doXMLRead(ptop,"gnu_name",prop.gnu_name,__PRETTY_FUNCTION__);
  doXMLRead(ptop,"mode",prop.mode,__PRETTY_FUNCTION__); 
  doXMLRead(ptop,"corr",prop.corr,__PRETTY_FUNCTION__);
}




//! Get a key/value
  template<typename K, typename V>
V printKeyValue(const K& ky, 
    AllConfStoreDB< SerialDBKey<K>,  SerialDBData<typename EnsemScalar<V>::Type_t> >& database)
{
  typedef typename EnsemScalar<V>::Type_t SV;

  SerialDBKey<K> key;
  key.key() = ky;

  std::vector< SerialDBData<SV> > vals;
  int ret;
  if ((ret = database.get(key, vals)) != 0)
  {
    std::cerr << __func__ << ": key not found\n" << ky;
    exit(1);
  }

  V eval;
  eval.resize(vals.size());
  eval.resizeObs(vals[0].data().numElem());

  for(int i=0; i < vals.size(); ++i)
  {
    SV sval = vals[i].data();
    pokeEnsem(eval, sval, i);
  }

  return eval;
}


void doPlot(const DATA &m_d, const XML_input_t &ini , const KEY &k)
{

  int nele = m_d.numElem();

  ENSEM::EnsemVectorReal real= ENSEM::real(m_d);
  ENSEM::EnsemVectorReal imag = ENSEM::imag(m_d); 

  ENSEM::VectorReal rm,rv;
  ENSEM::VectorReal im,iv;

  rm = ENSEM::mean(real);
  rv = ENSEM::variance(real);
  im = ENSEM::mean(imag); 
  iv = ENSEM::variance(imag); 


  std::stringstream ssi, ssr; 


  for(int i = 0; i < nele; ++i)
  {
    ssi << i << " " << ENSEM::toDouble(im.elem().elem(i)) << " " << sqrt(ENSEM::toDouble(iv.elem().elem(i))) << "\n";
    ssr << i << " " << ENSEM::toDouble(rm.elem().elem(i)) << " " << sqrt(ENSEM::toDouble(rv.elem().elem(i))) << "\n";
  }


  std::string rout ,iout;
  rout = ini.gnu_name + std::string("_real.dat"); 
  iout = ini.gnu_name + std::string("_imag.dat");


  std::ofstream out;
  out.open(rout.c_str());
  out << ssr.str();
  out.close();

  out.open(iout.c_str());
  out << ssi.str();
  out.close();


  std::string plot_file = ini.gnu_name + std::string("_plot.gp"); 
  out.open(plot_file.c_str());

  std::string title = Hadron::ensemFileName(k);


  out << "set title \""; 

  for(unsigned int c = 0; c < title.length(); ++c)
  {
    out << title[c];
    if(c != 0)
      if(c % 80 == 0)
        out << "\\n";
  }



  out << "\" \n";
  out << "plot '" << rout << "' using 1:2:3 with yerr ls 1 title 'Real' , \\\n"
    << "'" << iout << "' using 1:2:3 with yerr ls 3 title 'Imag' \n"; 

  out.close();

  std::string corrf = ini.gnu_name + std::string("_corr.jack"); 
  ENSEM::write(corrf,m_d);
  ENSEM::write(ini.gnu_name + std::string("_corr_real.jack") , real); 
  ENSEM::write(ini.gnu_name + std::string("_corr_imag.jack") , imag); 


  std::stringstream do_plot; 
  do_plot << "gnuplot -persist " << plot_file; 
  system( do_plot.str().c_str() ); 



}


void my_pause(void)
{
  std::cout << "Press ENTER to continue...";
  std::cin.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );
}



int main(int argc , char *argv[])
{

  if(argc < 2)
    return usage(argc,argv);


  XML_input_t ini;
  std::string xmlinifile;
  std::stringstream val(argv[1]);
  val >> xmlinifile;

  try
  {
    ADATXML::XMLReader xml(xmlinifile);
    read(xml,"/Params",ini);
  }
  catch(const std::string &e)
  {
    std::cerr << __func__ << ": ERROR: can't read xmlinifile (" << xmlinifile << "): " << e << std::endl;
    exit(1);
  }


  typedef EnsemScalar<DATA>::Type_t SV;

  // Open DB
  AllConfStoreDB< SerialDBKey<KEY>,  SerialDBData<SV> > database;

  if (database.open(ini.dbname, O_RDONLY, 0400) != 0)
  {
    std::cerr << __func__ << ": error opening dbase= " << ini.dbname << std::endl;
    exit(1);
  }

  if(ini.mode == "single")
  {
    DATA m_d = printKeyValue<KEY,DATA>(ini.corr,database); 
    doPlot(m_d,ini,ini.corr); 
  }
  else if (ini.mode == "all")
  {
    // Open DB
    AllConfStoreDB< SerialDBKey<KEY>,  SerialDBData<SV> > database;
    if (database.open(ini.dbname, O_RDONLY, 0400) != 0)
    {
      std::cerr << __func__ << ": error opening dbase= " << ini.dbname << std::endl;
      exit(1);
    }

    std::vector< SerialDBKey<KEY> > keys;
    database.keys(keys);

    std::vector<SerialDBKey<KEY> >::const_iterator it;
    for(it = keys.begin(); it != keys.end(); ++it)
    {
      doPlot(printKeyValue<KEY,DATA>(it->key(),database),ini, it->key());
      my_pause(); 
    }

  }
  else
  {
    std::cout << "Error: unrecognized mode(" << ini.mode << ") " << std::endl; 
    exit(1); 
  }

  return 0; 

}
