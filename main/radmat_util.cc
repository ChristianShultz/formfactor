/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : radmat.cc

 * Purpose :

 * Creation Date : 25-02-2013

 * Last Modified : Fri 17 Oct 2014 11:38:26 AM EDT

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/


#include <string>
#include <sstream>
#include <iostream>
#include <map>
#include <string>
#include <utility>
#include <exception>
#include <iostream>

#include "radmat/ff/lorentzff_canonical_rotations_checker.h"
#include "radmat/register_all/register_all.h"
#include "radmat/driver/radmat_driver.h"
#include "radmat/construct_data/lattice_multi_data_object.h"
#include "radmat/llsq/llsq_multi_data_serialize.h"
#include "radmat/llsq/llsq_solution.h"
#include "radmat/utils/splash.h"
#include "radmat/utils/handle.h"

#include "io/adat_io.h"
#include "io/adat_xmlio.h"

#include "jackFitter/three_point_fit_forms.h"



namespace
{
  template<typename T>
    void doXMLRead(ADATXML::XMLReader &ptop, const std::string &path, T &place, const char * f)
    {
      if(ptop.count(path) > 0)
        read(ptop,path,place);
      else
      {
        std::cerr << __PRETTY_FUNCTION__ << ": Error, called by " 
          << f << " trying to read path, " << path
          << ", path was empty, exiting" << std::endl;
        exit(1);
      }
    }

  template<typename T>
    void doXMLRead(ADATXML::XMLReader &ptop, const std::string &path, T &place, const char * f, const T &val)
    {
      if(ptop.count(path) > 0)
        read(ptop,path,place);
      else
      {
        std::cerr << __PRETTY_FUNCTION__ << ": Error, called by " 
          << f << " trying to read path, " << path
          << ", path was empty, reverting to default value " 
          << val << std::endl;
        place = val;
      }
    }

  struct SingleQ2Prop_t
  {
    std::string ff;                                                 // which form factor are we refitting
    ThreePointComparatorProps_t threePointComparatorProps;  // how are we fitting it
    FitParValue fitParameterValues; 
  };

  struct ArrSingleQ2Prop_t
  {
    ADATXML::Array<SingleQ2Prop_t>  ffs;                    // the list of ffs that we want to refit
    int tsrc;                                               // some duplicate info
    int tsnk;                                               
    std::string ff_dbfile;                                     // where does it live
    std::string s_dbfile; 
    std::string solnID;                                     // how are we inverting   
    double tolerance;                                       // tolerance
    ADATXML::Array<int> lat_elems;                          // which elements are we using in the refit
  };



  void read(ADATXML::XMLReader &xml, const std::string &path, SingleQ2Prop_t &p) 
  {
    ADATXML::XMLReader ptop(xml,path); 

    doXMLRead(ptop,"ff",p.ff,__PRETTY_FUNCTION__); 
    doXMLRead(ptop,"threePointComparatorProps",p.threePointComparatorProps,__PRETTY_FUNCTION__); 
    doXMLRead(ptop,"fitParameterValues",p.fitParameterValues,__PRETTY_FUNCTION__); 
  }


  void read(ADATXML::XMLReader &xml, const std::string &path, ArrSingleQ2Prop_t &p)
  {
    ADATXML::XMLReader ptop(xml,path); 

    doXMLRead(ptop,"ffs",p.ffs,__PRETTY_FUNCTION__); 
    doXMLRead(ptop,"tsrc",p.tsrc,__PRETTY_FUNCTION__); 
    doXMLRead(ptop,"tsnk",p.tsnk,__PRETTY_FUNCTION__); 
    doXMLRead(ptop,"ff_dbfile",p.ff_dbfile,__PRETTY_FUNCTION__); 
    doXMLRead(ptop,"s_dbfile",p.s_dbfile,__PRETTY_FUNCTION__); 
    doXMLRead(ptop,"solnID",p.solnID,__PRETTY_FUNCTION__); 
    doXMLRead(ptop,"tolerance",p.tolerance,__PRETTY_FUNCTION__); 
    doXMLRead(ptop,"lat_elems",p.lat_elems,__PRETTY_FUNCTION__); 
  }


  void prune_llsq_elems(radmat::rHandle< radmat::LLSQLatticeMultiData > &inout , const ArrSingleQ2Prop_t &p)
  {

    radmat::LLSQLatticeMultiData trim;

    std::cout << __func__ << ": pulled " << inout->tags().size() << " correlators from the dbase" << std::endl;

    // parse the list 
    if( p.lat_elems.size() > 0 )
    {
      for (int i = 0; i < p.lat_elems.size(); ++i)
      {
        trim.append_row_semble(inout->get_row_semble(p.lat_elems[i]),inout->get_tag(p.lat_elems[i])); 
      }

      *inout = trim;
    }

    std::cout << __func__ << " pulled :" << std::endl;
    inout->splash_tags();
  }

} // anonomyous 



// UTILITY FUNCTIONS
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////


  
// pull the arr_ini
//
/////////////////////////////////////////////////////
ArrSingleQ2Prop_t 
read_inifile(const std::string &xmlini)
{
// read the xml array of ffs that we want to refit
  ArrSingleQ2Prop_t arr_ini; 

  try
  {
    std::cout << "reading " << xmlini << std::endl;
    ADATXML::XMLReader xml(xmlini); 
    read(xml,"/props",arr_ini); 
  }
  catch( std::string &str) 
  {
    std::cout << "Error: " << str << std::endl; 
    exit(1); 
  }
  catch( std::exception &e)
  {
    std::cout << "excep: " << e.what() << std::endl; 
    exit(1);
  }
  catch ( ... ) 
  {
    SPLASH("An error occurred reading the inifile") ; 
    exit(1); 
  }

  return arr_ini;
}


// pull the llsq database struct 
//
/////////////////////////////////////////////////////
radmat::rHandle<radmat::LLSQLatticeMultiData> 
pull_llsq_data(const std::string &xmlini)
{
  // read the xml array of ffs that we want to refit
  ArrSingleQ2Prop_t arr_ini = read_inifile(xmlini);

  radmat::rHandle< radmat::LLSQLatticeMultiData > foo( new radmat::LLSQLatticeMultiData() ); 
  POW2_ASSERT( &*foo ) ; 

  try
  {
    ADATIO::BinaryFileReader bread(arr_ini.s_dbfile); 
    radmat::read(bread,*foo); 
  }
  catch ( ... ) 
  {
    SPLASH("An error occurred reading the state database");
    exit(1);
  }

  // get the lattice elems as a function of insertion time from 
  // the database that was saved in the orig run 
  prune_llsq_elems(foo,arr_ini); 

  return foo;
}


// pull the formfactor database struct 
//
/////////////////////////////////////////////////////
radmat::FormFacSolutions<std::complex<double> > 
pull_ff_soln(const std::string &xmlini)
{

  // read the xml array of ffs that we want to refit
  ArrSingleQ2Prop_t arr_ini = read_inifile(xmlini);

  radmat::FormFacSolutions<std::complex<double> > FF_of_t;

  try
  {
    ADATIO::BinaryFileReader bread(arr_ini.ff_dbfile); 
    radmat::read(bread,FF_of_t); 
  }
  catch ( ... ) 
  {
    SPLASH("An error occurred reading the state database");
    exit(1);
  }

  return FF_of_t;
}





//    generate redstar xml from the inifile 
//
/////////////////////////////////////////////////////
void gen_xml(int argc , char *argv[] )
{
  if(argc != 4)
  {
    std::cerr << "error: usage: radmat_util: gen_xml <xmlinifile> <mode> " << std::endl;
    exit(1); 
  }

  std::istringstream val(argv[2]); 
  std::string ini; 
  val >> ini; 

  std::istringstream val2(argv[3]); 
  std::string mode; 
  val2 >> mode; 

  radmat::RadmatDriver d; 
  d.xml_handler(ini,mode); 
}


//    just do a registration
//
/////////////////////////////////////////////////////
void do_registerAll(int argc , char *argv[] )
{}

//    use a statedatabase to check the rotation 
//    properties of the data 
//
/////////////////////////////////////////////////////
void rot_llsq(int argc, char *argv[])
{
  if(argc != 3)
  {
    std::cerr << "usage: radmat_util: rot_llsq <xmlinifile> " << std::endl;
    exit(1); 
  }


  // read the name of the ini file 
  std::string xmlini;
  std::istringstream val(argv[2]);
  val >> xmlini;

  radmat::rHandle< radmat::LLSQLatticeMultiData > foo = pull_llsq_data(xmlini);

  radmat::LatticeRotationRelationChecker bar;

  try
  {
    bar.check(foo); 
  }
  catch (std::string &s)
  {
    std::cout << __func__ << ": caught " << s << std::endl;
  }

}




//    use a statedatabase to resolve a llsq 
//
/////////////////////////////////////////////////////
void prune_llsq(int argc, char *argv[])
{
  if(argc != 3)
  {
    std::cerr << "usage: radmat_util:" << __func__ << " <xmlinifile>" << std::endl;
    exit(1); 
  }

  // read the name of the ini file 
  std::string xmlini;
  std::istringstream val(argv[2]);
  val >> xmlini;


  ArrSingleQ2Prop_t arr_ini = read_inifile(xmlini); 
  radmat::rHandle< radmat::LLSQLatticeMultiData > foo = pull_llsq_data(xmlini);
    

  radmat::RadmatSingleQ2Driver m_driver; 
  
  // use the back door
  m_driver.load_llsq(foo,arr_ini.tolerance,false); 

  // solve the llsq 
  m_driver.solve_llsq(arr_ini.solnID); 

  // grab the result 
  radmat::FormFacSolutions<std::complex<double> > solution_thing; 
  solution_thing = m_driver.grab_ff_solution(); 

  // save the state 
  std::stringstream ss; 
  ss << __func__ << ".ff_database.rad"; 
  std::cout << __func__ << ": saving the ff in " << ss.str() << std::endl;
  ADATIO::BinaryFileWriter bin(ss.str()); 
  write(bin,solution_thing); 
  bin.close(); 

  // now do the fits  
  for (int elem = 0; elem < arr_ini.ffs.size(); ++elem)
  {
    std::cout << "\n\n** refiting ff: " << arr_ini.ffs[elem].ff << std::endl;
    std::cout << "*********************************" << std::endl;
    // fit out the insertion time dependence
    m_driver.fit_and_dump_single_ffs(
        arr_ini.ffs[elem].threePointComparatorProps,
        solution_thing,
        arr_ini.tsrc,
        arr_ini.tsnk,
        arr_ini.ffs[elem].ff,
        arr_ini.ffs[elem].fitParameterValues);
  } // ff loop

}


//    go back and refit form factors 
//        -- leave the llsq alone here 
//
/////////////////////////////////////////////////////
void refit_ffs(int argc, char *argv[])
{
  if(argc != 3)
  {
    std::cerr << "usage: radmat_util: refit_ffs <xmlinifile> " << std::endl;
    exit(1); 
  }


  // read the name of the ini file 
  std::string xmlini;
  std::istringstream val(argv[2]);
  val >> xmlini;


  // read the xml array of ffs that we want to refit
  ArrSingleQ2Prop_t arr_ini = read_inifile(xmlini);
  radmat::FormFacSolutions<std::complex<double> > FF_of_t = pull_ff_soln(xmlini); 

  // driver
  radmat::RadmatSingleQ2Driver my_driver;

  
  for (int elem = 0; elem < arr_ini.ffs.size(); ++elem)
  {
    std::cout << "\n\n** refiting ff: " << arr_ini.ffs[elem].ff << std::endl;
    std::cout << "*********************************" << std::endl;
    // fit out the insertion time dependence
    my_driver.fit_and_dump_single_ffs(
        arr_ini.ffs[elem].threePointComparatorProps,
        FF_of_t,
        arr_ini.tsrc,
        arr_ini.tsnk,
        arr_ini.ffs[elem].ff,
        arr_ini.ffs[elem].fitParameterValues);
  } // ff loop

}


//      convert G_i to multipole -- did the fit too
//
/////////////////////////////////////////////////////
void conv_rho_multipole(int argc, char *argv[])
{
  if(argc != 3)
  {
    std::cerr << "usage: radmat_util: conv_rho_multipole <xmlinifile> " << std::endl;
    exit(1); 
  }


  // read the name of the ini file 
  std::string xmlini;
  std::istringstream val(argv[2]);
  val >> xmlini;


  // read the xml array of ffs that we want to refit
  ArrSingleQ2Prop_t arr_ini = read_inifile(xmlini);
  radmat::FormFacSolutions<std::complex<double> > FF_of_t = pull_ff_soln(xmlini); 

  // driver
  radmat::RadmatSingleQ2Driver my_driver;
  
  // map 
  std::map<std::string,ENSEM::EnsemReal> ff_map; 
  std::map<std::string,ENSEM::EnsemReal>::const_iterator ff_it; 
  ENSEM::EnsemReal Q2; 
 
  // run some fake fits  
  for (int elem = 0; elem < arr_ini.ffs.size(); ++elem)
  {
    std::cout << "\n\n** refiting ff: " << arr_ini.ffs[elem].ff << std::endl;
    std::cout << "*********************************" << std::endl;

    std::pair<ENSEM::EnsemReal,ENSEM::EnsemReal> bob; 

    // fit out the insertion time dependence
    bob = my_driver.fit_and_dump_single_ffs(
        arr_ini.ffs[elem].threePointComparatorProps,
        FF_of_t,
        arr_ini.tsrc,
        arr_ini.tsnk,
        arr_ini.ffs[elem].ff,
        arr_ini.ffs[elem].fitParameterValues);

    ff_map.insert(std::make_pair(arr_ini.ffs[elem].ff, bob.first)); 
    Q2 = bob.second; 
  } // ff loop

  ENSEM::EnsemReal GC,GQ,GM; 
  ENSEM::EnsemReal G1,G2,G3; 
  ENSEM::Real mrho(0.216292), one(1.), o6(1/6.), o4(1./4.);
  ENSEM::EnsemReal Qdm = Q2/mrho/mrho; 

  bool bG1(true),bG2(true),bG3(true); 

  ff_it = ff_map.find("RhoRhoG1");
  if(ff_it == ff_map.end())
    bG1 = false;
  else 
    G1 = ff_it->second; 

  ff_it = ff_map.find("RhoRhoG2");
  if(ff_it == ff_map.end())
    bG2 = false;
  else 
    G2 = ff_it->second; 

  ff_it = ff_map.find("RhoRhoG3");
  if(ff_it == ff_map.end())
    bG3 = false;
  else 
    G3 = ff_it->second; 


  if( bG1 && bG2 && bG3 )
  {
    GC = (one + Qdm*o6)*G1 
        - Qdm*o6*G2
        + Qdm*o6*(one + Qdm*o4)*G3; 
    std::string s("RhoRhoGC_fit.jack");
    ENSEM::write(s,GC); 
  } 

  if(  bG2  )
  {
    GM = G2;
    std::string s("RhoRhoGM_fit.jack");
    ENSEM::write(s,GM); 
  } 

  if( bG1 && bG2 && bG3 )
  {
    GQ = G1 
        - G2
        + (one + Qdm*o4)*G3; 
    std::string s("RhoRhoGQ_fit.jack");
    ENSEM::write(s,GQ); 
  } 
  
}


//
//
//    WORK HANDLER STUFF
//
//

// typedef 
typedef void (*fptr)(int argc , char *argv[]) ; 

// a map of operation names and function pointers
std::map<std::string , fptr> options; 

// init the map 
void init_options(void)
{
  options.insert(std::pair<std::string,fptr>("gen_xml",&gen_xml)); 
  options.insert(std::pair<std::string,fptr>("rot_llsq",&rot_llsq)); 
  options.insert(std::pair<std::string,fptr>("prune_llsq",&prune_llsq)); 
  options.insert(std::pair<std::string,fptr>("refit_ffs",&refit_ffs)); 
  options.insert(std::pair<std::string,fptr>("registerAll",&do_registerAll));
  options.insert(std::pair<std::string,fptr>("conv_rho_multipole",&conv_rho_multipole));
}

// pick appropriate function and pass on command line inputs 
void do_work(std::string &op, int argc,char *argv[])
{
  init_options(); 

  if(options.find(op) == options.end())
  {
    std::cerr << " unrecognized op " << op 
      << " more intelligent choices are " << std::endl; 
    std::map<std::string , fptr>::const_iterator it; 
    for(it = options.begin(); it != options.end(); ++it)
      std::cerr << it->first << std::endl; 
    exit(1); 
  }

  fptr foo = options[op];

  foo(argc,argv); 
}


// main program wrapper
int main(int argc, char *argv[])
{
  radmat::AllFactoryEnv::registerAll(); 

  // we will always have at least 2 , radmat_util operation_with_no_inputs
  if(argc < 2)
  {
    std::cerr << "usage: radmat_util : <operation> <op inputs ...> " << std::endl;
    exit(1); 
  }

  std::string op;
  std::istringstream opi(argv[1]); 
  opi >> op; 

  do_work(op,argc,argv); 

  return 0;
}
