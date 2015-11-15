/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : radmat_driver_props.cc

 * Purpose :

 * Creation Date : 29-11-2012

 * Last Modified : Sun 04 May 2014 02:53:47 PM EDT

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/

#include "radmat/driver/radmat_driver_props.h"

#include <sstream>


namespace
{

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




namespace radmat
{

  std::string toString (const RDriverProps_t &prop)
  {
    std::stringstream ss;
    ss << "version = " << prop.version; 
    ss << "\nbigPhase = " << prop.bigPhase; 
    ss << "\nthreePointComparatorProps = " << prop.threePointComparatorProps;
    ss << "\nthreePointIni = " << prop.threePointIni;
    ss << "\nchisq = " << prop.chisq; 
    ss << "\nmaxThread = " << prop.maxThread; 
    ss << "\npoleMass^2 = " <<  prop.poleMass; 
    ss << "\ntolerance = " << prop.tolerance; 
    return ss.str(); 
  }

  std::ostream& operator<<(std::ostream& o, const RDriverProps_t &p)
  {
    o << toString(p);
    return o; 
  }


  namespace
  {
    void check_version(const int v)
    {
      std::vector<int> version_list; 
      std::vector<int>::const_iterator it; 
      bool found = false; 
      version_list.push_back(3); 

      for( it = version_list.begin(); it != version_list.end(); ++it)
        if ( *it == v ) 
          found = true; 

      if ( !!! found ) 
      {
        std::cerr << "version " << v << "is not supported, must be at ";
        for( it = version_list.begin(); it != version_list.end(); ++it)
          std::cerr << *it << " ";
        std::cerr <<"try again later" << std::endl; 
        exit(1); 
      }
    }
  }


  void read(ADATXML::XMLReader &xml, const std::string &path, RDriverProps_t &prop)
  {
    ADATXML::XMLReader ptop(xml,path);
    doXMLRead(ptop,"version",prop.version,__PRETTY_FUNCTION__); 

    check_version(prop.version); 

    // default to no phase
    if(ptop.count("bigPhase") > 0)
      read(ptop,"bigPhase",prop.bigPhase);
    else
    {
      prop.bigPhase = 1; 
    }


    doXMLRead(ptop,"threePointComparatorProps",prop.threePointComparatorProps,__PRETTY_FUNCTION__);
    doXMLRead(ptop,"chisq",prop.chisq,__PRETTY_FUNCTION__);
    doXMLRead(ptop,"threePointIni",prop.threePointIni,__PRETTY_FUNCTION__);
    doXMLRead(ptop,"maxThread",prop.maxThread,__PRETTY_FUNCTION__);
    doXMLRead(ptop,"tolerance",prop.tolerance,__PRETTY_FUNCTION__); 

    // NB: we square the mass here!
    double pole_mass; 
    doXMLRead(ptop,"poleMass",pole_mass,__PRETTY_FUNCTION__); 
    prop.poleMass = pole_mass*pole_mass; 
  } 


} // namespace radmat 
