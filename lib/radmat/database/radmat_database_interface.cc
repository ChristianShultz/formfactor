/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : radmat_database_interface.cc

 * Purpose :

 * Creation Date : 04-12-2012

 * Last Modified : Fri 15 Nov 2013 04:50:23 PM EST

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/


#include "radmat_database_interface.h"
#include "io/adat_xmlio.h"
#include <sstream>


namespace radmat
{

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

    template<typename T>
      std::ostream& operator<<(std::ostream &o, const std::vector<T> &tt)
      {
        typename std::vector<T>::const_iterator it; 
        for (it = tt.begin(); it != tt.end(); ++it)
          o << *it << " "; 
        return o; 
      }
  } // namespace anonomyous 


  std::string toString(const dbProp_t &prop)
  {
    std::stringstream ss;
    ss << "dbname = " << prop.dbname << " badlist = " << prop.badlist;
    return ss.str();
  } 

  std::ostream& operator<<(std::ostream &o, const dbProp_t &prop)
  {
    o << toString(prop);
    return o;
  }

  void read(ADATXML::XMLReader &xml, const std::string &path, dbProp_t &prop)
  {
    ADATXML::XMLReader ptop(xml,path);
    doXMLRead(ptop,"dbname",prop.dbname,__PRETTY_FUNCTION__);
    doXMLRead(ptop,"badlist",prop.badlist,__PRETTY_FUNCTION__);
    if(prop.dbname.size() == 0)
      std::cerr <<__PRETTY_FUNCTION__ << ": warning, dbname had no elems" << std::endl; 
  }

  void write(ADATXML::XMLWriter &xml, const std::string &path, const dbProp_t &prop)
  {
    ADATXML::push(xml,path);
    write(xml,"dbname",prop.dbname);
    write(xml,"badlist",prop.badlist);
    ADATXML::pop(xml);
  }


  std::string toString(const radmatDBProp_t &prop)
  {
    std::stringstream ss;
    ss << "threePointDatabase = " << prop.threePointDatabase << "\nnormalizationDatabase = " << prop.normalizationDatabase;
    return ss.str();
  }

  std::ostream& operator<<(std::ostream &o, const radmatDBProp_t &prop)
  {
    o << toString(prop);
    return o;
  }

  void read(ADATXML::XMLReader &xml, const std::string &path, radmatDBProp_t &prop)
  {
    ADATXML::XMLReader ptop(xml,path);
    doXMLRead(ptop,"threePointDatabase",prop.threePointDatabase,__PRETTY_FUNCTION__);
    doXMLRead(ptop,"normalizationDatabase",prop.normalizationDatabase,__PRETTY_FUNCTION__);
    doXMLRead(ptop,"allow_daggering",prop.allow_daggering,__PRETTY_FUNCTION__);
    doXMLRead(ptop,"LG_symmetry",prop.LG_symmetry,__PRETTY_FUNCTION__);  
  }

  void write(ADATXML::XMLWriter &xml, const std::string &path, const radmatDBProp_t &prop)
  {
    ADATXML::push(xml,path);
    write(xml,"threePointDatabase",prop.threePointDatabase);
    write(xml,"normalizationDatabase",prop.normalizationDatabase);
    write(xml,"allow_daggering",prop.allow_daggering); 
    write(xml,"LG_symmetry",prop.LG_symmetry);
    ADATXML::pop(xml);
  }



} // namespace radmat

