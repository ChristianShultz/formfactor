/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : build_correlators_xml.cc

 * Purpose :

 * Creation Date : 25-04-2013

 * Last Modified : Wed 26 Mar 2014 11:15:11 AM EDT

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/


#include "construct_correlators_xml.h"
#include "radmat/utils/printer.h"
#include <exception>


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
          std::cerr << __PRETTY_FUNCTION__ << ": Error, called by " 
            << f << " trying to read path, " << path
            << ", path was empty, exiting" << std::endl;
          exit(1);
        }
      }


  } // namespace anonomyous 

  std::string toString(const ThreePointCorrXMLIni_t::RenormalizationProp &o)
  {
    std::stringstream ss; 
    ss << "RGE_prop = " << o.RGE_prop 
      << " Z_t = " << o.Z_t 
      << " Z_s = " << o.Z_s;

    return ss.str(); 
  }

  std::ostream & operator<<(std::ostream &os, 
      const ThreePointCorrXMLIni_t::RenormalizationProp &p)
  {
    return ( os << toString(p) ); 
  }


  //! xml reader
  void read(ADATXML::XMLReader &xml, 
      const std::string &path, 
      ThreePointCorrXMLIni_t::RenormalizationProp &prop)
  {
    ADATXML::XMLReader ptop(xml,path);
    doXMLRead(ptop,"RGE_prop",prop.RGE_prop,__PRETTY_FUNCTION__);
    doXMLRead(ptop,"Z_t",prop.Z_t,__PRETTY_FUNCTION__);
    doXMLRead(ptop,"Z_s",prop.Z_s,__PRETTY_FUNCTION__);
  }


  //! write it to a string
  std::string toString(const ThreePointCorrXMLIni_t &o)
  { 
    std::stringstream ss;
    ss << "redstar = " <<  o.redstar->type() 
      << " src = " << o.source_id << " snk " << o.sink_id 
      << " maSource = " << o.maSource << " maSink = " << o.maSink
      << "\nRenormalization = " << o.renormalization; 
    return ss.str();
  }

  //! stream it
  std::ostream& operator<<(std::ostream &o, const ThreePointCorrXMLIni_t &p)
  {
    o << toString(p);
    return o;
  }

  //! xml reader
  void read(ADATXML::XMLReader &xml, const std::string &path, ThreePointCorrXMLIni_t &prop)
  {
    ADATXML::XMLReader ptop(xml,path);

    // extra step since we are reading an abstraction and need to 
    // pull the object from the factory first 
    AbstractNamedObject<AbsRedstarXMLInterface_t> foo; 
    if( ptop.count("redstar") > 0 )
      read(ptop,"redstar",foo);
    else
    {
      printer_function<console_print>("missing redstar" + std::string(__PRETTY_FUNCTION__) ); 
    }
    prop.redstar = foo.param; 


    doXMLRead(ptop,"source_id",prop.source_id,__PRETTY_FUNCTION__);
    doXMLRead(ptop,"sink_id",prop.sink_id,__PRETTY_FUNCTION__);
    doXMLRead(ptop,"maSource",prop.maSource,__PRETTY_FUNCTION__);
    doXMLRead(ptop,"maSink",prop.maSink,__PRETTY_FUNCTION__); 
    doXMLRead(ptop,"renormalization",prop.renormalization,__PRETTY_FUNCTION__);
  }

  std::string toString(const ThreePointCorrIni_t &prop)
  {
    std::stringstream ss;
    ss << "threePointCorrXMLIni = " << prop.threePointCorrXMLIni
      << "\nradmatDBProp = "  << prop.radmatDBProp
      << "\nmatElemMode = " << prop.matElemMode
      << "\nmatElemID = " << prop.matElemID
      << " xi = " << prop.xi << " L_s = " << prop.L_s;  
    return ss.str();
  }

  std::ostream& operator<<(std::ostream &o, const ThreePointCorrIni_t &prop)
  {
    o << toString(prop);
    return o;
  }

  void read(ADATXML::XMLReader &xml, const std::string &path, ThreePointCorrIni_t &prop)
  { 
    ADATXML::XMLReader ptop(xml,path);
    doXMLRead(ptop,"threePointCorrXMLIni",prop.threePointCorrXMLIni,__PRETTY_FUNCTION__);
    doXMLRead(ptop,"radmatDBProp",prop.radmatDBProp,__PRETTY_FUNCTION__);
    doXMLRead(ptop,"matElemMode",prop.matElemMode,__PRETTY_FUNCTION__); 
    doXMLRead(ptop,"matElemID",prop.matElemID,__PRETTY_FUNCTION__);
    doXMLRead(ptop,"xi",prop.xi,__PRETTY_FUNCTION__);
    doXMLRead(ptop,"L_s",prop.L_s,__PRETTY_FUNCTION__);
  }

  //  void write(ADATXML::XMLWriter &xml, 
  //      const std::string &path, 
  //      const ThreePointCorrXMLIni_t::RenormalizationProp &prop)
  //  {
  //    ADATXML::push(xml,path);
  //    write(xml,"RGE_prop",prop.RGE_prop);
  //    write(xml,"Z_t",prop.Z_t);
  //    write(xml,"Z_s",prop.Z_s);
  //    ADATXML::pop(xml);
  //  }
  //
  //  //! xml writer
  //  void write(ADATXML::XMLWriter &xml, const std::string &path, const ThreePointCorrXMLIni_t &prop)
  //  {
  //    ADATXML::push(xml,path);
  //    write(xml,"continuumMatElemXML",prop.continuumMatElemXML);
  //    write(xml,"source_id",prop.source_id);
  //    write(xml,"sink_id",prop.sink_id);
  //    write(xml,"maSource",prop.maSource);
  //    write(xml,"maSink",prop.maSink); 
  //    write(xml,"gParitySymmetry",prop.gParitySymmetry); 
  //    write(xml,"cubicSymmetry",prop.cubicSymmetry); 
  //    write(xml,"renormalization",prop.renormalization); 
  //    ADATXML::pop(xml);
  //  }
  //
  //  void write(ADATXML::XMLWriter &xml, const std::string &path, const ThreePointCorrIni_t &prop)
  //  {
  //    ADATXML::push(xml,path);
  //    write(xml,"threePointCorrXMLIni",prop.threePointCorrXMLIni);
  //    write(xml,"radmatDBProp",prop.radmatDBProp);
  //    write(xml,"matElemID",prop.matElemID); 
  //    write(xml,"xi",prop.xi);
  //    write(xml,"L_s",prop.L_s); 
  //    ADATXML::pop(xml);
  //  }

} // namespace radmat


