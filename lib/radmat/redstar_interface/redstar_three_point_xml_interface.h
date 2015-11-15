#ifndef REDSTAR_THREE_POINT_XML_INTERFACE_H
#define REDSTAR_THREE_POINT_XML_INTERFACE_H 

#include "redstar_abstract_xml_interface.h"
#include "abstract_named_object.h"
#include "radmat/utils/stringify.h"


namespace radmat
{

  struct RedstarThreePointXML; 
  REGISTER_STRINGIFY_TYPE(RedstarThreePointXML); 

  struct RedstarThreePointXML
    : public AbsRedstarXMLInterface_t
  {
    virtual ~RedstarThreePointXML() {}

    virtual std::string type() const
    {
      return Stringify<RedstarThreePointXML>(); 
    }

    virtual void read(ADATXML::XMLReader &xml, 
        const std::string &path); 

    virtual std::string write() const; 

    virtual ADATXML::Array<AbstractNamedObject<AbsRedstarXMLInterface_t> >
      get_npt() const { return npt; }

    virtual std::string get_ensemble() const { return ensemble; }

    virtual const ADATXML::Array<int> & timeslice_info() const { return timeslice; }
    virtual void set_timeslice(const ADATXML::Array<int> &t) {timeslice = t;} 

    ADATXML::Array<AbstractNamedObject<AbsRedstarXMLInterface_t> > npt; 
    ADATXML::Array<int> timeslice; 
    std::string ensemble; 
  }; 

  void read(ADATXML::XMLReader &xml, 
      const std::string &path, 
      AbstractNamedObject<AbsRedstarXMLInterface_t> &); 

  namespace TheRedstarThreePointAbstractXMLFactoryEnv
  {
    rHandle<RedstarThreePointXML> callFactory(const std::string &id); 
  }

} // radmat 

#endif /* REDSTAR_THREE_POINT_XML_INTERFACE_H */
