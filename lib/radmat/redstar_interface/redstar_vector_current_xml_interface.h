#ifndef REDSTAR_VECTOR_CURRENT_XML_INTERFACE_H
#define REDSTAR_VECTOR_CURRENT_XML_INTERFACE_H 

#include "redstar_abstract_xml_interface.h"
#include "radmat/utils/stringify.h"


  //
  //
  // A BASIC gamma^mu VECTOR CURRENT XML INTERFACE
  //
  //

namespace radmat
{

  struct RedstarVectorCurrentXML;
  REGISTER_STRINGIFY_TYPE(RedstarVectorCurrentXML);


  struct RedstarVectorCurrentXML 
    : public AbsRedstarXMLInterface_t 
  {

    virtual void read(ADATXML::XMLReader &xml, 
        const std::string &path); 

    virtual std::string type(void) const
    {
      return std::string(Stringify<RedstarVectorCurrentXML>());
    }

    virtual std::string write(void) const; 

    // a photon fragment 
    //      -- this is a weight times the op_name 
    struct pfrag
    {
      double coeff_r; 
      double coeff_i; 
      std::string name; 
    };

    struct insertion
    {
      bool active; 
      bool smearedP; 
      ADATXML::Array<pfrag> photons;  
    };

    int pmin; 
    int pmax; 
    int t_slice; 
    insertion time; 
    insertion space; 
  };


  void read(ADATXML::XMLReader &xml, 
      const std::string &path, 
      RedstarVectorCurrentXML::pfrag &); 

  void read(ADATXML::XMLReader &xml, 
      const std::string &path, 
      RedstarVectorCurrentXML::insertion &); 

  void write(ADATXML::XMLWriter &xml, 
      const std::string &path, 
      RedstarVectorCurrentXML::pfrag &); 

  void write(ADATXML::XMLWriter &xml, 
      const std::string &path, 
      RedstarVectorCurrentXML::insertion &); 


}

#endif /* REDSTAR_VECTOR_CURRENT_XML_INTERFACE_H */
