#ifndef REDSTAR_IMPROVED_VECTOR_CURRENT_XML_INTERFACE_H
#define REDSTAR_IMPROVED_VECTOR_CURRENT_XML_INTERFACE_H 

#include "redstar_abstract_xml_interface.h"
#include "redstar_vector_current_xml_interface.h"
#include "radmat/utils/stringify.h"


  //
  //
  // A BASIC gamma^mu VECTOR CURRENT XML INTERFACE
  //
  //

namespace radmat
{

  struct RedstarImprovedVectorCurrentXML;
  REGISTER_STRINGIFY_TYPE(RedstarImprovedVectorCurrentXML);


  struct RedstarImprovedVectorCurrentXML 
    : public AbsRedstarXMLInterface_t 
  {
    typedef RedstarVectorCurrentXML::pfrag pfrag; 
    typedef RedstarVectorCurrentXML::insertion insertion; 

    virtual void read(ADATXML::XMLReader &xml, 
        const std::string &path); 

    virtual std::string type(void) const
    {
      return std::string(Stringify<RedstarImprovedVectorCurrentXML>());
    }

    virtual std::string write(void) const; 


    struct improvement
    {
     double Nu_s;      // tuning parameter 
     double Xi_0;      // bare anisotropy  
     std::string name; 
     double coeff_r; 
     double coeff_i; 
    };

    int pmin; 
    int pmax; 
    int t_slice; 
    insertion time; 
    insertion space; 
    improvement imp; 
  };


  void read(ADATXML::XMLReader &xml, 
      const std::string &path, 
      RedstarImprovedVectorCurrentXML::improvement &); 

  void write(ADATXML::XMLWriter &xml, 
      const std::string &path, 
      RedstarImprovedVectorCurrentXML::improvement &); 

}

#endif /* REDSTAR_IMPROVED_VECTOR_CURRENT_XML_INTERFACE_H */
