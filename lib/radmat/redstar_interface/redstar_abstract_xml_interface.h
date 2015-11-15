#ifndef REDSTAR_ABSTRACT_XML_INTERFACE_H
#define REDSTAR_ABSTRACT_XML_INTERFACE_H 

#include "io/adat_xmlio.h"
#include <string>

// set up an xml abstraction 

namespace radmat
{

  struct AbsRedstarXMLInterface_t
  {
    AbsRedstarXMLInterface_t() {}
    virtual ~AbsRedstarXMLInterface_t(){}

    virtual void read(ADATXML::XMLReader &xml, 
        const std::string &path) = 0;     
    virtual std::string write(void) const = 0; 
    virtual std::string type(void) const = 0; 
  }; 



} // radmat



#endif /* REDSTAR_ABSTRACT_XML_INTERFACE_H */
