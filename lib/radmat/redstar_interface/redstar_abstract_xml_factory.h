#ifndef REDSTAR_ABSTRACT_XML_FACTORY_H
#define REDSTAR_ABSTRACT_XML_FACTORY_H 



#include "adat/objfactory.h"
#include "adat/singleton.h"
#include "radmat/utils/handle.h"
#include "io/adat_xmlio.h"
#include "redstar_abstract_xml_interface.h"
#include <string>

namespace radmat
{

  typedef Util::SingletonHolder<
              Util::ObjectFactory<AbsRedstarXMLInterface_t, 
                                  std::string, 
                                  void,
                                  AbsRedstarXMLInterface_t* (*)(void), 
                                  Util::StringFactoryError> >
      TheRedstarAbstractXMLFactory; 

  namespace TheRedstarAbstractXMLFactoryEnv
  {
    bool registerAll(void); 
    rHandle<AbsRedstarXMLInterface_t> callFactory(const std::string &id); 
  }

} // radmat




#endif /* REDSTAR_ABSTRACT_XML_FACTORY_H */
