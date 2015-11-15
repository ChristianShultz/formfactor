#ifndef DATA_REPRESENTATION_FACTORY_H
#define DATA_REPRESENTATION_FACTORY_H 

#include "data_representation_primitive_rep.h"
#include "adat/singleton.h"
#include "adat/objfactory.h"
#include "radmat/utils/stringify.h"
#include "radmat/utils/handle.h"
#include <string>


namespace radmat
{


  // factory to hold them all
  typedef Util::SingletonHolder<
    Util::ObjectFactory<Rep_p,
			std::string,
			void,
			Rep_p* (*)(void),
			Util::StringFactoryError> >
  TheDataRepresentationFactory;

  namespace DataRepresentationFactoryEnv
  {
    bool registerAll(void);

    // nb keys are rep_id() -- RepEmbedType do not go into the factory
    rHandle<Rep_p> callFactory(const std::string &id);
  }





}


#endif /* DATA_REPRESENTATION_FACTORY_H */
