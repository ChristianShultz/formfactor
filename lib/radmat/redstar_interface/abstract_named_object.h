#ifndef ABSTRACT_NAMED_OBJECT_H
#define ABSTRACT_NAMED_OBJECT_H 

#include <iostream>
#include <string>
#include <algorithm>
#include <exception>
#include "radmat/utils/handle.h"

//   NB: the handles are NOT thread safe so this can 
//       only be used through single threadded sections
//       ie XML reading 


namespace radmat
{
  template<typename T> 
    struct AbstractNamedObject 
    {
      std::string object_name; 
      rHandle<T> param; 

      AbstractNamedObject(void) {}
      AbstractNamedObject(const AbstractNamedObject &o)
      {
        object_name = o.object_name; 
        param = o.param; 
      }

      AbstractNamedObject& operator=(const AbstractNamedObject &o)
      {
        if (this != &o)
        {
          AbstractNamedObject foo(o); 
          
          std::swap(object_name,foo.object_name); 
          std::swap(param,foo.param); 
        }

        return *this; 
      }

    };

}
#endif /* ABSTRACT_NAMED_OBJECT_H */
