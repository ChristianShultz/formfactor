#ifndef DATA_TAG_PRIMITIVE_H
#define DATA_TAG_PRIMITIVE_H 


#include <string>


// the base class tag type 

namespace radmat
{
  struct DataTagPrimitive
  {
    virtual ~DataTagPrimitive() {}
    virtual std::string type() const = 0; 
  };
}


#endif /* DATA_TAG_PRIMITIVE_H */
