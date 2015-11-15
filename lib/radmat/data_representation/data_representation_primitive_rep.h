#ifndef DATA_REPRESENTATION_PRIMITIVE_REP_H
#define DATA_REPRESENTATION_PRIMITIVE_REP_H 

#include "radmat/utils/stringify.h"
#include <string> 


namespace radmat
{

  // the things we care about live either 
  // in a lorentz rep or a cubic rep 

  struct Rep_p
  {
    virtual ~Rep_p() {}
    virtual std::string rep_type() const = 0; 
    virtual std::string rep_id() const = 0; 
  };  


  template<typename T> 
    struct Rep_t
    : public Rep_p
    {
      virtual ~Rep_t() {}
      virtual std::string rep_type(void) const { return Stringify<T>();}; 
    };


  struct CubicRep_t {}; 
  struct LorentzRep_t {}; 

  REGISTER_STRINGIFY_TYPE(CubicRep_t);
  REGISTER_STRINGIFY_TYPE(LorentzRep_t); 


  struct CubicRep_p 
  : public Rep_t<CubicRep_t> 
  {
    virtual ~CubicRep_p() {}
    virtual std::string rep_id() const = 0; 
  };

  struct LorentzRep_p 
  : public Rep_t<LorentzRep_t> 
  {
    virtual ~LorentzRep_p() {}
    virtual std::string rep_id() const = 0; 
  }; 

}


#endif /* DATA_REPRESENTATION_PRIMITIVE_REP_H */
