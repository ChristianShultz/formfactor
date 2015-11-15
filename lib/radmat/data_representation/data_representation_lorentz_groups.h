#ifndef DATA_REPRESENTATION_LORENTZ_GROUPS_H
#define DATA_REPRESENTATION_LORENTZ_GROUPS_H 

#include "data_representation_primitive_rep.h"
#include <string>
#include <vector>
#include "radmat/utils/handle.h"


namespace radmat
{
  template<typename T> 
    struct LorentzRepEmbedType
    : public LorentzRep_p
    { 
      virtual ~LorentzRepEmbedType() {}
      virtual std::string rep_id() const { return Stringify<T>(); }
    };

  struct J0p;  
  struct J0m;
  struct J1p;  
  struct J1m;
  struct J2p;  
  struct J2m;
  struct J3p;  
  struct J3m;
  struct J4p;  
  struct J4m;  

  REGISTER_STRINGIFY_TYPE( J0p ); 
  REGISTER_STRINGIFY_TYPE( J0m ); 
  REGISTER_STRINGIFY_TYPE( J1p ); 
  REGISTER_STRINGIFY_TYPE( J1m ); 
  REGISTER_STRINGIFY_TYPE( J2p ); 
  REGISTER_STRINGIFY_TYPE( J2m ); 
  REGISTER_STRINGIFY_TYPE( J3p ); 
  REGISTER_STRINGIFY_TYPE( J3m ); 
  REGISTER_STRINGIFY_TYPE( J4p ); 
  REGISTER_STRINGIFY_TYPE( J4m ); 

  struct J0p : public LorentzRepEmbedType<J0p> { enum{ SPIN = 0}; };  
  struct J0m : public LorentzRepEmbedType<J0m> { enum{ SPIN = 0}; };
  struct J1p : public LorentzRepEmbedType<J1p> { enum{ SPIN = 1}; };  
  struct J1m : public LorentzRepEmbedType<J1m> { enum{ SPIN = 1}; };
  struct J2p : public LorentzRepEmbedType<J2p> { enum{ SPIN = 2}; };  
  struct J2m : public LorentzRepEmbedType<J2m> { enum{ SPIN = 2}; };
  struct J3p : public LorentzRepEmbedType<J3p> { enum{ SPIN = 3}; };  
  struct J3m : public LorentzRepEmbedType<J3m> { enum{ SPIN = 3}; };
  struct J4p : public LorentzRepEmbedType<J4p> { enum{ SPIN = 4}; };  
  struct J4m : public LorentzRepEmbedType<J4m> { enum{ SPIN = 4}; };  


  // base class 
  struct LorentzRep
    : public LorentzRep_p
  {
    virtual ~LorentzRep() {}
    virtual std::string rep_id() const = 0; 
    virtual int rep_parity() const = 0; 
    virtual int rep_spin() const = 0;
    virtual int rep_eta_p() const = 0; 
    virtual LorentzRep * clone() const = 0; 
  }; 

  // template instantion of parameters 
  template<typename T, int parity> 
    struct LorentzRep_type
    : public LorentzRep
    {
      virtual ~LorentzRep_type() {}
      virtual std::string rep_id() const { return Stringify<T>(); }
      virtual int rep_parity() const { return parity; }
      virtual int rep_spin() const { return T::SPIN; }
      virtual int rep_eta_p() const 
      { return (rep_spin() % 2 == 0) ? rep_parity() : -rep_parity();}
      virtual LorentzRep* clone() const { return new LorentzRep_type(*this); }
    };


  struct J0pRep_t : public LorentzRep_type<J0p,1>  { virtual ~J0pRep_t() {} };
  struct J0mRep_t : public LorentzRep_type<J0m,-1> { virtual ~J0mRep_t() {} };

  struct J1pRep_t : public LorentzRep_type<J1p,1>  { virtual ~J1pRep_t() {} };
  struct J1mRep_t : public LorentzRep_type<J1m,-1> { virtual ~J1mRep_t() {} };

  struct J2pRep_t : public LorentzRep_type<J2p,1>  { virtual ~J2pRep_t() {} };
  struct J2mRep_t : public LorentzRep_type<J2m,-1> { virtual ~J2mRep_t() {} };

  struct J3pRep_t : public LorentzRep_type<J3p,1>  { virtual ~J3pRep_t() {} };
  struct J3mRep_t : public LorentzRep_type<J3m,-1> { virtual ~J3mRep_t() {} };

  struct J4pRep_t : public LorentzRep_type<J4p,1>  { virtual ~J4pRep_t() {} };
  struct J4mRep_t : public LorentzRep_type<J4m,-1> { virtual ~J4mRep_t() {} };


  namespace LorentzRepresentationFactoryEnv
  {
    bool registerAll(); 
    std::vector<std::string> spher_keys();
    rHandle<LorentzRep> callFactory(const std::string &id); 
  }


} 



#endif /* DATA_REPRESENTATION_LORENTZ_GROUPS_H */
