#ifndef DATA_REPRESENTATION_PRIM_H
#define DATA_REPRESENTATION_PRIM_H 


#include "data_representation_primitive_rep.h"
#include "data_representation_factory.h"
#include "data_representation_cubic_groups.h"
#include "data_representation_lorentz_groups.h"
#include "radmat/utils/handle.h"
#include "radmat/utils/stringify.h"
#include "io/adat_io.h"
#include "io/adat_xmlio.h"
#include <string>


namespace radmat
{

  // where does a matrix element live 
  //      primitive type 
  struct DataRepPrim
  {
    virtual ~DataRepPrim() {}

    virtual std::string type() const = 0; 

    virtual rHandle<Rep_p> call_rep_prim(const std::string &s) const 
    {
      return DataRepresentationFactoryEnv::callFactory(s); 
    }

    virtual rHandle<LorentzRep> call_rep_lor(const std::string &s) const
    {
      return LorentzRepresentationFactoryEnv::callFactory(s); 
    }

    virtual rHandle<CubicRep> call_rep_cub(const std::string &s) const 
    {
      return CubicRepresentationFactoryEnv::callFactory(s); 
    }

    bool is_cubic(const std::string &s) const
    {
      rHandle<Rep_p> r = call_rep_prim(s); 
      return r->rep_type() == Stringify<CubicRep_t>(); 
    }

    bool is_lorentz(const std::string &s) const
    {
      rHandle<Rep_p> r = call_rep_prim(s); 
      return r->rep_type() == Stringify<LorentzRep_t>(); 
    }
  };

}



#endif /* DATA_REPRESENTATION_PRIM_H */
