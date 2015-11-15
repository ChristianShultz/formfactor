#ifndef DATA_REPRESENTATION_THREE_POINT_H
#define DATA_REPRESENTATION_THREE_POINT_H 

#include "data_representation_prim.h"
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

  struct DataRep3pt; 
  REGISTER_STRINGIFY_TYPE(DataRep3pt); 


  // where does a matrix element live 
  struct DataRep3pt
    : public DataRepPrim
  {
    DataRep3pt() {}

    DataRep3pt( const std::string &ll,
        const std::string &gg, 
        const std::string &rr)
      : l(ll) , r(rr) , g(gg) 
    { }

    virtual ~DataRep3pt() {}

    virtual std::string type() const { return Stringify<DataRep3pt>(); }

    virtual rHandle<Rep_p> lefty() const { return this->call_rep_prim(l); }
    virtual rHandle<Rep_p> righty() const { return this->call_rep_prim(r); }
    virtual rHandle<Rep_p> gamma() const { return this->call_rep_prim(g); }

    std::string l,r,g; 
  };


  // binary read write for serialization routines
  void write(ADATIO::BinaryWriter &, const DataRep3pt &); 

  void read(ADATIO::BinaryReader &, DataRep3pt &); 


}


#endif /* DATA_REPRESENTATION_THREE_POINT_H */
