#ifndef LLSQ_MULTI_DATA_SERIALIZE_H
#define LLSQ_MULTI_DATA_SERIALIZE_H

#include "llsq_multi_data.h"
#include "ensem/ensem.h"
#include "semble/semble_meta.h"

namespace radmat
{

  // row ordered writing
  template<typename TT, typename ST>
    void write(ADATIO::BinaryWriter &bin, const LLSQMultiData<TT,ST> &d)
    {
      int nr = d.nrows(); 
      write(bin,nr); 
      for (int r = 0; r < nr; ++r)
      {
        write(bin,d.get_tag(r)); 
        ENSEM::write(bin,d.get_row_ensem(r)); 
      }
    }


  // row ordered reading
  template<typename TT, typename ST>
    void read(ADATIO::BinaryReader &bin, LLSQMultiData<TT,ST> &d)
    {
      int nr; 
      read(bin,nr); 
      LLSQMultiData<TT,ST> foo; 
      for(int r = 0; r < nr; ++r)
      {
        TT tag; 
        typename SEMBLE::PromoteEnsemVec<ST>::Type v;
        read(bin,tag); 
        ENSEM::read(bin,v); 
        foo.append_row_ensem(v,tag); 
      }
      d = foo; 
    }


} // radmat 

#endif /* LLSQ_MULTI_DATA_SERIALIZE_H */
