#ifndef LORENTZFF_CUBIC_REPS_H
#define LORENTZFF_CUBIC_REPS_H 

#include "radmat/utils/stringify.h"
#include "radmat/utils/handle.h"
#include "hadron/irrep_util.h"
#include "io/adat_xmlio.h"
#include "radmat/data_representation/data_representation.h"

namespace radmat
{

  struct RepHandler
  {

    rHandle<Rep_p> gen_rep(const ADATXML::Array<int> &p) const
    {
      std::string LG = Hadron::generateLittleGroup(p); 
      if(LG == "Oh")
        return rHandle<Rep_p>(new Oh() ); 
      else if (LG == "D2")
        return rHandle<Rep_p>(new D2() ); 
      else if (LG == "D3")
        return rHandle<Rep_p>(new D3() ); 
      else if (LG == "D4")
        return rHandle<Rep_p>(new D4() ); 
      else
      {
        std::cout << "LG " << LG << " not supported" << std::endl;
        exit(1); 
      }
      exit(1); 
    }

  };

  struct RepPair
  {
    RepPair(const ADATXML::Array<int> &left, const ADATXML::Array<int> &right)
      : l(left) , r(right)
    { 
      RepHandler R; 
      lefty = R.gen_rep(l);
      righty = R.gen_rep(r); 
    }

    ADATXML::Array<int> l;
    ADATXML::Array<int> r; 
    rHandle<Rep_p> lefty;
    rHandle<Rep_p> righty; 
  }; 

}
#endif /* LORENTZFF_CUBIC_REPS_H */
