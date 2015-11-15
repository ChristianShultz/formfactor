#ifndef WIGNER_D_MATRIX_FACTORY_H
#define WIGNER_D_MATRIX_FACTORY_H 


#include "radmat/utils/tensor.h"
#include "radmat/utils/handle.h"
#include "adat/singleton.h"
#include "adat/objfactory.h"
#include "io/adat_xmlio.h"
#include "radmat/redstar_interface/redstar_canonical_rotations.h"
#include <complex>
#include <map>
#include <string>

namespace radmat
{
  typedef Tensor<std::complex<double>,2> WignerMatrix_t; 

  struct WignerKey
  {
    WignerKey() {}
    WignerKey( const mom_t &p, const int JJ )
      :px(p[0]) , py(p[1]) , pz(p[2]) , J(JJ)
    {}
    int px,py,pz,J; 
  };

  struct WignerKeyClassComp
  {
    bool operator()(const WignerKey &l, const WignerKey &r)
    {
      if(l.px != r.px)
        return l.px < r.px;
      if(l.py != r.py) 
        return l.py < r.py; 
      if(l.pz != r.pz)
        return l.pz < r.pz; 
      return l.J < r.J; 
    }
  };

  namespace WignerDMatrixEnv
  {
    typedef Util::SingletonHolder<
      std::map<WignerKey,WignerMatrix_t,WignerKeyClassComp> >
      TheWignerDMatrixFactory;

    bool registerAll(const int Jmax); 

    // the most positive index is mapped to zero, next most is 1..
    WignerMatrix_t* call_factory(const mom_t &, const int); 

    template<int J> 
      WignerMatrix_t* call_factory(const mom_t &p)
      {
        return call_factory(p,J); 
      } 

  } // WignerDMatrixEnv 


} // radmat



#endif /* WIGNER_D_MATRIX_FACTORY_H */
