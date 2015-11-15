#ifndef SPECTRUM_FACTORY_H
#define SPECTRUM_FACTORY_H 

#include "spectrum_state.h"
#include "adat/singleton.h"
#include "adat/objfactory.h"
#include "radmat/utils/utils.h"


namespace radmat
{
  typedef Util::SingletonHolder<
    Util::ObjectFactory<SpectrumState_p,
    std::string, 
    void, 
    SpectrumState_p* (*)(void), 
    Util::StringFactoryError> > 
      TheSpectrumFactory; 


  namespace TheSpectrumFactoryEnv
  {
    bool registerAll(const int enumLattice); 
    rHandle<SpectrumState_p> callFactory(const std::string &id); 
    std::pair<double,double> pole_enhancement(const std::string &l, 
        const ADATXML::Array<int> &lmom, 
        const std::string &r, 
        const ADATXML::Array<int> &rmom, 
        const std::string &pole); 
  }

}

#endif /* SPECTRUM_FACTORY_H */
