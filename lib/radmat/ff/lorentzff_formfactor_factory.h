#ifndef LORENTZFF_FORMFACTOR_FACTORY_H_H_GUARD
#define LORENTZFF_FORMFACTOR_FACTORY_H_H_GUARD

#include "formfactor_abs_base_cfg.h"
#include <string>
#include <complex>
#include <vector>
#include "adat/singleton.h"
#include "adat/objfactory.h"
#include "radmat/utils/handle.h"

//
//  Hold the canonical_cont_spin_formfactors
//
//    these are then split off into helicity 
//    and cubic variants
//
//    helicities are a direct map back to cont_spin ffs
//    cubics are weighted sums of cont_spin ffs
//
//    everything in this factory is effectively hidden 
//    from the rest of radmat 
//


namespace radmat
{
  typedef Util::SingletonHolder<
    Util::ObjectFactory<FFAbsBase_t,
			std::string,
			void,
			FFAbsBase_t* (*)(void),
			Util::StringFactoryError> >
  TheLorentzffFormFactorDecompositionFactory;


  namespace LorentzffFormFactorDecompositionFactoryEnv
  {
    bool registerAll( void );
    rHandle<FFAbsBase_t> callFactory(const std::string &matElemID);
    std::vector<std::string> all_keys(void); 
  }

}

#endif
