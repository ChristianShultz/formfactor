#ifndef FORMFACTOR_FACTORY_H
#define FORMFACTOR_FACTORY_H 


#include "formfactor.h"
#include <string>
#include <complex>
#include <vector>
#include "adat/singleton.h"
#include "adat/objfactory.h"
#include "radmat/utils/handle.h"

//
//
//    helicities are a direct map back to cont_spin ffs
//    cubics are weighted sums of cont_spin ffs
//
//


namespace radmat
{

  struct PtrRecipeHolder
  {
    typedef std::map<std::string, FormFactorRecipe_t *> map_t; 
    ~PtrRecipeHolder() 
    {
      map_t::iterator it; 
      for(it = mappy.begin(); it != mappy.end(); ++it)
        delete it->second; 
    }

    map_t mappy; 
  };

  // use PtrRecipeHolder since this can be threadded 
  typedef Util::SingletonHolder< PtrRecipeHolder >
    TheFormFactorInjectionRecipeCookbook; 


  // inject recipes upon instantiation 
  typedef Util::SingletonHolder<
    Util::ObjectFactory<FormFactorBase_t,
    std::string,
    TYPELIST_1(const std::string &),
    FormFactorBase_t* (*)(const std::string &),
    Util::StringFactoryError> >
      TheFormFactorInjectionRecipeFactory;


  // return a recipe injected form factor
  namespace FormFactorDecompositionFactoryEnv
  {
    bool registerAll( void );
    rHandle<FormFactorBase_t> callFactory(const std::string &matElemID);
    std::vector<std::string> all_keys(void); 
  }

}



#endif /* FORMFACTOR_FACTORY_H */
