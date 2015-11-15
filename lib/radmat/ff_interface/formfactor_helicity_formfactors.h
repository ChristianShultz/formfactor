#ifndef FORMFACTOR_HELICITY_FORMFACTORS_H
#define FORMFACTOR_HELICITY_FORMFACTORS_H 




#include "formfactor.h"
#include "formfactor_factory.h"
#include "radmat/ff/lorentzff_formfactor_abs_base_cfg.h"
#include "radmat/utils/stringify.h"
#include "radmat/utils/handle.h"



namespace radmat
{

  struct HelicityFormFactorRecipe_t;
  REGISTER_STRINGIFY_TYPE( HelicityFormFactorRecipe_t ); 

  struct HelicityFormFactorRecipe_t
    : public FormFactorRecipe_t
  {
    typedef rHandle<LorentzFFAbsBase_t> mat_h; 
    typedef rHandle<LorentzRep> rep_h; 

    HelicityFormFactorRecipe_t(const mat_h &f, const rep_h &l, const rep_h &r)
    {
      mat = rHandle<LorentzFFAbsBase_t>(f->clone()); 
      lefty = rHandle<LorentzRep>(l->clone()); 
      righty = rHandle<LorentzRep>(r->clone()); 
    }

    HelicityFormFactorRecipe_t( const HelicityFormFactorRecipe_t &o)
    {
      mat = rHandle<LorentzFFAbsBase_t>(o.mat->clone()); 
      lefty = rHandle<LorentzRep>(o.lefty->clone()); 
      righty = rHandle<LorentzRep>(o.righty->clone()); 
    }

    HelicityFormFactorRecipe_t& operator=( const HelicityFormFactorRecipe_t &o)
    {
      if (this != &o) 
      {
        mat = rHandle<LorentzFFAbsBase_t>(o.mat->clone()); 
        lefty = rHandle<LorentzRep>(o.lefty->clone()); 
        righty = rHandle<LorentzRep>(o.righty->clone()); 
      }
      return *this;
    }

    virtual ~HelicityFormFactorRecipe_t() {}

    virtual int nFacs() const { return mat->nFacs(); }
    virtual std::string ff() const { return mat->id(); }
    virtual std::map<int,std::string> ff_ids() const { return mat->ff_ids(); }
    virtual rHandle<Rep_p> left_rep() const { return FormFactorRecipe_t::call(lefty->rep_id()); }
    virtual rHandle<Rep_p> right_rep() const { return FormFactorRecipe_t::call(righty->rep_id()); }
    virtual std::string id() const { return Stringify<HelicityFormFactorRecipe_t>(); }
    virtual std::string reg_id() const { return mat->reg_id(); } // use the same names as in lorentzff

    virtual FormFactorRecipe_t * clone() const 
    {
      return new HelicityFormFactorRecipe_t( mat,lefty,righty); 
    }

    rHandle<LorentzFFAbsBase_t> mat; 
    rHandle<LorentzRep> lefty, righty; 
  };


  struct HelicityFormFactor; 
  REGISTER_STRINGIFY_TYPE( HelicityFormFactor ); 

  struct HelicityFormFactor
    : public FormFactorBase_t
  {
    typedef FormFactorBase_t::recipe_h recipe_h;  

    HelicityFormFactor( const recipe_h & r)
      : FormFactorBase_t(r)
    { }

    HelicityFormFactor& operator=(const HelicityFormFactor &o)
    {
      FormFactorBase_t::operator=(o); 
    }

    HelicityFormFactor( const HelicityFormFactor &o)
      : FormFactorBase_t(o)
    { }

    virtual ~HelicityFormFactor() {}

    virtual itpp::Mat<std::complex<double> >
      generate_ffs( const MomRowPair_t &lefty, 
          const MomRowPair_t &righty, 
          const double mom_fac) const; 
  };


  namespace HelicityFormFactorDecompositionFactoryEnv
  {
    bool registerAll(void) ;
    std::string build_id( const rHandle<LorentzRep> &lefty,
        const rHandle<LorentzRep> &righty ); 
  }


} // radmat

#endif /* FORMFACTOR_HELICITY_FORMFACTORS_H */
