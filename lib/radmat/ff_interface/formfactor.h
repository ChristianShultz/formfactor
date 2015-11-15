#ifndef FORMFACTOR_H
#define FORMFACTOR_H 

#include "radmat/ff/formfactor_abs_base_cfg.h"
#include "radmat/data_representation/data_representation.h"
#include "radmat/utils/stringify.h"

namespace radmat
{
  // an abstract recipe that gets implemented
  // differently for helicity and cubic ffs
  //
  // this allows us to inject different parameters
  // to the form factors upon instantiation depending
  // on the "type" eg. helicity or cubic

  struct FormFactorRecipe_t; 
  REGISTER_STRINGIFY_TYPE( FormFactorRecipe_t ); 

  struct FormFactorRecipe_t
  {
    virtual ~FormFactorRecipe_t() {}
    virtual int nFacs() const = 0;
    virtual std::string ff() const = 0; 
    virtual std::map<int,std::string> ff_ids(void) const = 0; 
    virtual rHandle<Rep_p> left_rep() const = 0; 
    virtual rHandle<Rep_p> right_rep() const = 0; 
    virtual FormFactorRecipe_t* clone() const = 0; 
    virtual rHandle<Rep_p> call(const std::string &s) const
    {
      return ::radmat::DataRepresentationFactoryEnv::callFactory(s);
    }
    virtual std::string id() const { return Stringify<FormFactorRecipe_t>(); }
  };


  // an abstract base that allows us to inject a recipe 
  // upon instantiation 
  //
  // the recipe is essentially a list of coefficients and 
  // helicities that allow us to subduce to the cubic group
  //
  // for helicity form factors the list is a single element
  // with coefficient 1. corresponding to whatever we asked for
  //
  // for cubic form factors the list is the subduction recipe
  struct FormFactorBase_t;
  REGISTER_STRINGIFY_TYPE( FormFactorBase_t ); 

  struct FormFactorBase_t
    : public FFAbsBase_t 
  {
    typedef rHandle<FormFactorRecipe_t> recipe_h; 

    FormFactorBase_t( const recipe_h recipe_)
    {
      recipe = recipe_->clone(); 
    } 

    FormFactorBase_t(const FormFactorBase_t &o)
    {
      delete recipe; 
      recipe = o.recipe->clone(); 
    }

    virtual ~FormFactorBase_t() 
    {
      delete recipe; 
    }

    // we now have to sum over a list
    virtual itpp::Mat<std::complex<double> > 
      operator()(const MomRowPair_t &lefty, 
          const MomRowPair_t &righty, 
          const double mom_fac) const
      {
        return generate_ffs(lefty,righty,mom_fac); 
      }

    // may need to be re implemented in derived
    virtual int nFacs() const { return recipe->nFacs(); }
    virtual std::string ff() const { return recipe->ff(); }
    virtual std::map<int,std::string> ff_ids() const { return recipe->ff_ids(); }
    virtual std::string id() const { return Stringify<FormFactorBase_t>(); }
    virtual rHandle<Rep_p> left_rep() const { return recipe->left_rep(); }
    virtual rHandle<Rep_p> right_rep() const { return recipe->right_rep(); }

    // must be implemented in derived
    virtual itpp::Mat<std::complex<double> >
      generate_ffs( const MomRowPair_t &lefty, 
          const MomRowPair_t &righty, 
          const double mom_fac) const = 0; 

     recipe_h get_recipe() const { return recipe_h(recipe->clone()); }

     FormFactorRecipe_t* recipe; 
  };

} // radmat

#endif /* FORMFACTOR_H */
