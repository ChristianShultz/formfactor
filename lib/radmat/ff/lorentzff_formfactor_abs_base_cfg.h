#ifndef LORENTZFF_FORMFACTOR_ABS_BASE_CFG_H
#define LORENTZFF_FORMFACTOR_ABS_BASE_CFG_H 


#include "formfactor_abs_base_cfg.h"
#include "formfactor_utils.h"
#include "radmat/utils/stringify.h"
#include "radmat/utils/pow2assert.h"
#include "radmat/utils/tensor.h"
#include "radmat/utils/handle.h"
#include "ensem/ensem.h"
#include <complex>
#include <utility>
#include <sstream>
#include <string>
#include <list>


namespace radmat
{


  // generate some linear least squares system w/o having to know very much about the 
  // matrix element invariants. These are still at the configuration level. Just 
  // shove this into a SembleMatrix bin by bin to deal with the ensemble stats.

  // this will be polymorphic with different named classes corresponding to 
  // different quantum numbers, ie: PiPi  <-->   < 0- | j_mu | 0- >

  struct LorentzFFAbsBase_t; 
  REGISTER_STRINGIFY_TYPE( LorentzFFAbsBase_t );

  struct LorentzFFAbsBase_t
    : public FFAbsBase_t
  {
    // save some typing
    typedef FFAbsBlockBase_t<std::complex<double> > BBType;
    typedef rHandle< BBType > BBHandle_t;
    typedef std::list< BBHandle_t > LorentzFFAbs_list;

    // this will be useful when we derive
    LorentzFFAbsBase_t(const LorentzFFAbs_list& list)
      : m_list(list) 
    {  }

    LorentzFFAbsBase_t& operator=(const LorentzFFAbsBase_t &o)
    {
      if(this != &o)
      {
        m_list = o.m_list;
      }
      return *this;
    }

    LorentzFFAbsBase_t(const LorentzFFAbsBase_t &o)
      : m_list(o.m_list)
    {  }

    // needs to be present and virtual b/c we are using pointers to derived
    virtual ~LorentzFFAbsBase_t(void) {}

    // useful higher up
    virtual int nFacs(void) const {return m_list.size();}

    // generate some tex code corresponding to thestring of stuff we think this is making
    virtual std::string ff(void) const 
    {
      std::stringstream ss;
      LorentzFFAbs_list::const_iterator it;
      for (it = m_list.begin(); it != m_list.end(); it++)
        ss << (*it)->ff() << "  ";
      return ss.str();
    }

    virtual std::map<int,std::string> ff_ids(void) const
    {
      std::map<int,std::string> ret; 
      int index = 0; 
      LorentzFFAbs_list::const_iterator it; 
      for(it = m_list.begin(); it != m_list.end(); ++it)
      {
        ret[index] = (*it)->id(); 
        ++index; 
      } 
      return ret; 
    }

    // matrix elem id
    virtual std::string reg_id() const { return Stringify<LorentzFFAbsBase_t>(); }
    virtual int left_spin() const = 0; 
    virtual int right_spin() const = 0; 

    // generate the linear system based on the available set of kinematic factors
    virtual itpp::Mat<std::complex<double> > 
      operator()( const MomRowPair_t &lefty, 
          const MomRowPair_t &righty,
          const double mom_fac) const
      {
        itpp::Mat<std::complex<double> > ret;

        LorentzFFAbs_list::const_iterator it;
        for (it = m_list.begin(); it != m_list.end(); it++)
          ret.append_col(toItpp<std::complex<double> >((**it)(lefty,righty,mom_fac)));

        return ret;
      }

    virtual LorentzFFAbsBase_t * clone() const = 0; 

    protected:  // hide ctor
    LorentzFFAbsBase_t(void);

    // data store
    LorentzFFAbs_list m_list;
  };


} // radmat 



#endif /* LORENTZFF_FORMFACTOR_ABS_BASE_CFG_H */
