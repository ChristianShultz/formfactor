#ifndef LLSQ_FORMFACTOR_DATA_H
#define LLSQ_FORMFACTOR_DATA_H 


#include "semble/semble_semble.h"
#include <map>
#include <string>
#include "ensem/ensem.h"
#include "radmat/utils/printer.h"

namespace radmat
{

  ////////////////////////////////////////////////////
  ////////////////////////////////////////////////////
 
 template< typename T>  
  struct LLSQFormFactorData_t
  {
    typedef typename std::map<std::string,SEMBLE::SembleVector<T> > map_t; 
    typedef typename map_t::iterator iterator; 
    typedef typename map_t::const_iterator const_iterator; 
    typedef typename map_t::value_type value_type;

    LLSQFormFactorData_t() 
    {
      Qsq.resize(1); 
    }

    LLSQFormFactorData_t(const LLSQFormFactorData_t &d)
    { 
      Qsq = d.Qsq; 
      mappy = d.mappy; 
    }

    void insert( const std::string &s, const SEMBLE::SembleVector<T> &v)
    {

      if(mappy.find(s) != mappy.end())
      {
        std::cerr << __PRETTY_FUNCTION__ << "double insert error" << std::endl;
        exit(1); 
      }

      mappy.insert( value_type(s,v) ); 
    }

    iterator begin() { return mappy.begin(); }
    const_iterator begin() const { return mappy.begin(); }
    iterator end() { return mappy.end(); }
    const_iterator end() const { return mappy.end(); }
    ENSEM::EnsemReal Q2() const { return Qsq; }
    unsigned int size() const { return mappy.size(); }
    int esize() const { return mappy.begin()->second.getB(); }
    
    ENSEM::EnsemReal Qsq; 
    map_t mappy; 
  }; 

  ////////////////////////////////////////////////////
  ////////////////////////////////////////////////////

 typedef LLSQFormFactorData_t<std::complex<double> > LLSQComplexFormFactorData_t; 
 typedef LLSQFormFactorData_t<double> LLSQRealFormFactorData_t; 

 // there can be an arbitray overall phase, get rid of it
 LLSQRealFormFactorData_t 
   rephase_formfactor_data( const LLSQComplexFormFactorData_t &, 
       const int tlow, 
       const int thigh, 
       const std::string &filebase); 

} // radmat

#endif /* LLSQ_FORMFACTOR_DATA_H */
