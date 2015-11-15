#ifndef FIT_TINS_H_H_GUARD
#define FIT_TINS_H_H_GUARD

#include <utility>
#include <string>
#include <complex>
#include <vector>
#include "ensem/ensem.h"
#include "semble/semble_vector.h"
#include "radmat/utils/handle.h"
#include "radmat/llsq/llsq_formfactor_data.h"
#include "jackFitter/three_point_fit_forms.h"

namespace radmat
{


  // find the flat region, figure out if we want the real or imag part
  struct TinsFitter
  {

    // construct
    TinsFitter(void)
      : didFit(false) 
    { 
      Q2.resize(1);
      Q2 = ENSEM::toDouble(10000000.);
    }


    void fit( const std::string &filenameBase, 
        const LLSQComplexFormFactorData_t &data,
        const ThreePointComparatorProps_t &fitProps, 
        const int tsrc,
        const int tsnk,
        const bool usePars=false,
        const FitParValue fp=FitParValue())
    {
      fit(filenameBase, rephase_formfactor_data( data , fitProps.tlow, fitProps.thigh, filenameBase ), fitProps, tsrc, tsnk,usePars,fp); 
    }

    void fit( const std::string &filenameBase, 
        const LLSQRealFormFactorData_t &data,
        const ThreePointComparatorProps_t &fitProps, 
        const int tsrc,
        const int tsnk,
        const bool usePars=false,
        const FitParValue=FitParValue());

    // get the form factors at this q2
    std::pair<ENSEM::EnsemReal, SEMBLE::SembleVector<double> > fetchFF(void) const;

    // get the named form factors at this q2
    std::pair<ENSEM::EnsemReal, std::map<std::string,ENSEM::EnsemReal> > fetchNamedFF(void) const; 

    // get all of the ids
    std::vector<std::string> ff_ids(void) const; 

    // get a single form factor at this q2
    ENSEM::EnsemReal getFF(const std::string &ffid) const; 

    // get the fit associated with ffnum
    rHandle<FitThreePoint> getFit(const std::string &ffid) const;

    // get the ensemble value for q2
    ENSEM::EnsemReal getQ2(void) const {return Q2;}

    // how many form factors are there
    int nFF(void) const {return ff.getN();}

    // write out the fit log files 
    void writeFitLogs(const std::string &path) const;

    // write the fits with the components also plotted
    void writeFitPlotsWithComponents(const std::string &path) const; 

    private:

    // plays nicely with the "fit" function
    void doFit(const std::string &filenameBase, 
        const ENSEM::EnsemVectorReal &data, 
        const std::string &ffid,
        const ThreePointComparatorProps_t &fitProps,
        const int tsrc, 
        const int tsnk,
        const bool usePars=false,
        const FitParValue=FitParValue());

    // data store
    bool didFit;
    std::map<std::string,rHandle<FitThreePoint> > fitters; 
    std::map<std::string,int> fit_index; 
    ENSEM::EnsemReal Q2;
    SEMBLE::SembleVector<double> ff;
  };


} // namespace radmat




#endif
