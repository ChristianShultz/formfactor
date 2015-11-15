/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : fit_tins.cc

 * Purpose :

 * Creation Date : 01-08-2012

 * Last Modified : Mon 18 Aug 2014 06:56:58 PM EDT

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/


#include "fit_tins.h"
#include "ensem/ensem.h"
#include "semble/semble_semble.h"
#include "jackFitter/ensem_data.h"
#include "jackFitter/jackknife_fitter.h"
#include "jackFitter/ensem_data.h"
#include "jackFitter/three_point_fit_forms.h"
#include "jackFitter/plot.h"
#include <string>
#include <math.h>
#include <vector>
#include <iostream>
#include <complex>
#include "radmat/utils/splash.h"
#include "radmat/utils/pow2assert.h"
#include "radmat/utils/printer.h"
#include "adat/handle.h"




namespace radmat
{

  namespace 
  {

    template<typename T>
      std::string to_string( T t )
      {
        std::stringstream ss; 
        ss << t ;
        return ss.str(); 
      }


    struct to_ensem_printer
    {
      static void print(const std::string &s)
      {} //  { std::cout << s << std::endl; }
    };

    struct dimension_printer
    {
      static void print(const std::string &s)
      {} //  { std::cout << s << std::endl; }
    };

    template<typename T> 
      typename SEMBLE::PromoteEnsemVec<T>::Type
      to_ensem(const SEMBLE::SembleVector<T> &in)
      {
        typename SEMBLE::PromoteEnsemVec<T>::Type out;
        printer_function<to_ensem_printer>( "B=" + to_string(in.getB())); 
        printer_function<to_ensem_printer>( "N=" + to_string(in.getN())); 
        printer_function<to_ensem_printer>( "in.mean()=" + to_string(in.mean())); 
        out.resize(in.getB());
        out.resizeObs(in.getN());
        for(int elem = 0; elem < in.getN(); ++elem)
          ENSEM::pokeObs(out,in.getEnsemElement(elem),elem);

        return out;
      }
  }


  void TinsFitter::fit(const std::string &fname, 
      const LLSQRealFormFactorData_t &data, 
      const ThreePointComparatorProps_t &fitProps,
      const int tsrc, 
      const int tsnk,
      const bool usePars,
      const FitParValue fp)
  {
    const unsigned int sz = data.size();
    const int nbins = data.esize();

    printer_function<dimension_printer>(
        "fit dimensions, nff = " + to_string(sz)
        + " nbins = " + to_string(nbins) );   

    ff.reDim(nbins,sz);
    Q2 = data.Q2();

    LLSQRealFormFactorData_t::const_iterator it;

    for(it = data.begin(); it != data.end(); ++it)
      doFit(fname,to_ensem(it->second),it->first,fitProps,tsrc,tsnk,usePars,fp);

    didFit = true;
  }

  void TinsFitter::doFit(const std::string &filenameBase, 
      const ENSEM::EnsemVectorReal &data, 
      const std::string &ffid,
      const ThreePointComparatorProps_t &fitProps,
      const int tsrc,
      const int tsnk,
      const bool usePars,
      const FitParValue fp)
  {
    //  std::cout << __PRETTY_FUNCTION__ << ": entering " << std::endl;
    //  std::cout << "tsrc " << tsrc << "   tsnk " << tsnk << std::endl;

    std::stringstream ss,jack,ax,fit;
    ss << filenameBase;
    fit << ss.str() << ffid << "_fit.jack";
    jack << ss.str() << ffid << ".jack";
    ax << ss.str() << ffid << ".ax";

    std::vector<double> time; 
    const int Lt = data.numElem();

    for(int t = 0; t < Lt; t++)
      time.push_back(t);

    EnsemData corrData(time,data);
    corrData.setSVCutoff(fitProps.SVCutOff); 

    ADAT::Handle<FitComparator> fitComp = constructThreePointFitComparator(fitProps);
    // NB: I Have assumed that no chopping has gone on in the data
    rHandle<FitThreePoint> fitCorr (new FitThreePoint(corrData,tsnk,tsrc,
          fitProps.thigh,fitProps.tlow,fitComp,fitProps.minTSlice,fitProps.fit_type,usePars,fp));

    // send to files
    fitCorr->saveFitPlot(ax.str());
    write(jack.str(),data);

    // require unique fit names -- they are class names 
    // so this is a paranoia check -- we dont want to find it!!!
    POW2_ASSERT( fit_index.find(ffid) == fit_index.end());

    // this grows the map effectively incrementing the index!
    int index = fit_index.size(); 
    fit_index.insert(std::make_pair(ffid,index)); 

    printer_function<dimension_printer>(
        "inserting a fit into index " + to_string(index));

    //  ff.dimensions(__func__); 

    ff.loadEnsemElement(index,fitCorr->getFF());
    write(fit.str(),fitCorr->getFF()); 

    fitters.insert(std::map<std::string,rHandle<FitThreePoint> >::value_type(ffid,fitCorr));
  }


  void TinsFitter::writeFitLogs(const std::string &path) const
  {
    std::map<std::string,rHandle<FitThreePoint> >::const_iterator it;
    for(it = fitters.begin(); it != fitters.end(); ++it)
    {
      std::stringstream ss;
      ss << path << it->first << ".fitlog.summary";
      std::ofstream out(ss.str().c_str());
      out << it->second->getFitSummary();
      out.close();
    }

  }

  void TinsFitter::writeFitPlotsWithComponents(const std::string & path) const
  {
    std::map<std::string,rHandle<FitThreePoint> >::const_iterator it;
    for(it = fitters.begin(); it != fitters.end(); ++it)
    {
      std::stringstream ss; 
      ss << path << it->first << ".ax";
      std::ofstream out(ss.str().c_str());
      out << it->second->getFitPlotStringWithComponents(); 
      out.close();
    }

  }


  std::pair<ENSEM::EnsemReal,SEMBLE::SembleVector<double> > TinsFitter::fetchFF(void) const
  {
    POW2_ASSERT(didFit);

    return std::pair<ENSEM::EnsemReal,SEMBLE::SembleVector<double> >(Q2,ff);
  }

  std::pair<ENSEM::EnsemReal,std::map<std::string,ENSEM::EnsemReal> > TinsFitter::fetchNamedFF(void) const
  {
    POW2_ASSERT(didFit); 
    std::vector<std::string> ids = ff_ids(); 
    std::vector<std::string>::const_iterator it; 
    std::map<std::string,ENSEM::EnsemReal> mappy; 
    for(it = ids.begin(); it != ids.end(); ++it)
      mappy.insert( std::make_pair( *it, getFF(*it) ) ); 

    return std::make_pair( Q2, mappy ); 
  }


  ENSEM::EnsemReal TinsFitter::getFF(const std::string &id) const
  {
    return ff.getEnsemElement(fit_index.at(id)); // throws out_of_range
  }

  rHandle<FitThreePoint> TinsFitter::getFit(const std::string &id) const
  {
    return rHandle<FitThreePoint>(fitters.at(id)); // throws out_of_range
  }

  std::vector<std::string> TinsFitter::ff_ids(void) const
  {
    std::vector<std::string> k; 
    std::map<std::string,int>::const_iterator it;
    for(it = fit_index.begin(); it != fit_index.end(); ++it)
      k.push_back(it->first);
    return k; 
  }


} // namespace radmat
