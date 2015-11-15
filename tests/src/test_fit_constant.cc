// test_fit_constant.cc -
//
// Monday, July 30 2012
//

#include "../headers/test_common.h"
#include "../headers/tester.h"
#include "radmat/fitting/fit_three_point_ff.h"
#include "radmat/fitting/radmat_fit_forms.h"
#include "radmat/fake_data/covarrying_vectors.h"
#include "ensem/ensem.h"
#include "semble/semble_meta.h"
#include "semble/semble_vector.h"
#include "jackFitter/ensem_data.h"
#include "jackFitter/jackknife_fitter.h"
#include "adat/handle.h"
#include "itpp/itbase.h"
#include <stdlib.h>
#include <time.h>
#include <string>

namespace radmat
{

  bool test_ensem_equals(const ENSEM::EnsemReal &e, const double target, const double prec)
  {

    // even here the ensem stuff wants to be a pain.. sigh
    for(int cfg = 0; cfg < e.size(); ++cfg)
      if (fabs((SEMBLE::toScalar(e.elem(cfg)) - target)/target) > prec )
        return false;

    return true;
  }

  tester test_fit_constant(void)
  {
    tester m_test(__func__);

    srand(time(NULL));
    itpp::RNG_reset(rand()); // reseed itpp rng with time

    double mean = 4.;
    double varO = 0.01;
    const int ncfg = 300;
    const int nObs = 25;

    itpp::Vec<double> vmean(nObs);
    itpp::Vec<double> vvar = itpp::abs(itpp::randn(nObs)) * varO;
    corMat c;
    itpp::Mat<double> correlation = c.genMat<double>(nObs);

    vmean = mean;

    SEMBLE::SembleVector<double> ensemvec(genCovarryingDist(vmean,vvar,correlation,ncfg));


    ENSEM::EnsemVectorReal corr;
    corr.resize(ncfg);
    corr.resizeObs(nObs);

    std::vector<double> time;


    for(int i =0; i < nObs; i++)
    {
      ENSEM::pokeObs(corr,ensemvec.getEnsemElement(i),i);
      time.push_back(i);
    }

    EnsemData corrData(time,corr);
    ADAT::Handle<FitComparator> fitComp(new CompareFitsByChisqPerNDoF);


    TESTER_TEST(m_test,&*fitComp,"can't allocate a Constant");


    Fit3PtFF_Constant fitCorr(corrData,fitComp,0.1,5);

    ENSEM::EnsemReal ff = fitCorr.getFF();

    bool eq = test_ensem_equals(ff,mean,2.*varO);

    TESTER_TEST(m_test, eq, "couldn't fit a constant");

    fitCorr.saveFitPlot(std::string("test_fit_plot.ax"));

    return m_test;
  }



} // namespace radmat
