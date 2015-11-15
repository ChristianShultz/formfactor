// test_load_fake_data.cc -
//
// Thursday, July 19 2012
//

#include"../headers/test_common.h"
#include "../headers/tester.h"
#include "radmat/load_data/load_fake_data.h"
#include "radmat/load_data/build_q2_packs.h"
#include "radmat/load_data/three_point.h"
#include "radmat/fake_data/fake_data_ini.h"
#include "radmat/llsq/llsq_driver.h"
#include "radmat/llsq/llsq_q2_pack.h"
#include "jackFitter/plot.h"
#include "adat/handle.h"
#include <complex>
#include <vector>
#include <string>

namespace radmat
{

  tester test_load_fake_data(void)
  {
    tester m_test(__func__);

    FakeDataIni_t ini = makeFakeIni();

    LoadFake3pt<std::complex<double> > loader(ini);

    std::vector<ThreePointCorrelator<std::complex<double> > > C3 = loader.genData();

    BuildQ2Packs<std::complex<double> > Q2Builder;
    Q2Builder.load(C3);
    Q2Builder.strip_propagation_factor();

    std::vector<ADAT::Handle<LLSQDataPointQ2Pack> > junk = Q2Builder.getQ2Packs();

    LLSQRet_ff_Q2Pack<std::complex<double> > foobar;
    LLSQDriver_t<std::complex<double> > llsq_driver(std::string("SVDMakeSquare"));

    foobar = *(llsq_driver(junk[0]));
/*
  LLSQRet_ff_Q2Pack<std::complex<double> >::const_iterator it;


  for(it = foobar.begin(); it != foobar.end(); ++it)
  {
    std::stringstream ss,ssc,ssr;
    ss << "blah" << it->first << ".ax";
    ssr << "blah_real" << it->first << ".ax"; 
    ssc << "blah_complex" << it->first << ".ax";

    ENSEM::EnsemVectorReal cplx,rl;

    rl = ENSEM::real(it->second);
    cplx = ENSEM::imag(it->second);

    AxisPlot plot,Rplot,Cplot;

    plot.addEnsemData(rl,"\\sq",1);
    plot.addEnsemData(cplx,"\\sq",2);
    Rplot.addEnsemData(rl,"\\sq",1);
    Cplot.addEnsemData(cplx,"\\sq",1);

    plot.sendToFile(ss.str());
    Rplot.sendToFile(ssr.str());
    Cplot.sendToFile(ssc.str());
      
  }
*/

  

    TESTER_TEST(m_test,true,"foobizzle");
    return m_test;
  }

} // namespace radmat
 
