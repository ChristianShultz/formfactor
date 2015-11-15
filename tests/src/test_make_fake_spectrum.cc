// test_make_fake_spectrum.cc -
//
// Monday, July  2 2012
//

#include "../headers/tester.h"
#include "../headers/test_common.h"
#include "radmat/fake_data/fake_spectrum.h"
#include "radmat/fake_data/fake_data_ini.h"
#include "semble/semble_matrix.h"
#include <vector>
#include <string>
#include <iostream>


namespace radmat
{

//#define FAKE_SPECTRUM_DEBUG

  tester test_make_fake_spectrum(void)
  {
    tester m_test(__func__);
    TESTER_TEST(m_test,true,"foobar");

    Array<double> spectrum; 
    spectrum.resize(3);
    spectrum[0] = spectrum[1] = spectrum[2] = acos(-1.)/2.;


    FakeDataIni_t fakeIni;
    fakeIni.dataProps.ncfg = 50;
    fakeIni.timeProps.tsink = 5;
    fakeIni.timeProps.tsource = 0;
    fakeIni.stateProps.mProps.sinkUpdateCovariance = true;
    fakeIni.stateProps.mProps.sinkUpdateVariance = true;
    fakeIni.stateProps.mProps.sinkVarO = 0.1;
    fakeIni.stateProps.mProps.sinkMasses = spectrum;

    FakeSpectrum fakeSpec(fakeIni);

    std::vector<SEMBLE::SembleVector<double> > foo = fakeSpec.generate(std::string("sink"));


#ifdef FAKE_SPECTRUM_DEBUG 
  
    std::cout << "vector.size() " << foo.size() << std::endl; 
    std::cout << "getB() " << foo[0].getB() << std::endl;
    std::cout << "getN() " << foo[0].getN() << std::endl;
    std::vector<SEMBLE::SembleVector<double> >::const_iterator it;
    
  for(it = foo.begin(); it != foo.end(); it++)
    std::cout << it->mean() << std::endl;

#endif


    return m_test;
  }


#undef FAKE_SPECTRUM_DEBUG
}
