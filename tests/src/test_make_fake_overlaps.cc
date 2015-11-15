// test_make_fake_overlaps.cc -
//
// Monday, July  2 2012
//

#include "../headers/tester.h"
#include "../headers/test_common.h"
#include "radmat/fake_data/fake_overlaps.h"
#include "radmat/fake_data/fake_data_ini.h"
#include "semble/semble_matrix.h"
#include <complex>
#include <vector>
#include <string>
#include <iostream>

namespace radmat
{

//#define DEBUG_FAKE_LAPS_CJS 

  tester test_make_fake_overlaps(void) 
  {
    tester m_test(__func__);
    TESTER_TEST(m_test,true,"foobar"); // just need it to go into this function so a fake test is present

    FakeDataIni_t fakeIni;
    fakeIni.stateProps.zProps.overlapGenerator = std::string("unitary");
    fakeIni.dataProps.ncfg = 50;
    fakeIni.stateProps.zProps.updateCovariance = true;
    fakeIni.stateProps.zProps.updateVariance = true;
    fakeIni.stateProps.zProps.varianceOrder = 0.1;
    fakeIni.stateProps.readZ = false;
    fakeIni.timeProps.tsink = 5;
    fakeIni.timeProps.tsource = 0;

    FakeOverlaps fakeLaps(fakeIni);

    std::vector<SEMBLE::SembleMatrix<double> > doublefoo = fakeLaps.generate<double>(5,std::string("source"));
    std::vector<SEMBLE::SembleMatrix<std::complex<double> > > complexfoo = fakeLaps.generate<std::complex<double> >(5,std::string("source"));

    std::vector<SEMBLE::SembleVector<std::complex<double> > > bar = pullRow(complexfoo,0);

#ifdef DEBUG_FAKE_LAPS_CJS

    //  std::vector<SEMBLE::SembleMatrix<std::complex<double> > >::const_iterator bit;
    //  for(bit = complexfoo.begin(); bit != complexfoo.end(); bit++)
    //    std::cout << (*bit)[0] << std::endl;

    std::vector<SEMBLE::SembleVector<std::complex<double> > >::const_iterator it;
    for(it = bar.begin(); it != bar.end(); it++)
      std::cout << it->mean() << std::endl;
#endif

    return m_test;
  }



#undef DEBUG_FAKE_LAPS_CJS

}
