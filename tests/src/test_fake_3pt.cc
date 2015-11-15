/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

* File Name : test_fake_3pt.cc

* Purpose :

* Creation Date : 16-07-2012

* Last Modified : Mon Jul 16 12:08:10 2012

* Created By : shultz

_._._._._._._._._._._._._._._._._._._._._.*/


#include "headers/test_common.h"
#include "headers/tester.h"
#include "radmat/fake_data/fake_3pt_function_aux.h"
#include "radmat/fake_data/fake_3pt_function.h"
#include "radmat/fake_data/fake_data_ini.h"
#include "semble/semble_meta.h"
#include <complex>

namespace radmat
{

  tester test_fake_3pt(void)
  {
    tester m_test(__func__);

    FakeDataIni_t ini = makeFakeIni();


    ADAT::Handle<FakeDataInputs_p<std::complex<double> > > 
      orig = generateOriginalInputs<std::complex<double> >(ini);
    ADAT::Handle<FakeDataInputs<std::complex<double> > > input = copyFakeInput(orig);
    applyZSuppression(input->working);
    applyDispersion(input->working,ini.dataProps.momenta[0]);

    Fake3ptCorr<std::complex<double> > foo(input,0,0,0,ini.dataProps.momenta[0]);

// made it through 
    TESTER_TEST(m_test,true,"foo");

    return m_test;
  }


}
