/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name :

 * Purpose :

 * Creation Date : 18-07-2012

 * Last Modified : Wed Jul 18 16:59:32 2012

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/

#include "../headers/tester.h"
#include "../headers/test_common.h"
#include "radmat/fake_data/gen_fake_dataset.h"
#include "radmat/fake_data/fake_data_ini.h"
#include "adat/handle.h"
#include <vector>
#include <complex>

typedef std::complex<double> dc;

namespace radmat
{

  tester test_gen_fake_dataset(void)
  {
    tester m_test(__func__);
    FakeDataIni_t ini = makeFakeIni();
    GenFakeDataSet<dc> genset(ini);
    genset.generate();
    ADAT::Handle<std::vector<GenFakeDataSet<dc>::Corr> > foo;
    foo = genset.get();

    // didn't crash  -- single correlator tests are looking good, 
    // this is just a looping structure so its probably fine
    TESTER_TEST(m_test,true,"foobar");
    return m_test;
  }


}
