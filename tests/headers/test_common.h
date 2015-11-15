#ifndef TEST_COMMON_H_H_GUARD
#define TEST_COMMON_H_H_GUARD

#include "tester.h"
#include "radmat/fake_data/fake_data_ini.h"
#include <string>
#include <complex>

namespace radmat
{
  // utils
  tester test_polarisation_tensor(void);
  tester test_ff_fake(void);

  // formfacs
  namespace PiPi
  {
    tester test_ff_debug(void);
    tester test_ff(void);
    tester test_ffKinematicFactor(void);
  }

  // llsq 
  tester test_LLSQ_solver_SVDMakeSquare(void);

  tester test_solver_factory(void);
  tester test_mat_elem_factory(void);

  FakeDataIni_t makeFakeIni(void);

  // fake data
  struct test_covarrying_vectors
  {
    template<typename T>
      tester test(void) const;
  };

  template<> tester test_covarrying_vectors::test<double>(void) const;
  template<> tester test_covarrying_vectors::test<std::complex<double> >(void) const;

  tester test_minimal_fake_data(const std::string &matElemID);
  tester test_make_fake_overlaps(void);
  tester test_make_fake_spectrum(void);
  tester test_fake_3pt_aux(void); 
  tester test_fake_3pt(void);
  tester test_gen_fake_dataset(void);
  tester test_load_fake_data(void);


  // readers
  tester test_readers(void);
  tester test_invert_subduction(void);
  tester test_xml_to_redstar(void);
}
#endif
