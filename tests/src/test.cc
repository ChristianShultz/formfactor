// main.cc -
//
// Tuesday, May  1 2012
//

#include <omp.h>
#include <iostream>
#include "../headers/tester.h"
#include "../headers/test_tensor.h"
#include "../headers/test_common.h"

//#define DO_TESTS_PARALLEL

using namespace std;
using namespace radmat;

void dump_test(unsigned short &ct_err, unsigned short &ct_test,const tester &t)
{
#pragma omp critical
  {
    ct_err += t.ct_err;
    ct_test += t.ct_test;
    cout << t.test_unit << " " << t.ct_test - t.ct_err << " tests of " << t.ct_test << " possible passed" << endl;
  }
}


  int 
main(void)
{

  // NB: factories don't like to be parallelized, needs to be reorganized for can and can't parallel sections.. oops
#ifndef DO_TESTS_PARALLEL
  omp_set_num_threads(0);
#endif

  // test the tester
  tester  test_test(__func__);
  const char * msg = "tested false successfully";
  cout << "\n\nyou should see the following message\n\"" << __func__ << " : line # : "<< msg << "\""<< endl;
  TESTER_TEST(test_test,false,msg);

  // begin the actual tests
  cout << "\n\n" << "beginning tests.." << endl;
  unsigned short ct_err(0), ct_test(0);

#pragma omp parallel
  {  
#pragma omp sections
    {
      /* 
#pragma omp section
{
      // test the tensors
      dump_test(ct_err,ct_test,test_tensor(double(10)));
      dump_test(ct_err,ct_test,test_tensor(std::complex<double>(10,10)));
      }

#pragma omp section
{
      // test out polarisation tensors 
      dump_test(ct_err,ct_test,test_polarisation_tensor());
      }

#pragma omp section
{
      // test the PiPi form factor
      dump_test(ct_err,ct_test,PiPi::test_ff());
      dump_test(ct_err,ct_test,PiPi::test_ffKinematicFactor());
      }

#pragma omp section
{
      // test the factories
      dump_test(ct_err,ct_test,test_solver_factory());
      dump_test(ct_err,ct_test,test_mat_elem_factory());
      }

#pragma omp section
{
      // test the LLSQ framework
      dump_test(ct_err,ct_test,test_LLSQ_solver_SVDMakeSquare());
      }

    // fake data
#pragma omp section
{
test_covarrying_vectors cov;
dump_test(ct_err,ct_test,cov.test<double>());
dump_test(ct_err,ct_test,cov.test<std::complex<double> >());
}

#pragma omp section
{
dump_test(ct_err,ct_test,test_minimal_fake_data(std::string("PiPi")));
}

#pragma omp section
{
dump_test(ct_err,ct_test,test_make_fake_overlaps());
}

#pragma omp section
{
dump_test(ct_err,ct_test,test_make_fake_spectrum());
}

#pragma omp section
{
dump_test(ct_err,ct_test,test_fake_3pt_aux());
}

#pragma omp section
{
dump_test(ct_err,ct_test,test_fake_3pt());
}

#pragma omp section
{
dump_test(ct_err,ct_test,test_gen_fake_dataset()); //-- redundant, load fake data makes a fake set
}

#pragma omp section
{
  dump_test(ct_err,ct_test,test_load_fake_data());
}

#pragma omp section
{
  dump_test(ct_err,ct_test,test_fit_constant());
}
*/
#pragma omp section
{
//  dump_test(ct_err,ct_test,test_readers());
}

#pragma omp section
{
//  dump_test(ct_err, ct_test, test_invert_subduction());
}
#pragma omp section
{
  dump_test(ct_err, ct_test, test_xml_to_redstar());
}

} // omp sections
} // parallel 


// conclude testing
cout << "\n**********************************\n"
<< ct_test - ct_err << " of " << ct_test << " tests passed.\n" << endl;
cout << "finished testing" << endl;


if(ct_err != 0) 
  cout << "Warning: " << ct_err << " tests failed." << endl;

  return 0;
  }
