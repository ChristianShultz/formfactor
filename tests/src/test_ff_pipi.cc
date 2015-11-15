// test_ff_pipi.cc -
//
// Friday, June  1 2012
//


#include "../headers/tester.h"
#include "../headers/test_common.h"
#include "itpp/itbase.h"
#include "radmat/ff/lorentzff_PiPi.h"
#include "radmat/ff/ff_gen_llsq_row.h"
#include "radmat/ff/formfactor_factory.h"
#include "radmat/utils/tensor.h"
#include "radmat/utils/pow2assert.h"
#include <complex>
#include <utility>
#include <iostream>

namespace radmat
{
  namespace PiPi
  {
    
    
    tester test_ff(void)
    {
      return test_ff_debug(); // lazy
    }


    // --DEBUGGING TEST--
    tester test_ff_debug(void)
    {
      tester m_test(__func__);

      for(int i = 0; i < 200; i++)
	{
	  PiPi foobar;
	  Tensor<double,1> foo((TensorShape<1>())[4],0.) , bar((TensorShape<1>())[4],0.) ;
	  Tensor<std::complex<double> , 1> baz((TensorShape<1>())[4],0.);
      
	  for (idx_t i = 0; i < 4; i++)
	    {
	      double tmp = 2*(itpp::randu() - .5);
	      double tmp2 = 2*(itpp::randu() - .5);
	      foo[i] = tmp;
	      bar[i] = tmp2;
	      baz[i] = tmp + tmp2;
	    }

	  
	  itpp::Mat<std::complex<double> > bazMat = toItpp(baz);
     
    Tensor<std::complex<double> , 1> aa, bb;
    aa = convertTensorUnderlyingType<std::complex<double>, double, 1>(foo);
    bb = convertTensorUnderlyingType<std::complex<double>, double, 1>(bar);


	  bool isEqual = bazMat == foobar(foo,bar,1.);

	  if(!!!isEqual) 
	    std::cout << "you messed this up " << std::endl;

	  TESTER_TEST(m_test,isEqual,"PiPi ff is broken");
	}

      PiPi barbaz;
      TESTER_TEST(m_test,barbaz.nFacs() == 1,"nfacs() is broken");

      return m_test;
    }

// test ffKinematicFactors_t

  tester test_ffKinematicFactor(void)
{
  tester m_test(__func__);
  
  ffKinematicFactors_t<std::complex<double> >
   foo(FormFactorDecompositionFactoryEnv::callFactory(std::string("PiPi")));

  TESTER_TEST(m_test,foo.nFacs() == 1, "nFacs() broken");


  return m_test;
}


  } // close PiPi

} // close radmat
