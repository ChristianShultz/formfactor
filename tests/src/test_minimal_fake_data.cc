// test_minimal_fake_data.cc -
//
// Monday, June  4 2012
//

#include "../headers/test_common.h"
#include "../headers/tester.h"
#include "semble/semble_matrix.h"
#include <string>
#include <complex>
#include <algorithm>
#include "radmat/ff/formfactor_factory.h"
#include "radmat/fake_data/minimal_fake_data.h"

namespace radmat
{


  tester test_minimal_fake_data(const std::string &matElemID)
  {
    tester m_test(__func__);
    
    // test the functor part for a pi pi mat elem before testing the interesting bits.. 
    // need to get this right to be able to do the next part anyways
    genMinimalFakeData_c matElem( std::string("PiPi") );      
    itpp::Mat<std::complex<double> > pplus(4,1);
    Tensor<double,1> p_f((TensorShape<1>())[4],0.) , p_i((TensorShape<1>())[4],0.);

    for(int i = 0; i < 4; i++)
      {
	p_f[i] = itpp::randn();
	p_i[i] = itpp::randn();
	pplus(i,0) = std::complex<double>(p_f[i]+p_i[i] , 0.);
      }


    bool minimalPiPiFunctorOp = ( matElem(p_f,p_i) ) == pplus;
    TESTER_TEST(m_test,minimalPiPiFunctorOp,"minimal functor broken for pipi");

    matElem.swapMatElem(matElemID);
    double precision = 1e-3;
    itpp::Mat<std::complex<double> > shouldBeTheMean = matElem(p_f,p_i);
    itpp::Mat<std::complex<double> > actualMean = (matElem(p_f,p_i,1000,precision,true)).mean();
    
    for(int row = 0; row < shouldBeTheMean.rows(); row++)
      for(int col = 0; col < shouldBeTheMean.cols(); col++)
	{
	  bool approxEqual = (norm(shouldBeTheMean(row,col) - actualMean(row,col)) // ~ (A - B)/max(A,B) < prec
			      /std::max(norm(shouldBeTheMean(row,col)),
					norm(actualMean(row,col)))
			      < precision);
	  TESTER_TEST(m_test,approxEqual,"generation of fake data is failing");
	}
    

    return m_test;
  }

} // close namespace
