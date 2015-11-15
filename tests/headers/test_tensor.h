#ifndef TEST_TENSOR_H_H_GUARD
#define TEST_TENSOR_H_H_GUARD

#include "tester.h"
#include <algorithm>
#include <functional>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "radmat/utils/tensor.h"
#include "radmat/utils/pow2assert.h"


namespace radmat
{
  template<typename T>
  tester test_tensor(const T &)
  {
    tester m_test(__PRETTY_FUNCTION__);

    T factor1 = T(0);
    Tensor<T,4> test_resize;

    bool resizeable;

    resizeable = test_resize.resize((TensorShape<4>())[4][4][4][4], factor1);
    TESTER_TEST(m_test,resizeable,"resize failed");

    TensorShape<0> foo = TensorShape<4>()[4][4][4][4];
    Tensor<T,4> test_shape1((TensorShape<4>())[4][4][4][4]);
    Tensor<T,4> test_shape2(foo);

    T factor2 = T(10);
    Tensor<T,4> test_init((TensorShape<4>())[4][4][4][4], factor2);
    int init = std::count_if(test_init.begin(),test_init.end(),
			     std::bind1st(std::equal_to<T>(),factor2));
    TESTER_TEST(m_test,init == 4*4*4*4,"init const failed");

    Tensor<T,4> test_assignment(test_init);
    bool assignment = std::equal(test_init.begin(), test_init.end(), test_assignment.begin());
    TESTER_TEST(m_test,assignment,"assignment operator failed");

    TESTER_TEST(m_test,test_init == test_assignment,"equal failed");
  
    bool _equals = equals<T>(test_assignment, test_init, T(1e-14));
    TESTER_TEST(m_test,_equals,"approx equals failed");
  
    TESTER_TEST(m_test,test_resize != test_init, "!= failed");

    bool idx_correct = (std::vector<bool>(4,true) == test_assignment.getRaised());
    TESTER_TEST(m_test,idx_correct,"initial index pos failed");

    Tensor<T,4> test_empty;
    TESTER_TEST(m_test,test_empty.empty(),"empty failed");
    
    resizeable = test_empty.resize((TensorShape<4>())[4][4][4][4],0);
    TESTER_TEST(m_test,!test_empty.empty(),"empty failed");

    T factor = T(2.71828183);
    Tensor<T,4> test_algebra_pluseq((TensorShape<4>())[4][4][4][4],0);
    test_algebra_pluseq += factor;
    bool test_op_pluseq = test_algebra_pluseq.end() == std::find_if(test_algebra_pluseq.begin(), 
								    test_algebra_pluseq.end(),
								    std::bind1st(std::not_equal_to<T>(),factor));
    TESTER_TEST(m_test,test_op_pluseq,"+= factor failed");

    Tensor<T,4> test_algebra_minuseq(test_algebra_pluseq);
    test_algebra_minuseq -= factor;
    bool test_op_minuseq = test_algebra_minuseq.end() == std::find_if(test_algebra_minuseq.begin(),
								      test_algebra_minuseq.end(),
								      std::bind1st(std::not_equal_to<T>(),0));
    TESTER_TEST(m_test,test_op_minuseq,"-= factor failed");

    Tensor<T,4> test_algebra_teq((TensorShape<4>())[4][4][4][4],1);
    test_algebra_teq *= factor;
    bool test_op_teq = test_algebra_teq.end() == std::find_if(test_algebra_teq.begin(), test_algebra_teq.end(),
							      std::bind1st(std::not_equal_to<T>(),factor));
    TESTER_TEST(m_test,test_op_teq,"*= factor failed");
    
    test_algebra_teq /= factor;
    bool test_op_deq = test_algebra_teq.end() == std::find_if(test_algebra_teq.begin(), test_algebra_teq.end(),
							      std::bind1st(std::not_equal_to<T>(),1));
    TESTER_TEST(m_test,test_op_deq,"/= factor failed");

    Tensor<T,4> test_algebra_plus_eq_t((TensorShape<4>())[4][4][4][4],0);
    test_algebra_plus_eq_t += test_algebra_teq;
    bool test_op_peq_t = std::equal(test_algebra_plus_eq_t.begin(), test_algebra_plus_eq_t.end(), test_algebra_teq.begin());
    TESTER_TEST(m_test,test_op_peq_t,"+= tensor failed");

    Tensor<T,4> test_algebra_plus = test_algebra_pluseq + test_algebra_pluseq;
    bool test_op_plus = test_algebra_plus.end() == std::find_if(test_algebra_plus.begin(),
								test_algebra_plus.end(),
								std::bind1st(std::not_equal_to<T>(),2.*factor));
    TESTER_TEST(m_test,test_op_plus,"operator + failed");

    test_algebra_plus -= test_algebra_pluseq;
    bool test_op_meq = std::equal(test_algebra_plus.begin(), test_algebra_plus.end(), test_algebra_pluseq.begin());
    TESTER_TEST(m_test,test_op_meq,"-= tensor failed");

    test_algebra_plus = test_algebra_pluseq - test_algebra_pluseq;
    bool test_op_m = test_algebra_plus.end() == std::find_if(test_algebra_plus.begin(), test_algebra_plus.end(),
						       std::bind1st(std::not_equal_to<T>(),0));
    TESTER_TEST(m_test,test_op_m,"operator minus failed");
   
    idx_t bound4 = 4;
    Tensor<T,4> test_access_copy((TensorShape<4>())[bound4][bound4][bound4][bound4],0.);
    
    for(idx_t i = 0; i < bound4; ++i)
      for(idx_t j = 0; j < bound4; ++j)
	for(idx_t k = 0; k < bound4; ++k)
	  for(idx_t l = 0; l < bound4; ++l)
	    test_access_copy[i][j][k][l] = i*(j + 5*k) - sqrt(7*l) * l + 1./(double(k + i + 1)) + pow(2,j);
    
    Tensor<T,4> test_access_copy2 = test_access_copy;
    bool test_access_copy_success = true;
    for(idx_t i = 0; i < bound4; ++i)
      for(idx_t j = 0; j < bound4; ++j)
	for(idx_t k = 0; k < bound4; ++k)
	  for(idx_t l = 0; l < bound4; ++l)
	    if(test_access_copy2[i][j][k][l] != i*(j + 5*k) - sqrt(7*l) * l + 1./(double(k + i + 1)) + pow(2,j))
	      {
		test_access_copy_success = false;
		break;
	      }

    TESTER_TEST(m_test,test_access_copy_success,"access copy success failed");
 
    Tensor<T,4> contractA((TensorShape<4>())[4][4][4][4], 0.);
    Tensor<T,2> contract_mI4((TensorShape<2>())[4][4],0.);
    Tensor<T,4> contractB((TensorShape<4>())[4][4][4][4],0.);

    contract_mI4.setRaised(std::vector<bool>(2,false));
    contractB.lower_index(3);

    for(idx_t i = 0; i < 4; ++i)
      { 
	contract_mI4[i][i] = -1.;
	for(idx_t j = 0; j < 4; ++j)
	  for(idx_t k = 0; k < 4; ++k)
	    for(idx_t l = 0; l < 4; ++l)
	      {
		contractA[i][j][k][l] = i*(j + 5*k) - sqrt(7*l) * l + 1./(double(k + i + 1)) + pow(2,j);
		contractB[i][j][k][l] = contractA[i][j][k][l] == 0. ? 0 : 1./contractA[i][j][k][l];
	      }
      }
    
    //T^{(4)}_{abcd}T^{(2)}_{ea} = T^{(4)}_{bcde}
    Tensor<T,4> ctrackA0_I4_1 = contract(contractA,contract_mI4,0,1);
    Tensor<T,4> ctrackA0_I4_0 = contract(contractA,contract_mI4,0,0);

    //T^{(4)}_{abcd}T^{(2)}_{be} = T^{(4)}_{acde}
    Tensor<T,4> ctrackB1_I4_0 = contract(contractB,contract_mI4,1,0);
    Tensor<T,4> ctrackB1_I4_1 = contract(contractB,contract_mI4,1,1);

    //T^{(4)}_{abcd}T^{(4)}_{efgc} = T^{(6)}_{abdefg}
    Tensor<T,6> ctrackA2_B3 = contract(contractA,contractB,2,3);

    // contract_mI4 is symmetric
    TESTER_TEST(m_test,ctrackA0_I4_1 == ctrackA0_I4_0,"contraction with minus identitiy failed");
    TESTER_TEST(m_test,ctrackB1_I4_1 == ctrackB1_I4_1,"contraction with minus identitiy failed");
    
    bool AcB(true), AcI(true), BcI(true);

    Tensor<T,6> _AcB((TensorShape<6>())[4][4][4][4][4][4],0.);
    Tensor<T,4> _AcI((TensorShape<4>())[4][4][4][4],0.), _BcI((TensorShape<4>())[4][4][4][4],0.);

    // sigh -- this is why we need contract to work, these large sum loops are a pain to write
    for(idx_t i = 0; i < 4; ++i)
      for(idx_t j = 0; j < 4; ++j)
	for(idx_t k = 0; k < 4; ++k)
	  for(idx_t l = 0; l < 4; ++l)
	    for(idx_t ii = 0; ii < 4; ++ii)
		{
		  _AcI[j][k][l][ii] += contractA[i][j][k][l]*contract_mI4[ii][i];
		  _BcI[i][k][l][ii] += contractB[i][j][k][l]*contract_mI4[j][ii];

		  for(idx_t jj = 0; jj < 4; ++jj)
		    for(idx_t kk = 0; kk < 4; ++kk)
		     _AcB[i][j][l][ii][jj][kk] += contractA[i][j][k][l] * contractB[ii][jj][kk][k];
		}

    // indicies need to be kept track of.. contract does this for us.
    _AcI.lower_index(3);
    _BcI.lower_index(3);
    _BcI.lower_index(2);

    // this should be strictly equal, not any of that approximate nonsense
    AcI = (_AcI == ctrackA0_I4_0);
    BcI = (_BcI == ctrackB1_I4_1);


    TESTER_TEST(m_test,AcI,"contraction failed AcI");
    TESTER_TEST(m_test,BcI,"contraction failed BcI");

    AcB = equals(ctrackA2_B3,_AcB,T(1e-14));

    TESTER_TEST(m_test,AcB,"contraction failed AcB");

    //! now test the clear function
    contractA.clear();
    contractB.clear();
    _AcB.clear();

    TESTER_TEST(m_test,contractA.empty(),"clear failed");

    contractA.resize((TensorShape<4>())[4][4][4][4],0.);
    contractB = contractA;
    _AcB.resize((TensorShape<6>())[4][4][4][4][4][4],0.);

    srand( time(NULL) );
    
    // uniform rand in (-1,1)
    for(idx_t i = 0; i < 4; ++i)
      for(idx_t j = 0; j < 4; ++j)
	for(idx_t k = 0; k < 4; ++k)
	  for(idx_t l = 0; l < 4; ++l)
	    {
	      contractA[i][j][k][l] = 2.*( (double)rand()/((double)RAND_MAX) - 0.5);
	      contractB[i][j][k][l] = 2.*( (double)rand()/((double)RAND_MAX) - 0.5);
	    }

    contractB.lower_index(0);
    _AcB = contract(contractA,contractB,3,0);

    Tensor<T,6> __AcB((TensorShape<6>())[4][4][4][4][4][4],0.);
    
    for(idx_t i = 0; i < 4; ++i)
      for(idx_t j = 0; j < 4; ++j)
	for(idx_t k = 0; k < 4; ++k)
	  for(idx_t l = 0; l < 4; ++l)
	    for(idx_t ii = 0; ii < 4; ++ii)
	      for(idx_t jj = 0; jj < 4; ++jj)
		for(idx_t kk = 0; kk < 4; ++kk)
		  __AcB[i][j][k][ii][jj][kk] += contractA[i][j][k][l] * contractB[l][ii][jj][kk];

    TESTER_TEST(m_test,__AcB == _AcB, "contraction failed");


    return m_test;
  }

}

#endif
