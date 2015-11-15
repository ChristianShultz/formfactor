// test_factories.cc -
//
// Monday, June  4 2012
//

#include "../headers/test_common.h"
#include "../headers/tester.h"
#include "radmat/utils/pow2assert.h"
#include "radmat/llsq/llsq_solvers.h"
#include "adat/handle.h"
#include <string>


namespace radmat
{

 
  typedef std::complex<double> cd;

  // just check that we can get something
  tester test_solver_factory(void)
  {
    tester m_test(__func__);

    LLSQBaseSolver_t<cd> *ptr;
    // call the factory and get a handle -- the call Factory method takes care of registration the first time its called
    ADAT::Handle<LLSQBaseSolver_t<cd> > test_handle = LLSQSolverFactoryEnv::callFactory(std::string("LU"));

    // get the address of the dereferenced pointer inside the handle then check to make sure it isn't null
    ptr = &*test_handle;    
    TESTER_TEST(m_test,ptr,"factory call returning handle failed");

    // call the factory, get a ptr, check that it isnt null then delete the ptr
    LLSQBaseSolver_t<cd> *ptr2 = TheLLSQSolverFactory::Instance().createObject( std::string("LU") );
    TESTER_TEST(m_test,ptr2,"factory call returning new failed");
    delete ptr2;

    return m_test;
  }

  tester test_mat_elem_factory(void)
  {
    tester m_test(__func__);

    ffBase_t<cd> *ptr;
    ADAT::Handle<ffBase_t<cd> > test_handle = 
      FormFactorDecompositionFactoryEnv::callFactory(std::string("PiPi") );
    ptr = &*test_handle;
    TESTER_TEST(m_test, ptr,"factory call returning handle failed");

    ffBase_t<cd> *ptr2 = 
      TheFormFactorDecompositionFactory::Instance().createObject( std::string("PiPi") );
    TESTER_TEST(m_test,ptr2,"factory call returning new failed");
    delete ptr2;

    return m_test;
  }



}
