// test_LLSQ_solver.cc -
//
// Monday, June  4 2012
//

#include "../headers/test_common.h"
#include "../headers/tester.h"
#include "radmat/utils/pow2assert.h"
#include "radmat/llsq/llsq_solvers.h"
#include "radmat/llsq/llsq_gen_system.h"
#include "radmat/llsq/llsq_solver.h"
#include "adat/handle.h"
#include <string>
#include <complex>
#include "semble/semble_semble.h"


namespace radmat
{
 
  typedef std::complex<double> dc;
  typedef ADAT::Handle<LLSQBaseSolver_t<dc> > csolver_h;
  typedef ADAT::Handle<LLSQRetTypeBase_t<dc> > cret_h;


  tester test_LLSQ_solver_SVDMakeSquare(void)
  {
    tester m_test(__func__);
    
    csolver_h svd_h = LLSQSolverFactoryEnv::callFactory(std::string("SVDMakeSquare"));
    TESTER_TEST(m_test,&*svd_h,"factory call broken");

    const int ncfg(500);
    const int nmatElems(4);
    const int nFF(2);
    SEMBLE::SembleMatrix<dc> A(ncfg,nmatElems,nFF);
    SEMBLE::SembleVector<dc> b(ncfg,nmatElems),x;

    itpp::Mat<dc> iA = itpp::randn_c(nmatElems,nFF);
    itpp::Vec<dc> ix = itpp::randn_c(nFF);

    A = iA;		  
    b = iA*ix;

    LLSQInputType_t<dc> *inputPtr = new LLSQInputType_t<dc>(A,b);
    POW2_ASSERT(inputPtr);
    ADAT::Handle<LLSQInputType_t<dc> > inputHandle(inputPtr);
   
    cret_h solnHandle = (*svd_h)(inputHandle,4);
    x = solnHandle->m_FF;
   
    SEMBLE::SembleMatrix<dc> U,V;
    SEMBLE::SembleVector<double> s;
    std::string log = SEMBLE::svd(SEMBLE::adj(A)*A,U,s,V);
    SEMBLE::SembleVector<dc> xx;
    SEMBLE::pseudoInvert(s,s.getN(),true);
    xx = V * SEMBLE::diagAsym<dc,double>(s) * SEMBLE::adj(U) *SEMBLE::adj(A) * b;

    /*
      std::cout << SEMBLE::mean(xx) << std::endl;
      std::cout << SEMBLE::mean(SEMBLE::adj(A) * b) << std::endl;
      std::cout << "x:\n" << ix << std::endl;
      std::cout << "\n<x>:\n" << SEMBLE::mean(x) << std::endl;
    */

    TESTER_TEST(m_test,x.getB() == ncfg,"soln ncfg broken");
    TESTER_TEST(m_test,x.getN() == nFF,"soln nFF broken");
    
    itpp::Vec<dc> zero(nFF);
    zero.zeros();
    bool solverLLSQ = true;
    for(int i = 0; i < ncfg; i ++)
      	if( (xx[i] - x[i] != zero)
	    || (itpp::round_to_zero(x[i] - ix,1e-10) != zero)) 
	  {
	    solverLLSQ = false;
	    break;
	  }

      
    
    TESTER_TEST(m_test,solverLLSQ,"couldn't solve a linear system correctly");

    return m_test;
  }



} // close namespace
