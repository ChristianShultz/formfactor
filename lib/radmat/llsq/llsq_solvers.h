#ifndef LLSQSOLVERS_H_H_GUARD
#define LLSQSOLVERS_H_H_GUARD


#include "llsq_gen_system.h"
#include <string>
#include <complex>
#include "adat/singleton.h"
#include "adat/objfactory.h"
#include "radmat/utils/handle.h"
#include "semble/semble_vector.h"
#include "semble/semble_matrix.h"
#include "semble/semble_algebra.h"
#include "semble/semble_linear_algebra.h"

namespace radmat
{

  // factory
  typedef Util::SingletonHolder<
    Util::ObjectFactory<LLSQBaseSolver_t<std::complex<double> >,
    std::string,
    void,
    LLSQBaseSolver_t<std::complex<double> >* (*)(),
    Util::StringFactoryError > >
      TheLLSQSolverFactory;

  namespace LLSQSolverFactoryEnv
  {
    bool registerAll();
    rHandle<LLSQBaseSolver_t<std::complex<double> > > callFactory(const std::string &solnID);
  }


  /////////////////////////////////////////////////////////////////////////////

  // solver functors

  // solve a square system via LU decomposition to get the inverse
  // -- this is junk and is only here as a placeholder for testing
  //    we should be using some kind of fancy SVD on non-square 
  //    overdetermined linear systems so don't use this one!!!
  template<typename T> 
    struct LLSQSolverLU_t : public LLSQBaseSolver_t<T>
  {
    typedef typename LLSQBaseSolver_t<T>::LLSQRetTypeBase_h LLSQRetTypeBase_h;
    typedef typename LLSQBaseSolver_t<T>::LLSQInputType_h LLSQInputType_h;

    LLSQRetTypeBase_h operator()(const LLSQInputType_h & input, const int t_ins) 
    {
      print_llsq_system(input, t_ins);
      LLSQRetTypeBase_t<T> *foo = new LLSQRetTypeBase_t<T>();
      POW2_ASSERT(foo);                                              // check pointer alloc 
      POW2_ASSERT(input->m_KFacs.getN() == input->m_KFacs.getM());   // check LLSQ is not underdetermined -- it would still work
      SEMBLE::SembleMatrix<T> Kinv = this->inv(input->m_KFacs); 
      foo->m_FF = Kinv * (input->m_MatElems);
      LLSQRetTypeBase_h ret(foo);
      print_llsq_soln(input,ret,t_ins);

      return ret;
    }


    bool invertable(void) const {return true;}

    SEMBLE::SembleMatrix<T> inv(const SEMBLE::SembleMatrix<T> &in) 
    {
      SEMBLE::SembleMatrix<T> out;
      SEMBLE::inv(in,out);
      return out;     
    }


    std::string echo(void) const {return std::string("LLSQSolverLU_t");}
  };

  // solve A*x = b for a non_square A via
  // Adag*A*x = Adag*b -- makes a square linear system
  // 
  // TO DO: write a version of this that does the whole nonsense with the residual?
  //        our system shouldn't have null space so maybe this is a waste of time..
  //

  template<typename T> 
    struct LLSQSolverSVDMakeSquare_t : public LLSQBaseSolver_t<T>
  {
    typedef typename LLSQBaseSolver_t<T>::LLSQRetTypeBase_h LLSQRetTypeBase_h;
    typedef typename LLSQBaseSolver_t<T>::LLSQInputType_h LLSQInputType_h;

    LLSQRetTypeBase_h operator()(const LLSQInputType_h &input, const int t_ins) 
    {
      print_llsq_system(input, t_ins);
      LLSQRetTypeBase_t<T> *foo = new LLSQRetTypeBase_t<T>();
      POW2_ASSERT(foo);                                              // check pointer alloc 
      POW2_ASSERT(input->m_KFacs.getN() >= input->m_KFacs.getM());   // check LLSQ is not underdetermined -- it would still work
      SEMBLE::SembleMatrix<T> Kinv = this->inv(input->m_KFacs); 
      foo->m_FF = Kinv * (input->m_MatElems);
      LLSQRetTypeBase_h ret(foo);
      print_llsq_soln(input,ret,t_ins);

      return ret;
    };

    bool invertable(void) const {return true;}

    // Ax = b, A'Ax = A'b , x = (A'A)^-1 * A'b  prime means dagger
    SEMBLE::SembleMatrix<T> inv(const SEMBLE::SembleMatrix<T> &in) 
    {
      SEMBLE::SembleMatrix<T> out;
      SEMBLE::SembleMatrix<T> KinvDag_Kinv = SEMBLE::adj(in)*(in),U,V;
      SEMBLE::SembleVector<double> s;
      set_solution_log( SEMBLE::svd(KinvDag_Kinv,U,s,V) );
      SEMBLE::pseudoInvert(s,s.getN(),true); // s -> 1/s
      out = V * (SEMBLE::diagAsym<T,double>(s) )* (SEMBLE::adj(U) ) * SEMBLE::adj(in) ;   
      return out;     
    }


    std::string echo(void) const {return std::string("LLSQSolverSVDMakeSquare_t");}

  };


  // solve A*x = b for non-square A via svd decomp of a non-square matrix
  // to do -- also return residual of the solution 


  template<typename T>
    struct LLSQSolverSVDNonSquare_t : public LLSQBaseSolver_t<T>
  {

    typedef typename LLSQBaseSolver_t<T>::LLSQRetTypeBase_h LLSQRetTypeBase_h;
    typedef typename LLSQBaseSolver_t<T>::LLSQInputType_h LLSQInputType_h;

    LLSQRetTypeBase_h operator()(const LLSQInputType_h &input, const int t_ins) 
    {
      print_llsq_system(input, t_ins);
      LLSQRetTypeBase_t<T> *foo = new LLSQRetTypeBase_t<T>();
      POW2_ASSERT(foo);                                              // check pointer alloc 
      POW2_ASSERT(input->m_KFacs.getN() >= input->m_KFacs.getM());   // check LLSQ is not underdetermined -- it would still work
      SEMBLE::SembleMatrix<T> Kinv = this->inv(input->m_KFacs);

      foo->m_FF = Kinv *  (input->m_MatElems);

      LLSQRetTypeBase_h ret(foo);
      print_llsq_soln(input,ret,t_ins);

      return ret;

    };

    bool invertable(void) const {return true;}

    SEMBLE::SembleMatrix<T> inv(const SEMBLE::SembleMatrix<T> &in) 
    {
      SEMBLE::SembleMatrix<T> K(in); 
      SEMBLE::SembleMatrix<T> Kinv; 
      SEMBLE::svdNonSquare(K,Kinv); 
      return Kinv; 
    }

    std::string echo(void) const {return std::string("LLSQSolverSVDNonSquare_t");}

  };



  namespace FancyChisqExtremization 
  {
    void doExtremization(const SEMBLE::SembleVector<double> &C, 
        const SEMBLE::SembleMatrix<double> &K, 
        SEMBLE::SembleVector<double> &FF); 
    void doExtremization(const SEMBLE::SembleVector<std::complex<double> > &C, 
        const SEMBLE::SembleMatrix<std::complex<double> > &K ,
        SEMBLE::SembleVector<double> &FF);
  }


  // transform the llsq and return the soln that minimizes chisq 
  template<typename T> 
    struct LLSQExtremizeChisq_t : public LLSQBaseSolver_t<T>
  {
    // technically we will always be returning real types b/c the solution method 
    // explicitly requires it but we can get away free by just returning ensemble zero 
    // for the imaginary part
    typedef typename LLSQBaseSolver_t<T>::LLSQRetTypeBase_h LLSQRetTypeBase_h;
    typedef typename LLSQBaseSolver_t<T>::LLSQInputType_h LLSQInputType_h;


    LLSQRetTypeBase_h operator()(const LLSQInputType_h &input, const int t_ins) 
    {
      print_llsq_system(input, t_ins);

      LLSQRetTypeBase_t<T> *foo = new LLSQRetTypeBase_t<T>();
      POW2_ASSERT(foo);                                              // check pointer alloc 
      POW2_ASSERT(input->m_KFacs.getN() >= input->m_KFacs.getM());   // check LLSQ is not underdetermined -- it would still work

      FancyChisqExtremization::doExtremization(input->m_MatElems,input->m_KFacs,foo->m_FF);
      LLSQRetTypeBase_h ret(foo);
      print_llsq_soln(input,ret,t_ins);

      return ret;

    };

    bool invertable(void) const {return false;}


    SEMBLE::SembleMatrix<T> inv(const SEMBLE::SembleMatrix<T> &in) 
    {
      std::cerr << __func__ << ": error: unsupported operation " << std::endl;
      exit(1);
    }

    std::string echo(void) const {return std::string("LLSQExtremizeChisq_t");}
  };





}

#endif
