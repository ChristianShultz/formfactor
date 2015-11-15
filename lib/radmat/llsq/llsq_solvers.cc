// LLSQSolvers.cc -
//
// Saturday, June  2 2012
//

#include"llsq_solvers.h"
#include <string>
#include <complex>
#include <exception>
#include "adat/singleton.h"
#include "adat/objfactory.h"
#include "radmat/utils/pow2assert.h"
#include "semble/semble_linear_algebra.h"
#include "radmat/utils/printer.h"

namespace FacEnv = radmat::LLSQSolverFactoryEnv;
typedef radmat::TheLLSQSolverFactory Factory;

namespace radmat
{


  namespace FancyChisqExtremization 
  {



    namespace 
    {

      inline double doConj(const double &d)
      {
        return d; 
      }

      inline std::complex<double> doConj(const std::complex<double> &d)
      {
        return std::complex<double>(d.real(),-d.imag()); 
      }


      template<typename T>
        itpp::Mat<T> doReset(const itpp::Vec<double> &s, const double tol)
        {
          itpp::Mat<T> SS(s.size(),s.size());
          SS.zeros();

          for(int i = 0; i < s.size(); ++i)
          {
            if(s(i) < tol)
              SS(i,i) = T(0.);
            else
              SS(i,i) = T(1./s(i));
          }
          return SS;
        }


      template<typename T>
        itpp::Mat<T> InverseCovariance(const SEMBLE::SembleVector<T> &c, const double tol=0.0000001)
        {
          itpp::Vec<T> cbar = SEMBLE::mean(c);
          SEMBLE::SembleVector<T> v,vdag; 

          v = (c - cbar); 

          const int nb = c.getB();  
          const int nele = c.getN();
          itpp::Mat<T> cov(nele,nele);
          cov.zeros();

          for(int row = 0; row < nele; ++ row)
            for(int col = row; col < nele; ++col)
            {
              T accum(0.); 
              for(int cfg = 0; cfg < nb; ++cfg)
                accum += (v[cfg][row])*doConj(v[cfg][col]);

              accum *= 1./double(nb);

              cov(row,col) = accum;
              if(row != col)
                cov(col,row) = doConj(accum);
            }

          itpp::Mat<T> U,V,S;
          itpp::Vec<double> s; 
          itpp::svd(cov,U,s,V);
          S = doReset<T>(s,tol);

          return V*S*itpp::hermitian_transpose(U);

        }



    } // namespace anonomyous








    void doExtremization(const SEMBLE::SembleVector<double> &C, 
        const SEMBLE::SembleMatrix<double> &K, 
        SEMBLE::SembleVector<double> &FF)
    {
      itpp::Mat<double> Cinv = InverseCovariance(C); 
      SEMBLE::SembleMatrix<double> A = SEMBLE::transpose(K)*Cinv*K; 
      SEMBLE::SembleVector<double> b = SEMBLE::transpose(K)*Cinv*C; 
      SEMBLE::SembleMatrix<double> Ainv;

      SEMBLE::pseudoInvertCond(A,Ainv);
      FF = Ainv * b; 
    } 






    void doExtremization(const SEMBLE::SembleVector<std::complex<double> > &C, 
        const SEMBLE::SembleMatrix<std::complex<double> > &K ,
        SEMBLE::SembleVector<std::complex<double> > &FF)
    {
      itpp::Mat<std::complex<double> > CinvStar,Cinv = InverseCovariance(C); 
      SEMBLE::SembleVector<std::complex<double> > CStar,b,bStar;
      SEMBLE::SembleMatrix<std::complex<double> > KStar, A,AStar; 
      CStar = C;
      CStar.conj();
      KStar = K;
      KStar.conj();
      CinvStar = itpp::transpose(Cinv); // hermitian matrix

      b = SEMBLE::transpose(KStar)*Cinv*C;
      bStar = SEMBLE::transpose(K)*CinvStar*CStar; 
      A = SEMBLE::transpose(K)*CinvStar*KStar;
      AStar = SEMBLE::transpose(KStar)*Cinv*K;


      SEMBLE::SembleVector<std::complex<double> > bb;
      SEMBLE::SembleMatrix<std::complex<double> > AA;
      
      bb = b + bStar;
      AA = A + AStar;


      std::cout << __func__  << " the folowing should be real within roundoff" << std::endl;
      std::cout << " b = " << std::endl;
      std::cout << SEMBLE::mean(bb) << std::endl;
      std::cout << " A = " << std::endl;
      std::cout << SEMBLE::mean(AA) << std::endl;

      
      SEMBLE::SembleMatrix<std::complex<double> > AAinv;
      SEMBLE::pseudoInvertCond(AA,AAinv);
      FF = AAinv * bb;
    }

  } // namespace FancyChisqExtremization








  namespace LLSQSolverFactoryEnv
  {

    namespace
    {
      template<class T, class U>
        T* upCast(void)
        {
          T *t = new U();
          POW2_ASSERT(t);
          return t;
        }

      volatile bool registered = false;
    } // close anonymous 



    bool registerAll(void)
    {
      printer_function<debug_print_reg_all>(__PRETTY_FUNCTION__); 

      bool success = true;

        if(!!!registered)
        {
          success &= Factory::Instance().registerObject(std::string("LU"),
              FacEnv::upCast<LLSQBaseSolver_t<std::complex<double> >, 
              LLSQSolverLU_t<std::complex<double> > > );
          success &= Factory::Instance().registerObject(std::string("SVDMakeSquare"),
              FacEnv::upCast<LLSQBaseSolver_t<std::complex<double> >, 
              LLSQSolverSVDMakeSquare_t<std::complex<double> > > );
          success &= Factory::Instance().registerObject(std::string("SVDNonSquare"),
              FacEnv::upCast<LLSQBaseSolver_t<std::complex<double> >, 
              LLSQSolverSVDNonSquare_t<std::complex<double> > > );

          registered = true;
        }

      if(!!! success)
      {   
        throw std::string("reg error in the LLSQFactoryEnv"); 
      }   


      return success;
    }


    // interface between the factory and the outside world for ease of use
    rHandle<LLSQBaseSolver_t<std::complex<double> > > callFactory(const std::string &solnID)
    {
      LLSQBaseSolver_t<std::complex<double> > *foo;
      foo = NULL;
      try
      {
        foo = Factory::Instance().createObject(solnID);
      }
      catch(...)
      {
        POW2_ASSERT(false);
      }

      POW2_ASSERT(foo);
      return rHandle<LLSQBaseSolver_t<std::complex<double> > >(foo);
    }


  } // Env

} // radmat
