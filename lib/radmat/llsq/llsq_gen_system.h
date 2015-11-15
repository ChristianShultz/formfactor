#ifndef LLSQ_GEN_SYSTEM_H_H_GUARD
#define LLSQ_GEN_SYSTEM_H_H_GUARD

#include "radmat/ff_interface/formfactor_factory.h"
#include "radmat/ff_interface/formfactor_kinematic_factors.h"
#include "io/adat_xmlio.h"
#include "radmat/utils/splash.h"
#include "semble/semble_vector.h"
#include "semble/semble_matrix.h"
#include "semble/semble_algebra.h"
#include "semble/semble_file_management.h"
#include "radmat/utils/handle.h"
#include <vector>
#include <string>
#include <sstream>
#include <iostream>

using namespace ADATXML;
using namespace ADATIO;

namespace radmat
{

  struct LLSQDataPoint
  {
    // need to provide the ensems with default starting values so that they can be used
    // in the stl containers b/c the ensem copy constructor doesnt allow one to 
    // copy around empty ensembles

    LLSQDataPoint(void)
    {
      ENSEM::EnsemComplex Ezero;
      Ezero.resize(1); 
      Ezero = SEMBLE::toScalar(std::complex<double>(0.,0.));

      zero.first = false;
      zero.second = Ezero;
      one = zero;
      two = one;
      three = two;

      E_f = ENSEM::real(Ezero);
      E_i = E_f; 
    }

    LLSQDataPoint(const LLSQDataPoint &o)
      : matElemID(o.matElemID) , zero(o.zero) , one(o.one) , two(o.two) , three(o.three),
      p_f(o.p_f) , p_i(o.p_i) , E_f(o.E_f) , E_i(o.E_i) , mom_fac(o.mom_fac)
    {  }

    ENSEM::EnsemReal Q2(void) const
    {
      double pp(0); 
      pp = mom_fac*mom_fac*((p_f[0] - p_i[0])*(p_f[0] - p_i[0])
          + (p_f[1] - p_i[1])*(p_f[1] - p_i[1])
          + (p_f[2] - p_i[2])*(p_f[2] - p_i[2]));


      return ( - (E_f-E_i)*(E_f-E_i) + SEMBLE::toScalar(pp));
    }

    std::string toString(void) const
    {
      std::stringstream ss;

      ss << matElemID << "_pf" << p_f[0] << p_f[1] << p_f[2] << "_pi" << p_i[0] << p_i[1] << p_i[2]; 
      ss << "_t" << zero.first << "_" << SEMBLE::toScalar(ENSEM::mean(zero.second));
      ss << "_x" << one.first << "_" << SEMBLE::toScalar(ENSEM::mean(one.second));
      ss << "_y" << two.first << "_" << SEMBLE::toScalar(ENSEM::mean(two.second));
      ss << "_z" << three.first << "_" << SEMBLE::toScalar(ENSEM::mean(three.second));

      return ss.str(); 
    }


    std::string matElemID;                      
    std::pair<bool,ENSEM::EnsemComplex> zero;   // lorentz index of measurements
    std::pair<bool,ENSEM::EnsemComplex> one;    // bool is if we want to use it
    std::pair<bool,ENSEM::EnsemComplex> two;
    std::pair<bool,ENSEM::EnsemComplex> three;
    Array<int> p_f;                              // the momenta used
    Array<int> p_i;
    ENSEM::EnsemReal E_f;                              // these are in flight energies 
    ENSEM::EnsemReal E_i;                              //  ie: sqrt(m*m + p*p)
    int h_f;
    int h_i; 
    double mom_fac;                           // (1/ xi) * 2pi/L_s -- the "unit" size  

  };




  /////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////

  // deriving from this polymorphicall will inevitably 
  // give us some freedom later, return handles from the solver routines
  template<typename T>
    struct LLSQRetTypeBase_t
    {
      typedef typename SEMBLE::SembleVector<T> FFType; 

      LLSQRetTypeBase_t(const FFType &FFSoln)
        : m_FF(FFSoln)
      {  }

      LLSQRetTypeBase_t(void) 
        : m_FF(SEMBLE::SembleVector<T>())
      {  }

      virtual ~LLSQRetTypeBase_t(void) {}

      virtual std::string echo(void) const {return std::string("LLSQRetTypeBase_t");}

      public: 
      FFType m_FF;
    };


  ///////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////
  template<typename T>
    struct LLSQInputType_t
    {
      typedef typename SEMBLE::SembleVector<T> LatticeMatrixElements;
      typedef typename SEMBLE::SembleMatrix<T> KinematicFactors;

      LLSQInputType_t(const KinematicFactors &KFacs, const LatticeMatrixElements &MatElems)
        : m_KFacs(KFacs) , m_MatElems(MatElems) , qsq(-1000000.)
      {  }


      LLSQInputType_t(const KinematicFactors &KFacs, const LatticeMatrixElements &MatElems,
          const double q2)
        : m_KFacs(KFacs) , m_MatElems(MatElems) , qsq(q2)
      {  }


      virtual std::string echo(void) const {return std::string("LLSQInputType_t");}

      private:
      LLSQInputType_t(void);

      public:
      KinematicFactors m_KFacs;
      LatticeMatrixElements m_MatElems;
      double qsq; 
    };


  // put the guys with active data to the left
  std::vector<LLSQDataPoint> left_sort(const std::vector<LLSQDataPoint> &indata); 


  ///////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////
//
//  template<typename T>
//    rHandle<LLSQInputType_t<T> >
//    generateLLSQSystem(const std::vector<LLSQDataPoint> &unsorted_data, 
//        const int t_ins)
//    {
//
//      typename LLSQInputType_t<T>::KinematicFactors K,KWork;
//      std::vector<ENSEM::EnsemComplex> vectorData;
//      std::vector<LLSQDataPoint>::const_iterator it;
//      std::vector<LLSQDataPoint> data = left_sort(unsorted_data); 
//      bool initK = false;
//
//      it = data.begin();
//      FFKinematicFactors_t genK(FormFactorDecompositionFactoryEnv::callFactory(it->matElemID));
//
//      do {
//        KWork = genK.genFactors(makeMomInvariants(it->E_f,it->E_i,it->p_f,it->p_i,it->mom_fac));
//
//        if(it->zero.first)
//        {
//          initK = true;
//          K.reDim(KWork.getB(),1,KWork.getM());
//          for(int i = 0; i < KWork.getM(); i++)
//            K.loadEnsemElement(0,i,KWork(0,i));
//          vectorData.push_back(it->zero.second);
//
//          if(it->one.first)
//          {
//            K.append_row(KWork.getRow(1));
//            vectorData.push_back(it->one.second);
//          }
//          if(it->two.first)
//          {
//            K.append_row(KWork.getRow(2));
//            vectorData.push_back(it->two.second);
//          }
//          if(it->three.first)
//          {
//            K.append_row(KWork.getRow(3));
//            vectorData.push_back(it->three.second);
//          }
//        }
//        else if(it->one.first)
//        {
//          initK = true;
//          K.reDim(KWork.getB(),1,KWork.getM());
//          for(int i = 0; i < KWork.getM(); i++)
//            K.loadEnsemElement(0,i,KWork(1,i));
//          vectorData.push_back(it->one.second);
//
//          if(it->two.first)
//          {
//            K.append_row(KWork.getRow(2));
//            vectorData.push_back(it->two.second);
//          }
//          if(it->three.first)
//          {
//            K.append_row(KWork.getRow(3));
//            vectorData.push_back(it->three.second);
//          }
//        }
//        else if(it->two.first)
//        {
//          initK = true;
//          K.reDim(KWork.getB(),1,KWork.getM());
//          for(int i = 0; i < KWork.getM(); i++)
//            K.loadEnsemElement(0,i,KWork(2,i));
//          vectorData.push_back(it->two.second);
//
//          if(it->three.first)
//          {
//            K.append_row(KWork.getRow(3));
//            vectorData.push_back(it->three.second);
//          }
//        }
//        else if(it->three.first)
//        {
//          initK = true;
//          K.reDim(KWork.getB(),1,KWork.getM());
//          for(int i = 0; i < KWork.getM(); i++)
//            K.loadEnsemElement(0,i,KWork(3,i));
//          vectorData.push_back(it->three.second);
//        }
//        else
//        {
//          SPLASH("there was a data error in this context, basically christian is a moron");
//          exit(1);
//        }
//      } while(false); // execute once to force a scoped loop --> should be rewritten to break on initK? probably
//
//      do
//      {
//        it ++;
//        if(it == data.end())
//          break;
//
//        ffKinematicFactors_t<T> genK(FormFactorDecompositionFactoryEnv::callFactory(it->matElemID));
//        KWork = genK.genFactors(makeMomInvariants(it->E_f,it->E_i,it->p_f,it->p_i,it->mom_fac));
//
//        if(it->zero.first)
//        {
//          K.append_row(KWork.getRow(0));
//          vectorData.push_back(it->zero.second);
//        }
//        if(it->one.first)
//        {
//          K.append_row(KWork.getRow(1));
//          vectorData.push_back(it->one.second);
//        }
//        if(it->two.first)
//        {
//          K.append_row(KWork.getRow(2));
//          vectorData.push_back(it->two.second);
//        }
//        if(it->three.first)
//        {
//          K.append_row(KWork.getRow(3));
//          vectorData.push_back(it->three.second);
//        }
//
//      } while(true);
//
//      SEMBLE::SembleVector<std::complex<double> >
//        matElems(KWork.getB(),vectorData.size());
//
//      for(unsigned int i =0; i < vectorData.size(); i++)
//        matElems.loadEnsemElement(i,vectorData[i]);
//
//
//      double q2 = SEMBLE::toScalar(ENSEM::mean(data.begin()->Q2())); 
//
//      return rHandle<LLSQInputType_t<T> >( new LLSQInputType_t<T>(K,matElems,q2));
//    }

  ///////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////
  template<typename T>
    struct LLSQBaseSolver_t
    {
      typedef rHandle< LLSQRetTypeBase_t<T> > LLSQRetTypeBase_h;
      typedef rHandle< LLSQInputType_t<T> > LLSQInputType_h;
      virtual void 
        print_llsq_system(const LLSQInputType_h &d,
            const int t_ins) const
        {
          std::string path = SEMBLE::SEMBLEIO::getPath(); 
          path += std::string("llsq");
          SEMBLE::SEMBLEIO::makeDirectoryPath(path);
          path += std::string("/system");
          SEMBLE::SEMBLEIO::makeDirectoryPath(path);  
          std::stringstream Q2;
          Q2 << "/Q2_" << d->qsq; 
          path += Q2.str(); 
          SEMBLE::SEMBLEIO::makeDirectoryPath(path); 
          std::stringstream ss; 
          ss << path << "/Q2_" << d->qsq << "__t_ins_" << t_ins <<"__"; 
          std::string A = ss.str() + std::string("A.txt");
          std::string b = ss.str() + std::string("b.txt"); 

          std::ofstream AA, bb;
          AA.open(A.c_str()); 
          AA << d->m_KFacs.mean(); 
          AA.close();
          bb.open(b.c_str()); 
          bb << d->m_MatElems.mean(); 
          bb.close(); 
        }
      virtual void
        print_llsq_soln(const LLSQInputType_h &d,
            const LLSQRetTypeBase_h &soln, 
            const int t_ins) const
        {
          std::string path = SEMBLE::SEMBLEIO::getPath(); 
          path += std::string("llsq");
          SEMBLE::SEMBLEIO::makeDirectoryPath(path); 
          path += std::string("/system");
          SEMBLE::SEMBLEIO::makeDirectoryPath(path);  
          std::stringstream Q2;
          Q2 << "/Q2_" << d->qsq; 
          path += Q2.str(); 
          SEMBLE::SEMBLEIO::makeDirectoryPath(path); 
          std::stringstream ss; 
          ss << path << "/Q2_" << d->qsq << "__t_ins_" << t_ins <<"__"; 

          std::string x = ss.str() + std::string("x.txt");
          std::ofstream xx;
          xx.open(x.c_str());
          xx << soln->m_FF.mean();
          xx.close();
        }
      virtual LLSQRetTypeBase_h 
        operator()(const LLSQInputType_h &, 
            const int t_ins) = 0;

      virtual bool invertable(void) const = 0;  
      virtual SEMBLE::SembleMatrix<T> 
        inv(const SEMBLE::SembleMatrix<T> &in) = 0;

      virtual std::string echo(void) const = 0;
      virtual ~LLSQBaseSolver_t(void) {}

      virtual std::string solution_log(void) const
      { return my_solution_log; }

      virtual void set_solution_log(const std::string &s) 
      { my_solution_log = s; }

      std::string my_solution_log; 
    };



}

#endif
