#ifndef LLSQ_SOLVER_H_H_GUARD
#define LLSQ_SOLVER_H_H_GUARD

#include "llsq_solvers.h"
#include "llsq_gen_system.h"
#include "radmat/utils/handle.h"
#include <string>
#include <complex>

namespace radmat
{


  template<typename T>
  struct LLSQSolver_t
  {
    // save some typing
    typedef typename LLSQBaseSolver_t<T>::LLSQRetTypeBase_h LLSQRetTypeBase_h;
    typedef typename LLSQBaseSolver_t<T>::LLSQInputType_h LLSQInputType_h;
    typedef rHandle<LLSQBaseSolver_t<T> > LLSQBaseSolver_h; 

    LLSQSolver_t(void) {}

    LLSQRetTypeBase_h operator()(const std::vector<LLSQDataPoint> &data, const std::string &solverID, const int t_ins) const
    {
      LLSQBaseSolver_h solver = LLSQSolverFactoryEnv::callFactory(solverID);
      return (*solver)(generateLLSQSystem<T>(data, t_ins), t_ins);
    }

  };

}

#endif
