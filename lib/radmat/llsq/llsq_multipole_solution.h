#ifndef LLSQ_MULTIPOLE_SOLUTION_H
#define LLSQ_MULTIPOLE_SOLUTION_H 

#include "llsq_solution.h"
#include "radmat/ff_interface/ff_interface.h"

template<typename T> 
  FormFacSolutions<T> 
convert_to_multipole(const FormFacSolutions &arb,
    const std::string &matrix_element_id )
{
  // pull the data out of the solution
  SEMBLE::SembleMatrix<T> ff = arb.FF_t; 
  std::vector<ThreePointDataTag> ingredients = arb.Ingredients; 
  std::vector<std::string> names = arb.names; 

  // get the decomposition
  rHandle<FormFactorBase_t> decomp;
  decomp = FormFactorDecompositionFactoryEnv::callFactory(matrix_element_id); 

}


// args are 1*spin
double Wigner3J(const int J1, const int J2, const int J3 , 
                const int m1, const int m2, const int m3);

#endif /* LLSQ_MULTIPOLE_SOLUTION_H */
