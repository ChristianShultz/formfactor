#ifndef LLSQ_CHISQ_H
#define LLSQ_CHISQ_H


#include "semble/semble_semble.h"
#include <complex>


namespace radmat
{

  std::pair<double,int> ChisqAndDoF(const SEMBLE::SembleVector<std::complex<double> > &data,
      const SEMBLE::SembleVector<std::complex<double> > &thy, double tol=1e-6);

}


 
#endif /* LLSQ_CHISQ_H */
