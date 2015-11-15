/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : test_cov_vecs.cc

 * Purpose :

 * Creation Date : 16-07-2012

 * Last Modified : Tue Jul 17 14:57:51 2012

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/



#include "../headers/test_common.h"
#include "../headers/tester.h"
#include "radmat/fake_data/covarrying_vectors.h"
#include "itpp/itbase.h"
#include <iostream>
#include <complex>
#include <vector>
#include <stdlib.h>
#include <time.h>
#include "radmat/fake_data/fake_overlaps.h"
#include "radmat/fake_data/fake_data_ini.h"
#include "semble/semble_matrix.h"



using namespace radmat;


  int 
main(void)
{
  
     srand(time(NULL));
     itpp::RNG_reset(rand()); // reseed itpp rng with time


     FakeDataIni_t fakeIni;
     fakeIni.stateProps.zProps.overlapGenerator = std::string("unitary");
     fakeIni.dataProps.ncfg = 50;
     fakeIni.stateProps.zProps.updateCovariance = true;
     fakeIni.stateProps.zProps.updateVariance = true;
     fakeIni.stateProps.zProps.varianceOrder = 0.1;
     fakeIni.timeProps.tsink = 5;
     fakeIni.timeProps.tsource = 0;

     FakeOverlaps fakeLaps(fakeIni);

     std::vector<SEMBLE::SembleMatrix<double> > doublefoo = fakeLaps.generate<double>(5,std::string("source"));
     std::vector<SEMBLE::SembleMatrix<double> >::const_iterator dit;

     for(dit = doublefoo.begin(); dit != doublefoo.end(); dit++)
     std::cout << dit->mean() << std::endl;

   /*

  int n = 5;
  corMat genCorrelationMatrix;
  itpp::Mat<std::complex<double> > C(genCorrelationMatrix.genMat<std::complex<double> >(n));

  std::cout << "correlation matrix in \n" << C << std::endl; 

  wrapRandN rng;
  itpp::Vec<std::complex<double> > mean(rng.wrap<std::complex<double> >(n));
  itpp::Vec<double> var(rng.wrap<double>(n));

  std::cout << "mean in \n" << mean << std::endl;
std::cout << "var in \n" << var << std::endl;

*/
  return 0;
}





