// test_breit.cc -
//
// Thursday, April 26 2012
//

#include "radmat/utils/breit_frame.h"
#include "itpp/itbase.h"
#include <iostream> 

using namespace radmat::breit; 
using namespace std;
using namespace itpp; 


void set_mass(Vec<double> &p, const double m)
{
  p[0] = sqrt(m*m + p[1]*p[1] + p[2]*p[2] + p[3]*p[3]);
}

double get_mass(const Vec<double> pmu)
{
  return  sqrt(pmu * (gmunu() *pmu) ); 
}

int 
main(void)
{
  double tol = 1e-14;
  itpp::Vec<double> pmu(4) , pmup(4); 

  pmu = randu(4);
  pmup = randu(4);



  pmu[1] = 0.397959;
  pmu[2] = 0.618845;
  pmu[3] = 0.837407;
  pmup[1] = 0.732995;
  pmup[2] = 0.922663;
  pmup[3] = 0.60049; 

  set_mass(pmu, 0.5); 
  set_mass(pmup, 0.2); 

  // cout << "pA " << pmu << endl;
  // cout << "pB " << pmup << endl;

  LorentzTransform lt = genBreitLT(pmup,pmu); 

  /*
  std::cout << "G \n" << lt.lambda << std::endl;
  std::cout << "det(G) = " << det(lt.lambda) << std::endl;

  std::cout << "pmu_breit " << lt.lambda * pmu << std::endl; 
  std::cout << "pmup_breit " << lt.lambda * pmup << std::endl; 
  std::cout << get_mass(lt.lambda*pmu) << std::endl; 
  std::cout << get_mass(lt.lambda*pmup) << std::endl; 
  */
  LorentzTransform LT = genRotMat(lt,pmup,pmu); 

  /*
  std::cout << " R_G \n" << round_to_zero(LT.lambda,tol) << std::endl; 
  std::cout << "R_G_inv \n" << round_to_zero(LT.lambda_inv,tol) << std::endl;


  std::cout << round_to_zero(LT.lambda * transpose(LT.lambda),tol) << std::endl; 
  std::cout << round_to_zero(LT.lambda_inv * transpose(LT.lambda_inv),tol) << std::endl;
  */

  return 0;

}

