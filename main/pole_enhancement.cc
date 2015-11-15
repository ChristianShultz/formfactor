/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

* File Name : pole_enhancement.cc

* Purpose :

* Creation Date : 06-05-2014

* Last Modified : Tue 06 May 2014 06:56:12 PM EDT

* Created By : shultz

_._._._._._._._._._._._._._._._._._._._._.*/

#include "radmat/spectrum/spectrum.h"
#include "radmat/register_all/register_all.h"
#include <map>
#include <string>
#include <sstream>
#include "radmat/redstar_interface/redstar_canonical_rotations.h"
#include "radmat/redstar_interface/redstar_canonical_lattice_rotations.h"
#include "radmat/ff/ff.h"
#include "radmat/utils/levi_civita.h"
#include "hadron/irrep_util.h"
#include "formfac/formfac_qsq.h"


typedef radmat::mom_t mom_t; 


void read_momentum(const int st, mom_t &p, char *argv[])
{
  for(int i = st; i < st+3; ++i)
  {
    std::istringstream val(argv[i]);
    val >> p[i-st];
  }
}


int main(int argc, char *argv[])
{

  if( argc != 9 )
  {
    std::cout << "error usage: test_rotations [" << __func__ << "] "
      << "<left_id> <mom1> <right_id> <mom2> " << std::endl; 
    exit(1);
  }

  radmat::AllFactoryEnv::registerAll(); 
  
  std::string l,r,p; 
  p = "rho";

  { std::istringstream val(argv[1]); val  >> l ; }
  { std::istringstream val(argv[5]); val  >> r ; }


  mom_t m1 = radmat::gen_mom<0,0,0>(); 
  mom_t m2 = radmat::gen_mom<0,0,0>(); 

  read_momentum(2,m1,argv);
  read_momentum(6,m2,argv);


  std::pair<double,double> polish = radmat::TheSpectrumFactoryEnv::pole_enhancement(l,m1,r,m2,p); 

  std::cout << " Q2 for the transition is " << polish.second << std::endl;
  std::cout << " pole enhancement factor is " << polish.first << std::endl;

  return 0; 
}
