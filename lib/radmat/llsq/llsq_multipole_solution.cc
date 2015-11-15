/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : llsq_multipole_solution.cc

 * Purpose :

 * Creation Date : 16-07-2014

 * Last Modified : Wed 16 Jul 2014 03:31:20 PM EDT

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/

#include "llsq_multipole_solution.h"
#include "hadron/clebsch.h"

namespace radmat
{

  namespace
  {
    struct non_zero_wigner
    {
      non_zero_wigner(const int hf , const int hg, 
          const int hi, const double val)
        : lf(hf) , lg(hg) , li(hi) , v(val)
      { }

      int lf,lg,li; 
      double v; 
    };

    std::vector<non_zero_wigner> 
      non_zero_wigner_coeffs(const int Jf, const int Ji)
      {
        std::vector<non_zero_wigner> ret; 
        for(int hf = -Jf; hf <= Jf; ++hf)
          for(int hg = -1; hg <= 1; ++hg)
            for(int hi = -Ji; hi <= Ji; ++hi)
            {
              double CG = Wigner3J(Jf,1,Ji,hf,hg,hi); 
              if( fabs(CG) < 1e-6) 
                continue; 

              ret.push_back(non_zero_wigner(hf,hg,hi,CG)); 
            }

        return ret; 
      }


  }

  double Wigner3J(const int J1, const int J2, const int J3 , 
      const int m1, const int m2, const int m3)
  {
    double prefactor = ((J1 -J2 -m3) %2 == 0) ? 1. : -1.; 
    prefactor *= 1./(sqrt(2.*J3 +1)); 

    return prefactor * Hadron::clebsch(2*J1,2*m1,2*J2,2*m2,2*J3,-2*m3); 
  }


}


