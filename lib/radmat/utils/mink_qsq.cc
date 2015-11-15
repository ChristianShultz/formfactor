/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : mink_qsq.cc

 * Purpose :

 * Creation Date : 01-11-2013

 * Last Modified : Fri 21 Feb 2014 09:58:21 AM EST

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/

#include "mink_qsq.h"
#include <math.h>

namespace radmat
{

  double Mink_qsq(const ADATXML::Array<int> &pf, const double mf, 
      const ADATXML::Array<int> &pi, const double mi,
      const double factor)
  {
    ADATXML::Array<double> p_f,p_i; 
    p_f.resize(3);
    p_i.resize(3); 


    double p_i_sq(0.);
    double p_f_sq(0.); 
    double q_space_sq(0.);

    for(int i = 0; i < 3; ++i)
    {        
      p_f[i] = factor*double(pf[i]);
      p_i[i] = factor*double(pi[i]); 
      p_f_sq += p_f[i]*p_f[i];
      p_i_sq += p_i[i]*p_i[i];
      q_space_sq += (p_i[i] - p_f[i])*(p_i[i] - p_f[i]);
    }


    double q_time =  sqrt(mi*mi + p_i_sq) - sqrt(mf*mf + p_f_sq);  

    return -((q_time*q_time) - q_space_sq); 
  }


  double mom_factor(const double xi, const int L_s) 
  {
    return 2.*acos(-1.)/xi/double(L_s);
  }


} // radmat
