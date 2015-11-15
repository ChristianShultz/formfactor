// harmonic.cc -
//
// Wednesday, June 20 2012
//

#include"harmonic.h"
#include <math.h>


namespace radmat
{

 // force a lookup table to be computed at compile time..
  int FactorialTable(const int n)
  {
    switch(n)
      {
	// these get overflow..
	// case 20: return Factorial<20>::value;
	// case 19: return Factorial<19>::value;
	// case 18: return Factorial<18>::value;
	// case 17: return Factorial<17>::value;
	// case 16: return Factorial<16>::value;
	// case 15: return Factorial<15>::value;
	// case 14: return Factorial<14>::value;
	// case 13: return Factorial<13>::value;
      case 12: return Factorial<12>::value;
      case 11: return Factorial<11>::value;
      case 10: return Factorial<10>::value;
      case 9: return Factorial<9>::value;
      case 8: return Factorial<8>::value;
      case 7: return Factorial<7>::value;
      case 6: return Factorial<6>::value;
      case 5: return Factorial<5>::value;
      case 4: return Factorial<4>::value;
      case 3: return Factorial<3>::value;
      case 2: return Factorial<2>::value;
      case 1: return Factorial<1>::value;
      case 0: return Factorial<0>::value;
      default: return 0;
      }
  }
 
  double HermiteP(const int n, const double x)
  {
    if(n == 0)
      return 1.;
    if(n == 1)
      return 2.*x;
    return 2.*(x*HermiteP(n-1, x) - double(n -1)*HermiteP(n-2,x));    
  }

  namespace
  {
    const double pi = 2.*acos(0);
  }

  double HarmonicNormalization(const int n)
  {
    return 1./(
	       sqrt(
		    double(FactorialTable(n))
		    *pow(2,n)
		    *sqrt(pi)
		    )
	       );      
  }

  double Harmonic(const int n, const double x)
  {
    return HermiteP(n,x)*HarmonicNormalization(n)*exp(-x*x/2.);
  }

  double HarmonicPlusOne(const int n, const double x)
  {
    return Harmonic(n,x) + 1;
  }

  double One(const int, const double)
  {
  return double(1);
  }

}
