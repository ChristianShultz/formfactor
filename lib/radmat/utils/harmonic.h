#ifndef HARMONIC_H_H_GUARD
#define HARMONIC_H_H_GUARD


namespace radmat
{
  
  template <int N>
  struct Factorial {
    enum { value = N * Factorial<N - 1>::value };
  };
 
  template <>
  struct Factorial<0> {
    enum { value = 1 };
  };

 
  double HermiteP(const int n, const double x);
  double Harmonic(const int n, const double x);
  double HarmonicPlusOne(const int n, const double x);
  double One(const int , const double); // returns one
      
}
#endif
