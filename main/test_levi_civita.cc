/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : test_levi_civita.cc

 * Purpose :

 * Creation Date : 11-12-2013

 * Last Modified : Thu 12 Dec 2013 04:00:31 PM EST

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/

#include "radmat/utils/tensor.h"
#include "radmat/utils/levi_civita.h"
#include <complex>
#include <iostream>

int nchk(0); 
int nonz(0); 

  template<typename T>
void check_exit( const radmat::Tensor<T,4> &t, int a, int b, int c, int d, T expected)
{
  if ( t[a][b][c][d] != expected )
  {
    std::cout << __func__ << ": error " << a << b << c << d 
      << " expected " << expected << " got " << t[a][b][c][d] << std::endl;
    exit(1); 
  }

  if( expected != T(0) )
  {
    std::cout << a << b << c << d << "->" << expected << std::endl;
    ++nonz;
  }
  ++nchk;
}

  template<typename T> 
void test_levi(const radmat::Tensor<T,4> &l)
{
  T p(1);
  T m(-1); 
  T z(0); 

  //*
  check_exit(l,0,0,0,0,z); 
  check_exit(l,1,0,0,0,z); 
  check_exit(l,2,0,0,0,z); 
  check_exit(l,3,0,0,0,z); 
  check_exit(l,0,1,0,0,z); 
  check_exit(l,1,1,0,0,z); 
  check_exit(l,2,1,0,0,z); 
  check_exit(l,3,1,0,0,z); 
  check_exit(l,0,2,0,0,z); 
  check_exit(l,1,2,0,0,z); 
  check_exit(l,2,2,0,0,z); 
  check_exit(l,3,2,0,0,z); 
  check_exit(l,0,3,0,0,z); 
  check_exit(l,1,3,0,0,z); 
  check_exit(l,2,3,0,0,z); 
  check_exit(l,3,3,0,0,z); 
  //
  check_exit(l,0,0,1,0,z); 
  check_exit(l,1,0,1,0,z); 
  check_exit(l,2,0,1,0,z); 
  check_exit(l,3,0,1,0,z); 
  check_exit(l,0,1,1,0,z); 
  check_exit(l,1,1,1,0,z); 
  check_exit(l,2,1,1,0,z); 
  check_exit(l,3,1,1,0,z); 
  check_exit(l,0,2,1,0,z); 
  check_exit(l,1,2,1,0,z); 
  check_exit(l,2,2,1,0,z); 
  check_exit(l,3,2,1,0,p); //  
  check_exit(l,0,3,1,0,z); 
  check_exit(l,1,3,1,0,z); 
  check_exit(l,2,3,1,0,m); //
  check_exit(l,3,3,1,0,z); 
  //
  check_exit(l,0,0,2,0,z); 
  check_exit(l,1,0,2,0,z); 
  check_exit(l,2,0,2,0,z); 
  check_exit(l,3,0,2,0,z); 
  check_exit(l,0,1,2,0,z); 
  check_exit(l,1,1,2,0,z); 
  check_exit(l,2,1,2,0,z); 
  check_exit(l,3,1,2,0,m); //
  check_exit(l,0,2,2,0,z); 
  check_exit(l,1,2,2,0,z); 
  check_exit(l,2,2,2,0,z); 
  check_exit(l,3,2,2,0,z); 
  check_exit(l,0,3,2,0,z); 
  check_exit(l,1,3,2,0,p); //
  check_exit(l,2,3,2,0,z); 
  check_exit(l,3,3,2,0,z); 
  //
  check_exit(l,0,0,3,0,z); 
  check_exit(l,1,0,3,0,z); 
  check_exit(l,2,0,3,0,z); 
  check_exit(l,3,0,3,0,z); 
  check_exit(l,0,1,3,0,z); 
  check_exit(l,1,1,3,0,z); 
  check_exit(l,2,1,3,0,p); //
  check_exit(l,3,1,3,0,z); 
  check_exit(l,0,2,3,0,z); 
  check_exit(l,1,2,3,0,m); //
  check_exit(l,2,2,3,0,z); 
  check_exit(l,3,2,3,0,z); 
  check_exit(l,0,3,3,0,z); 
  check_exit(l,1,3,3,0,z); 
  check_exit(l,2,3,3,0,z); 
  check_exit(l,3,3,3,0,z); 
  //*
  check_exit(l,0,0,0,1,z); 
  check_exit(l,1,0,0,1,z); 
  check_exit(l,2,0,0,1,z); 
  check_exit(l,3,0,0,1,z); 
  check_exit(l,0,1,0,1,z); 
  check_exit(l,1,1,0,1,z); 
  check_exit(l,2,1,0,1,z); 
  check_exit(l,3,1,0,1,z); 
  check_exit(l,0,2,0,1,z); 
  check_exit(l,1,2,0,1,z); 
  check_exit(l,2,2,0,1,z); 
  check_exit(l,3,2,0,1,m); //
  check_exit(l,0,3,0,1,z); 
  check_exit(l,1,3,0,1,z); 
  check_exit(l,2,3,0,1,p); //
  check_exit(l,3,3,0,1,z); 
  //
  check_exit(l,0,0,1,1,z); 
  check_exit(l,1,0,1,1,z); 
  check_exit(l,2,0,1,1,z); 
  check_exit(l,3,0,1,1,z); 
  check_exit(l,0,1,1,1,z); 
  check_exit(l,1,1,1,1,z); 
  check_exit(l,2,1,1,1,z); 
  check_exit(l,3,1,1,1,z); 
  check_exit(l,0,2,1,1,z); 
  check_exit(l,1,2,1,1,z); 
  check_exit(l,2,2,1,1,z); 
  check_exit(l,3,2,1,1,z); 
  check_exit(l,0,3,1,1,z); 
  check_exit(l,1,3,1,1,z); 
  check_exit(l,2,3,1,1,z); 
  check_exit(l,3,3,1,1,z); 
  //
  check_exit(l,0,0,2,1,z); 
  check_exit(l,1,0,2,1,z); 
  check_exit(l,2,0,2,1,z); 
  check_exit(l,3,0,2,1,p); //
  check_exit(l,0,1,2,1,z); 
  check_exit(l,1,1,2,1,z); 
  check_exit(l,2,1,2,1,z); 
  check_exit(l,3,1,2,1,z); 
  check_exit(l,0,2,2,1,z); 
  check_exit(l,1,2,2,1,z); 
  check_exit(l,2,2,2,1,z); 
  check_exit(l,3,2,2,1,z); 
  check_exit(l,0,3,2,1,m); //
  check_exit(l,1,3,2,1,z); 
  check_exit(l,2,3,2,1,z); 
  check_exit(l,3,3,2,1,z); 
  //
  check_exit(l,0,0,3,1,z); 
  check_exit(l,1,0,3,1,z); 
  check_exit(l,2,0,3,1,m); //
  check_exit(l,3,0,3,1,z); 
  check_exit(l,0,1,3,1,z); 
  check_exit(l,1,1,3,1,z); 
  check_exit(l,2,1,3,1,z); 
  check_exit(l,3,1,3,1,z); 
  check_exit(l,0,2,3,1,p); //
  check_exit(l,1,2,3,1,z); 
  check_exit(l,2,2,3,1,z); 
  check_exit(l,3,2,3,1,z); 
  check_exit(l,0,3,3,1,z); 
  check_exit(l,1,3,3,1,z); 
  check_exit(l,2,3,3,1,z); 
  check_exit(l,3,3,3,1,z); 
  //*
  check_exit(l,0,0,0,2,z); 
  check_exit(l,1,0,0,2,z); 
  check_exit(l,2,0,0,2,z); 
  check_exit(l,3,0,0,2,z); 
  check_exit(l,0,1,0,2,z); 
  check_exit(l,1,1,0,2,z); 
  check_exit(l,2,1,0,2,z); 
  check_exit(l,3,1,0,2,p); //
  check_exit(l,0,2,0,2,z); 
  check_exit(l,1,2,0,2,z); 
  check_exit(l,2,2,0,2,z); 
  check_exit(l,3,2,0,2,z); 
  check_exit(l,0,3,0,2,z); 
  check_exit(l,1,3,0,2,m); //
  check_exit(l,2,3,0,2,z); 
  check_exit(l,3,3,0,2,z); 
  //
  check_exit(l,0,0,1,2,z); 
  check_exit(l,1,0,1,2,z); 
  check_exit(l,2,0,1,2,z); 
  check_exit(l,3,0,1,2,m); //
  check_exit(l,0,1,1,2,z); 
  check_exit(l,1,1,1,2,z); 
  check_exit(l,2,1,1,2,z); 
  check_exit(l,3,1,1,2,z); 
  check_exit(l,0,2,1,2,z); 
  check_exit(l,1,2,1,2,z); 
  check_exit(l,2,2,1,2,z); 
  check_exit(l,3,2,1,2,z); 
  check_exit(l,0,3,1,2,p); //
  check_exit(l,1,3,1,2,z); 
  check_exit(l,2,3,1,2,z); 
  check_exit(l,3,3,1,2,z); 
  //
  check_exit(l,0,0,2,2,z); 
  check_exit(l,1,0,2,2,z); 
  check_exit(l,2,0,2,2,z); 
  check_exit(l,3,0,2,2,z); 
  check_exit(l,0,1,2,2,z); 
  check_exit(l,1,1,2,2,z); 
  check_exit(l,2,1,2,2,z); 
  check_exit(l,3,1,2,2,z); 
  check_exit(l,0,2,2,2,z); 
  check_exit(l,1,2,2,2,z); 
  check_exit(l,2,2,2,2,z); 
  check_exit(l,3,2,2,2,z); 
  check_exit(l,0,3,2,2,z); 
  check_exit(l,1,3,2,2,z); 
  check_exit(l,2,3,2,2,z); 
  check_exit(l,3,3,2,2,z); 
  //
  check_exit(l,0,0,3,2,z); 
  check_exit(l,1,0,3,2,p); //
  check_exit(l,2,0,3,2,z); 
  check_exit(l,3,0,3,2,z); 
  check_exit(l,0,1,3,2,m); //
  check_exit(l,1,1,3,2,z); 
  check_exit(l,2,1,3,2,z); 
  check_exit(l,3,1,3,2,z); 
  check_exit(l,0,2,3,2,z); 
  check_exit(l,1,2,3,2,z); 
  check_exit(l,2,2,3,2,z); 
  check_exit(l,3,2,3,2,z); 
  check_exit(l,0,3,3,2,z); 
  check_exit(l,1,3,3,2,z); 
  check_exit(l,2,3,3,2,z); 
  check_exit(l,3,3,3,2,z); 
  //*
  check_exit(l,0,0,0,3,z); 
  check_exit(l,1,0,0,3,z); 
  check_exit(l,2,0,0,3,z); 
  check_exit(l,3,0,0,3,z); 
  check_exit(l,0,1,0,3,z); 
  check_exit(l,1,1,0,3,z); 
  check_exit(l,2,1,0,3,m); //
  check_exit(l,3,1,0,3,z); 
  check_exit(l,0,2,0,3,z); 
  check_exit(l,1,2,0,3,p); //
  check_exit(l,2,2,0,3,z); 
  check_exit(l,3,2,0,3,z); 
  check_exit(l,0,3,0,3,z); 
  check_exit(l,1,3,0,3,z); 
  check_exit(l,2,3,0,3,z); 
  check_exit(l,3,3,0,3,z); 
  //
  check_exit(l,0,0,1,3,z); 
  check_exit(l,1,0,1,3,z); 
  check_exit(l,2,0,1,3,p); //
  check_exit(l,3,0,1,3,z); 
  check_exit(l,0,1,1,3,z); 
  check_exit(l,1,1,1,3,z); 
  check_exit(l,2,1,1,3,z); 
  check_exit(l,3,1,1,3,z); 
  check_exit(l,0,2,1,3,m); //
  check_exit(l,1,2,1,3,z); 
  check_exit(l,2,2,1,3,z); 
  check_exit(l,3,2,1,3,z); 
  check_exit(l,0,3,1,3,z); 
  check_exit(l,1,3,1,3,z); 
  check_exit(l,2,3,1,3,z); 
  check_exit(l,3,3,1,3,z); 
  //
  check_exit(l,0,0,2,3,z); 
  check_exit(l,1,0,2,3,m); //
  check_exit(l,2,0,2,3,z); 
  check_exit(l,3,0,2,3,z); 
  check_exit(l,0,1,2,3,p); //
  check_exit(l,1,1,2,3,z); 
  check_exit(l,2,1,2,3,z); 
  check_exit(l,3,1,2,3,z); 
  check_exit(l,0,2,2,3,z); 
  check_exit(l,1,2,2,3,z); 
  check_exit(l,2,2,2,3,z); 
  check_exit(l,3,2,2,3,z); 
  check_exit(l,0,3,2,3,z); 
  check_exit(l,1,3,2,3,z); 
  check_exit(l,2,3,2,3,z); 
  check_exit(l,3,3,2,3,z); 
  //
  check_exit(l,0,0,3,3,z); 
  check_exit(l,1,0,3,3,z); 
  check_exit(l,2,0,3,3,z); 
  check_exit(l,3,0,3,3,z); 
  check_exit(l,0,1,3,3,z); 
  check_exit(l,1,1,3,3,z); 
  check_exit(l,2,1,3,3,z); 
  check_exit(l,3,1,3,3,z); 
  check_exit(l,0,2,3,3,z); 
  check_exit(l,1,2,3,3,z); 
  check_exit(l,2,2,3,3,z); 
  check_exit(l,3,2,3,3,z); 
  check_exit(l,0,3,3,3,z); 
  check_exit(l,1,3,3,3,z); 
  check_exit(l,2,3,3,3,z); 
  check_exit(l,3,3,3,3,z); 


  if (nchk != 4*4*4*4)
  {
    std::cout << __func__ << " incomplete test " << std::endl;
  }

  if ( nonz != 4*3*2*1 )
  {
    std::cout << __func__ << " missing permutations " << std::endl;
    std::cout << __func__ << " nperm = " << nonz << std::endl;
  }

  nchk = 0; 
  nonz = 0; 
}


template<typename T>
void test_levi(void)
{
  radmat::Tensor<T,4> l = radmat::levi_civita<T,4>(); 
  test_levi(l); 
}


int main(void)
{
  test_levi<double>();
  test_levi<std::complex<double> >(); 
  test_levi<int>();  

  return 0; 
}
