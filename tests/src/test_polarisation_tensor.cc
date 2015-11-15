// test_polarisation_tensor.cc -
//
// Wednesday, May  9 2012
//

#include "../headers/test_common.h"
#include "../headers/tester.h"
#include "radmat/utils/polarisation_tensors.h"
#include "radmat/utils/tensor.h"
#include "xml_array.h"
#include "hadron/clebsch.h"
#include "radmat/utils/pow2assert.h"
#include "itpp/itbase.h"
#include <iostream>

typedef XMLArray::Array<int> Aint;

double genE(const Aint mom, double m)
{
  return sqrt(m*m + mom[0]*mom[0] + mom[1]*mom[1] + mom[2]*mom[2]);
}

radmat::Tensor<std::complex<double>, 1 > genEps(const double E, const short h, const Aint mom)
{
  radmat::Tensor<std::complex<double> , 1> eps((radmat::TensorShape<1>())[4],0.);
  radmat::Tensor<double,2> R = radmat::genRotationMatrix(mom);
  if(h == 1)
  {
    eps[1] = -1./sqrt(2.);
    eps[2] = std::complex<double>(0.,-1./sqrt(2));
    return R*eps;
  }
  if(h == -1)
  {
    eps[1] = 1./sqrt(2.);
    eps[2] = std::complex<double>(0.,-1./sqrt(2));
    return R*eps;  
  }
  if(h == 0)
  {
    double m = sqrt(E*E - (mom[0]*mom[0] + mom[1]*mom[1] + mom[2]*mom[2]));
    double p = sqrt(mom[0]*mom[0] + mom[1]*mom[1] + mom[2]*mom[2]);
    eps[0] = p/m;
    eps[3] = E/m;
    return R*eps;
  }

  POW2_ASSERT(false);
  return eps; // dummy to get the compiler to shut up
}

// Lists of momentum directions
const int momListD2[][3] = {{1, 1, 0}, {0, 1, 1}, {1, 0, 1}, {1, -1, 0}, 
  {0, 1, -1}, {-1, 0, 1}, {-1, 1, 0}, {0, -1, 1}, {1, 0, -1}, {-1, -1, 0},
  {0, -1, -1}, {-1, 0, -1}};
const int dimMomListD2 = 12;

const int momListD3[][3] = {{1, 1, 1}, {-1, 1, 1}, {1, -1, 1}, {1, 1, -1},
  {-1, -1, 1}, {1, -1, -1}, {-1, 1, -1}, {-1, -1, -1}};
const int dimMomListD3 = 8;

const int momListC4nm0[][3] = {{0, 1, 2}, {1, 2, 0}, {2, 0, 1}, {0, 2, 1},
  {2, 1, 0}, {1, 0, 2}, {0, 1, -2}, {1, -2, 0}, {-2, 0, 1}, {0, -2, 1},
  {-2, 1, 0}, {1, 0, -2}, {0, -1, 2}, {-1, 2, 0}, {2, 0, -1}, {0, 2, -1},
  {2, -1, 0}, {-1, 0, 2}, {0, -1, -2}, {-1, -2, 0}, {-2, 0, -1}, {0, -2, -1}, 
  {-2, -1, 0}, {-1, 0, -2}};
const int dimMomListC4nm0 = 24;

const int momListC4nnm[][3] = {{1, 1, 2}, {1, 2, 1}, {2, 1, 1}, {-1, 1, 2},
  {-1, 2, 1}, {1, -1, 2}, {1, 2, -1}, {2, -1, 1}, {2, 1, -1}, {1, 1, -2},
  {1, -2, 1}, {-2, 1, 1}, {-1, -1, 2}, {-1, 2, -1}, {2, -1, -1}, {-1, 1, -2},
  {-1, -2, 1}, {1, -1, -2}, {1, -2, -1}, {-2, -1, 1}, {-2, 1, -1}, {-1, -1, -2},
  {-1, -2, -1}, {-2, -1, -1}};
const int dimMomListC4nnm = 24;


struct triplet
{
  int x,y,z;
  triplet(int a, int b, int c)
    : x(a) , y(b) , z(c)
  {  }
};


std::vector<triplet> genMomList(void)
{
  std::vector<triplet> ret;
  for(int i = 0; i < dimMomListD2; ++i)
    ret.push_back(triplet(momListD2[i][0],momListD2[i][1],momListD2[i][2]));
  for(int i = 0; i < dimMomListD3; ++i)
    ret.push_back(triplet(momListD3[i][0],momListD3[i][1],momListD3[i][2]));
  for(int i = 0; i < dimMomListC4nm0; ++i)
    ret.push_back(triplet(momListC4nm0[i][0],momListC4nm0[i][1],momListC4nm0[i][2]));
  for(int i = 0; i < dimMomListC4nnm; ++i)
    ret.push_back(triplet(momListC4nnm[i][0],momListC4nnm[i][1],momListC4nnm[i][2]));

  return ret;
}


namespace radmat
{
  // something of an inductive test, we test the base (J=1) and then test J=2,J=3 and assumme since it is a recursive 
  // pattern that it will work for all J
  tester test_polarisation_tensor(void)
  {
    tester m_test(__func__);


    Aint mom(3);
    double precision = 1e-10;

    std::vector<triplet> moms = genMomList();
    std::vector<triplet>::const_iterator it;

    // start at J = 1
    for(it = moms.begin(); it != moms.end(); it++)
    {
      mom[0] = it->x;
      mom[1] = it->y;
      mom[2] = it->z;

      double m = itpp::randu() + 0.01;
      genPolTens<1> first_rank(mom);

      for(short i = -1; i < 2; ++i)
      {
        bool eq = equals(first_rank(genE(mom,m),i,1.), 
            genEps(genE(mom,m),i,mom),
            std::complex<double>(precision,precision)) ;
        TESTER_TEST(m_test,eq,"construction of base failed");

        if(!!!eq) 
          std::cout << first_rank(genE(mom,m),i,1.)- genEps(genE(mom,m),i,mom) << std::endl;
      }
    }

    Tensor<std::complex<double> , 2> coupled((TensorShape<2>())[4][4],0.);
    std::complex<double> zero(0.,0.);

    // check J = 2 explicitly
    for(it = moms.begin(); it != moms.end(); ++it)
    {
      mom[0] = it->x;
      mom[1] = it->y;
      mom[2] = it->z;

      double m = itpp::randu() + 0.01;
      genPolTens<1> first_rank(mom);
      genPolTens<2> second_rank(mom);

      for(short i = -2; i <3; ++i)
      {
        coupled.fill(zero);
        for(short h1 = -1; h1 < 2; ++h1)
          for(short h2 = -1; h2 < 2; ++h2)
            coupled += Hadron::clebsch(2,2*h1,2,2*h2,4,2*i) * 
              (first_rank(genE(mom,m),h1,1.)  ^ first_rank(genE(mom,m),h2,1.)) ;  // ^ is tensor product

        bool eq = equals(coupled,second_rank(genE(mom,m),i,1.),std::complex<double>(precision,precision));
        TESTER_TEST(m_test,eq,"construction of rank 2 failed");
      }
    }

    // check J = 3 explicitly
    Tensor<std::complex<double> , 3> coupled3((TensorShape<3>())[4][4][4],0.);
    for(it = moms.begin(); it != moms.end(); ++it)
    {
      mom[0] = it->x;
      mom[1] = it->y;
      mom[2] = it->z;

      double m = itpp::randu() + 0.01;
      genPolTens<1> first_rank(mom);
      genPolTens<2> second_rank(mom);
      genPolTens<3> third_rank(mom);

      for(short i = -3; i <4; ++i)
      {
        coupled3.fill(zero);
        for(short h1 = -1; h1 < 2; ++h1)
          for(short h2 = -2; h2 < 3; ++h2)
            coupled3 += Hadron::clebsch(2,2*h1,4,2*h2,6,2*i) * 
              (first_rank(genE(mom,m),h1,1.)  ^ second_rank(genE(mom,m),h2,1.)) ; // ^ is tensor product

        bool eq = equals(coupled3,third_rank(genE(mom,m),i,1.),std::complex<double>(precision,precision));
        TESTER_TEST(m_test,eq,"construction of rank 3 failed");
      }
    }

    return m_test;
  }




}
