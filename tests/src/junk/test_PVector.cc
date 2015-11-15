// test_PVector.cc -
//
// Monday, April  2 2012
//

#include "radmat/polarisation/polarisationbase.h"
#include "radmat/tensor/tensorbase.h"
#include "radmat/tensor/test_utils.h"
#include "radmat/utils/pow2assert.h"
#include "radmat/utils/breit_frame.h"
#include "itpp/itbase.h"
#include <iostream>


using namespace radmat;

int 
C1(void)
{
  // test a whole bunch of random momentum using the Coupler<1> class
  unsigned int ct(0);
  const double precision = 1e-13;
  for(int count = -50; count < 51; ++count)
    {
      for(short h = -1; h < 2; ++h)
	{
	  ++ct;
	  Tensor<double,1> pmu, p_z;
	  std::vector<idx_t> dim(1,4);
	  pmu.create(&dim[0]);

	  pmu[1] = itpp::randu()*count + itpp::randu();
	  pmu[2] = itpp::randu()*count + itpp::randu();
	  pmu[3] = itpp::randu()*count + itpp::randu();
	  pmu[0] = sqrt(pmu[1]*pmu[1] + 
			pmu[2]*pmu[2] + 
			pmu[3]*pmu[3]) + itpp::randu(); // keep the vector timelike

	  // doesn't make sense
	  if((pmu[1] == pmu[2]) && (pmu[2] == pmu[3]) && (pmu[3] == 0))
	    {
	      std::cout << "skipping" << std::endl;
	      --ct;
	      continue;
	    }

	  p_z = pmu;
	  p_z[1] = 0.;
	  p_z[2] = 0.;
	  p_z[3] = sqrt(pmu[1]*pmu[1] + 
			pmu[2]*pmu[2] + 
			pmu[3]*pmu[3]
			);

	  itpp::Mat<double> R = breit::rodRotMat(p_z,pmu);

	  for(idx_t i = 0; i < 3; ++i)
	    {
	      POW2_ASSERT( fabs(pmu[i+1] - R(i,2)*p_z[3]) < precision);
	    }
	  pFacKey k(pmu,1,h), k_z(p_z,1,h);
	  pFac factory;

	  Coupler<1> c(k,&factory), c_z(k_z,&factory);
      
	  Tensor<std::complex<double> , 1> *eps , *eps_z;
      
	  // need to be called in this order to avoid a leak with valgrind b/c of 
	  // operator() definition in Coupler<1>
	  eps_z = c_z();
	  eps = c();
      
	  Tensor<std::complex<double> ,1 > e,e_z, er;
	  e = * eps;
	  e_z = *eps_z;

	  delete eps;
	  delete eps_z;

	  er = e_z;

	  for(idx_t i = 0; i < 3; ++i)
	    {
	      er[i+1] = 0;
	      er[i+1] += R(i,0)*e_z[1];
	      er[i+1] += R(i,1)*e_z[2];
	      er[i+1] += R(i,2)*e_z[3];
	    }

	  for(idx_t i = 0; i < 4; ++i)
	    POW2_ASSERT(fabs(e[i] - er[i]) < precision);
	}
    }

  std::cout << ct <<  " J = 1 Coupler<J> polarisation vector tests passed" << std::endl;

  return 0;
}

int
factory(void)
{
  pFac factory;
  unsigned int ct(0);
  // test a whole bunch of random momentum using the pFac factory class
  const double precision = 1e-13;
  for(double count = -50; count < 51; ++count)
    {
      for(short h = -1; h < 2; ++h)
	{
	  ++ct;
	  Tensor<double,1> pmu, p_z;
	  std::vector<idx_t> dim(1,4);
	  pmu.create(&dim[0]);

	  pmu[1] = itpp::randu()*count + itpp::randu();
	  pmu[2] = itpp::randu()*count + itpp::randu();
	  pmu[3] = itpp::randu()*count + itpp::randu();
	  pmu[0] = sqrt(pmu[1]*pmu[1] + 
			pmu[2]*pmu[2] + 
			pmu[3]*pmu[3]) + itpp::randu();


			
	  // doesn't make sense
	  if((pmu[1] == pmu[2]) && (pmu[2] == pmu[3]) && (pmu[3] == 0))
	    {
	      std::cout << "skipping" << std::endl;
	      --ct;
	      continue;
	    }

	  p_z = pmu;
	  p_z[1] = 0.;
	  p_z[2] = 0.;
	  p_z[3] = sqrt(pmu[1]*pmu[1] + 
			pmu[2]*pmu[2] + 
			pmu[3]*pmu[3]
			);

	  itpp::Mat<double> R = breit::rodRotMat(p_z,pmu);

	  for(idx_t i = 0; i < 3; ++i)
	    {
	      POW2_ASSERT( fabs(pmu[i+1] - R(i,2)*p_z[3]) < precision);
	    }
	  pFacKey k(pmu,1,h), k_z(p_z,1,h);

      
	  Tensor<std::complex<double> ,1 > e,e_z, er;
	  e = downcastAndDeref<Tensor<std::complex<double> ,1 >, TensorImplBase>( factory.get(k) );
	  e_z = downcastAndDeref<Tensor<std::complex<double> ,1 >, TensorImplBase>( factory.get(k_z) );
	  er = e_z;

	  for(idx_t i = 0; i < 3; ++i)
	    {
	      er[i+1] = 0;
	      er[i+1] += R(i,0)*e_z[1];
	      er[i+1] += R(i,1)*e_z[2];
	      er[i+1] += R(i,2)*e_z[3];
	    }

	  for(idx_t i = 0; i < 4; ++i)
	    POW2_ASSERT(fabs(e[i] - er[i]) < precision);
	}
    }

  std::cout << ct << " J = 1 factory polarisation vector tests passed" << std::endl;

  // polarisation::pFacInv::dumpInventory();

  return 0;

}



int 
main(void)
{
  return C1() + factory();
}
