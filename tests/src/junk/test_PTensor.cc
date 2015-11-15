// test_PTensor.cc -
//
// Tuesday, April  3 2012
//

#include "radmat/tensor/tensorbase.h"
#include "radmat/polarisation/polarisationbase.h"
#include "radmat/utils/breit_frame.h"
#include "radmat/utils/pow2assert.h"
#include "hadron/clebsch.h"
#include "radmat/tensor/test_utils.h"

using namespace radmat;
using namespace radmat::breit;

int 
main(void)
{
  
  double precision = 1e-14;

  pFac factory;

  //short HELICITY = 2;
  for(short HELICITY = -2; HELICITY < 3; ++HELICITY)
    {

      Tensor<double,1> pmu, p_z;
      std::vector<idx_t> dim(1,4);
      pmu.create(&dim[0]);

      pmu[1] = itpp::randu() + itpp::randu();
      pmu[2] = itpp::randu() + itpp::randu();
      pmu[3] = itpp::randu() + itpp::randu();
      pmu[0] = sqrt(pmu[1]*pmu[1] + 
		    pmu[2]*pmu[2] + 
		    pmu[3]*pmu[3]) + itpp::randu();

      p_z = pmu;
      p_z[1] = 0.;
      p_z[2] = 0.;
      p_z[3] = sqrt(pmu[1]*pmu[1] + 
		    pmu[2]*pmu[2] + 
		    pmu[3]*pmu[3]
		    );

      itpp::Mat<double> R = rodRotMat(p_z,pmu);

      for(idx_t i = 0; i < 3; ++i)
	{
	  POW2_ASSERT( fabs(pmu[i+1] - R(i,2)*p_z[3]) < precision);
	}


      std::map<short,Tensor<std::complex<double> , 1>* > h_map;
      for(short h = -1; h < 2; ++h)
	{
	  pFacKey key(p_z,1,h); 
	  h_map[h] = dynamic_cast<Tensor<std::complex<double> ,1>* >( factory.get( key ) );
	}

      // get the one from the factory
      pFacKey P2_z(p_z,2,HELICITY);
      Tensor<std::complex<double> , 2> p2z = downcastAndDeref<Tensor<std::complex<double>,2 >, TensorImplBase>( factory.get(P2_z) );
      Tensor<std::complex<double> , 2> fake;
      std::vector<idx_t> dimfake(2,4);
      fake.create(&dimfake[0]);

      // cook up a fake one using cg
      for(short h1 = -1; h1 < 2; ++h1)
	for(short h2 = -1; h2 < 2; ++h2)
	  for(idx_t mu = 0; mu < 4; ++mu)
	    for(idx_t nu =0; nu < 4; ++nu)
	      {
		fake[mu][nu] += (Hadron::clebsch(2,2*h1,2,2*h2,4,2*HELICITY)
				 *(*h_map[h1])[mu] 
				 *(*h_map[h2])[nu]
				 );

	      }   

      // check that they are the same
      for(idx_t mu = 0; mu < 4; ++mu)
	for(idx_t nu =0; nu < 4; ++nu)
	  POW2_ASSERT( fabs(fake[mu][nu] - p2z[mu][nu]) < precision);

    }

  return 0;
}
