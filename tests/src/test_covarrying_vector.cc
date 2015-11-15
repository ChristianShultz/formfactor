// test_covarrying_vector.cc -
//
// Tuesday, June  5 2012
//

#include "../headers/test_common.h"
#include "../headers/tester.h"
#include "radmat/fake_data/covarrying_vectors.h"
#include "itpp/itbase.h"
#include <iostream>
#include <complex>
#include <vector>
#include <stdlib.h>
#include <time.h>


#define SIZE 2                       // the problem gets more difficult for larger rank..
#define NCFG 50000
#define PRECISION_CHECK_MEAN 1e-2    // allow for 1% error in each in 50k draws
#define PRECISION_CHECK_COR 1e-2     // ~ this method gives a realistic data set but
#define NCHECK 10                    // it converges to the target painfully slowly 

/*
  NB: for O(500) cfg this will only loosely resemble the target distribution.
    Another thing to think about is to draw in this manner using some sort of markov 
    chain where the probability of accepting the drawn vector is related to the 
    covariance of the current set. This would converge faster at the expense of 
    autocorrelation and randomness (this might actually resemble real data more since
    we have some ammount of autocorrelation coming from gauge configurations..).
    Possibly set some long burn time to reduce this.. 
 */

namespace radmat
{


  template<> 
  tester test_covarrying_vectors::test<double>(void) const
  {
    tester m_test(__PRETTY_FUNCTION__);
    srand(time(NULL));
    itpp::RNG_reset(rand()); // reseed itpp rng with time

    for(int i = 0; i < NCHECK; i ++)
      {
	// wrapRandN wrap;
	// std::cout << wrap.wrap<double>(SIZE) << std::endl;
	// std::cout << wrap.wrap<std::complex<double> >(SIZE) << std::endl;

	corMat c;
  
	// std::cout << itpp::eig_sym(c.genMat<double>(SIZE)) << std::endl;
	// std::cout << itpp::eig_sym(c.genMat<cd::dc>(SIZE)) << std::endl;

	itpp::Vec<double> mean = itpp::randu(SIZE);
	itpp::Vec<double> var = itpp::randu(SIZE);
	itpp::Mat<double> corr = c.genMat<double>(SIZE);
  
	for(int i = 0; i < SIZE; i++)
	  var[i] *= mean[i] *0.2;     // 20% max variance


	itpp::Vec<double> rootVar(SIZE);
	for(int i = 0; i < SIZE; i++)
	  rootVar(i) = std::sqrt(var(i));
	itpp::Mat<double> covIn = itpp::diag(rootVar)* corr * itpp::diag(rootVar);  

	itpp::Vec<double> dmean(SIZE),dvar(SIZE);
	dmean.zeros();
	dvar.zeros();

	std::vector<itpp::Vec<double> > cv = genCovarryingDist(mean,var,corr,NCFG);

	std::vector<itpp::Vec<double> >::const_iterator it;
	for(it = cv.begin(); it != cv.end(); it++)
	  dmean += *it;
  
	dmean /= double(NCFG);

	for(it = cv.begin(); it != cv.end(); it++)
	  dvar += itpp::elem_mult((*it -dmean),(*it -dmean));

	dvar /= double(NCFG);
  
	wrapCor wCor;
	wrapCov wCov;
	// std::cout << "target:\nmean: " << mean << std::endl;
	// std::cout << "var:  " << var << std::endl;
	// std::cout << "drawn:\nmean: " << dmean << std::endl;
	// std::cout << "var:  " << dvar << std::endl;
	// std::cout << "\n correlation in: \n" << corr << std::endl;
	// std::cout << "\n covariance in: \n" << covIn << std::endl;
	// std::cout << "\n covariance out: \n" << wCov.cov<double>(cv) << std::endl;
	// std::cout << "\n correlation out: \n" <<  wCor.cor<double>(cv) << std::endl;


	itpp::Vec<double> rmean = itpp::round_to_zero(mean - dmean,PRECISION_CHECK_MEAN);
	itpp::Vec<double> zero(mean);
	zero.zeros();
	TESTER_TEST(m_test, zero == rmean, "didn't draw the correct mean");
    
	itpp::Vec<double> rvar = itpp::round_to_zero(var - dvar, PRECISION_CHECK_MEAN);
	TESTER_TEST(m_test, zero == rvar, "didn't draw the correct variance");

	itpp::Mat<double> rcor = itpp::round_to_zero(corr - wCor.cor<double>(cv), PRECISION_CHECK_COR);
	itpp::Mat<double> Mzero(rcor);
	Mzero.zeros();
	TESTER_TEST(m_test, Mzero == rcor, "didn't draw the correct correlation matrix");

	itpp::Mat<double> rcov = itpp::round_to_zero(covIn - wCov.cov<double>(cv),PRECISION_CHECK_COR);
	TESTER_TEST(m_test,Mzero == rcov, "didn't draw the correct covariance");

	if(rcov != Mzero)
	  {
	    std::cout << covIn - wCov.cov<double>(cv) << "\n"
		      << itpp::round_to_zero(covIn - wCov.cov<double>(cv), PRECISION_CHECK_COR) 
		      << std::endl;
	    std::cout << "\n covariance in: \n" << covIn << std::endl;
	    std::cout << "\n covariance out: \n" << wCov.cov<double>(cv) << std::endl;
	  }

	if(rcor != Mzero)
	  {
	    std::cout << corr - wCor.cor<double>(cv) << "\n"
		      << itpp::round_to_zero(corr - wCor.cor<double>(cv) )
		      << std::endl;
	    std::cout << "\n correlation in: \n" << corr << std::endl;
	    std::cout << "\n correlation out: \n" <<  wCor.cor<double>(cv) << std::endl;
	  }

      }
    return m_test;
  }


  template<>
  tester test_covarrying_vectors::test<std::complex<double> >(void) const
  {
    tester m_test(__PRETTY_FUNCTION__);
    srand(time(NULL));
    itpp::RNG_reset(rand()); // reseed itpp rng with time

    for(int i = 0; i < NCHECK; i++)
      {
	wrapRandN wrap;
	corMat c;
	itpp::Vec<cd::dc> mean = wrap.wrap<cd::dc>(SIZE);
	itpp::Vec<double> var = itpp::randu(SIZE);
	itpp::Mat<cd::dc> corr = c.genMat<cd::dc>(SIZE);

	for(int i =0; i < SIZE; i++)
	  var[i] *= std::sqrt(std::norm(mean[i]))*0.2;


	itpp::Vec<double> rootVar(SIZE);
	for(int i = 0; i < SIZE; i++)
	  rootVar(i) = std::sqrt(var(i));
	itpp::Mat<cd::dc> covIn = itpp::diag(rootVar)*corr*itpp::diag(rootVar);

	itpp::Vec<cd::dc> dmean(SIZE), dcvar(SIZE);
	itpp::Vec<double> dvar(SIZE);
	dmean.zeros(); 
	dvar.zeros();
	dcvar.zeros();

	std::vector<itpp::Vec<cd::dc> > cv = genCovarryingDist(mean,var,corr,NCFG);
	std::vector<itpp::Vec<cd::dc> >::const_iterator it;
	for(it = cv.begin(); it != cv.end(); it++)
	  dmean += *it;

	dmean /= double(NCFG);

	for(it = cv.begin(); it != cv.end(); it++)
	  dcvar += itpp::elem_mult((*it - dmean),itpp::conj(*it - dmean));

	dcvar /= double(NCFG);

	for(int i =0; i < SIZE; i++)
	  {
	    double real = std::real(dcvar(i));
	    double complex = std::imag(dcvar(i));
	    TESTER_TEST(m_test, fabs(complex/real) < PRECISION_CHECK_MEAN,
			"phasing issues are occurring");
	    dvar(i) = real;
	  }

	wrapCor wCor;
	wrapCov wCov;
	
	itpp::Vec<cd::dc> rmean = itpp::round_to_zero(mean - dmean, PRECISION_CHECK_MEAN);
	itpp::Vec<cd::dc> czero(rmean);
	czero.zeros();
	TESTER_TEST(m_test,czero == rmean, "didn't draw the correct mean");
	
	itpp::Vec<double> rvar = itpp::round_to_zero(var - dvar, PRECISION_CHECK_MEAN);
	itpp::Vec<double> rzero(rvar);
	rzero.zeros();
	TESTER_TEST(m_test, rzero == rvar, "didn't draw the correct variance");

	itpp::Mat<cd::dc> rcor = itpp::round_to_zero(corr - wCor.cor<cd::dc>(cv),PRECISION_CHECK_COR);
	itpp::Mat<cd::dc> Mczero(rcor);
	Mczero.zeros();
	TESTER_TEST(m_test,Mczero == rcor, "didn't draw the correct correlation matrix");

	itpp::Mat<cd::dc> rcov = itpp::round_to_zero(covIn - wCov.cov<cd::dc>(cv), PRECISION_CHECK_COR);
	TESTER_TEST(m_test,rcov == Mczero, "didn't draw the correct covariance matrix");

	if(rcov != Mczero)
	  {
	    std::cout << covIn - wCov.cov<cd::dc>(cv) << "\n"
		      << itpp::round_to_zero(covIn - wCov.cov<cd::dc>(cv), PRECISION_CHECK_COR) 
		      << std::endl;
	    std::cout << "\n covariance in: \n" << covIn << std::endl;
	    std::cout << "\n covariance out: \n" << wCov.cov<cd::dc>(cv) << std::endl;
	  }

	if(rcor != Mczero)
	  {
	    std::cout << corr - wCor.cor<cd::dc>(cv) << "\n"
		      << itpp::round_to_zero(corr - wCor.cor<cd::dc>(cv), PRECISION_CHECK_COR)
		      << std::endl;
	    std::cout << "\n correlation in: \n" << corr << std::endl;
	    std::cout << "\n correlation out: \n" << wCor.cor<cd::dc>(cv) << std::endl;
	  }

	/*
	  std::cout << "target:\nmean: " << mean << std::endl;
	  std::cout << "var:  " << var << std::endl;
	  std::cout << "drawn:\nmean: " << dmean << std::endl;
	  std::cout << "var:  " << dvar << std::endl;
	  std::cout << "\n correlation in: \n" << corr << std::endl;
	  std::cout << "\n covariance in: \n" << covIn << std::endl;
	  std::cout << "\n covariance out: \n" << wCov.cov<cd::dc>(cv) << std::endl;
	  std::cout << "\n correlation out: \n" <<  wCor.cor<cd::dc>(cv) << std::endl;
	*/
      }

    return m_test;
  }
  





  /*
    tester test_covarrying_vector(void)
    {
    tester m_test(__func__);


    srand(time(NULL));
    itpp::RNG_reset(rand()); // reseed itpp rng with time

    
    double val,val2,v;
    std::complex<double> valc,valc2,vc;
    wrapRandN w;
    val = val2 = 0;
    valc = valc2 = 0;

    itpp::Vec<double> vv = w.wrap<double>(NCFG);
    itpp::Vec<std::complex<double> > vvcc = w.wrap<std::complex<double> >(NCFG);

    for(int i = 0; i < NCFG; i++)
    {
    v = vv[i];
    vc = vvcc[i];
    val += v;
    val2 +=v*v;
    valc += vc;
    valc2 += vc * std::conj(vc);
    }

    val /= double(NCFG);
    val2 /= double(NCFG);
    valc /= double(NCFG);
    valc2 /= double(NCFG);
  
    std::cout << NCFG << "draws :" << std::endl;
    std::cout << "<val> = " << val << std::endl;
    std::cout << "<val^2> = " << val2 << std::endl;
    std::cout << "<valc> = " << valc << std::endl;
    std::cout << "<valc^2> = " << valc2 << std::endl;
    

    wrapRandN wrap;
    // std::cout << wrap.wrap<double>(SIZE) << std::endl;
    // std::cout << wrap.wrap<std::complex<double> >(SIZE) << std::endl;

    corMat c;
  
    // std::cout << itpp::eig_sym(c.genMat<double>(SIZE)) << std::endl;
    // std::cout << itpp::eig_sym(c.genMat<cd::dc>(SIZE)) << std::endl;

    itpp::Vec<double> mean = itpp::randu(SIZE);
    itpp::Vec<double> var = itpp::randu(SIZE);
    itpp::Mat<double> corr = c.genMat<double>(SIZE);
  
    for(int i = 0; i < SIZE; i++)
    var[i] *= mean[i] *0.2;     // 20% max variance


    itpp::Vec<double> rootVar(SIZE);
    for(int i = 0; i < SIZE; i++)
    rootVar(i) = std::sqrt(var(i));
    itpp::Mat<double> covIn = itpp::diag(rootVar)* corr * itpp::diag(rootVar);  

    itpp::Vec<double> dmean(SIZE),dvar(SIZE);
    dmean.zeros();
    dvar.zeros();

    std::vector<itpp::Vec<double> > cv = genCovarryingDist(mean,var,corr,NCFG);

    std::vector<itpp::Vec<double> >::const_iterator it;
    for(it = cv.begin(); it != cv.end(); it++)
    dmean += *it;
  
    dmean /= double(NCFG);

    for(it = cv.begin(); it != cv.end(); it++)
    dvar += itpp::elem_mult((*it -mean),(*it -mean));

    dvar /= double(NCFG);
  
    
    std::cout << "target:\nmean: " << mean << std::endl;
    std::cout << "var:  " << var << std::endl;

    std::cout << "drawn:\nmean: " << dmean << std::endl;
    std::cout << "var:  " << dvar << std::endl;

    wrapCor wCor;
    wrapCov wCov;
    std::cout << "\n correlation in: \n" << corr << std::endl;
    std::cout << "\n covariance in: \n" << covIn << std::endl;
    std::cout << "\n covariance out: \n" << wCov.cov<double>(cv) << std::endl;
    std::cout << "\n correlation out: \n" <<  wCor.cor<double>(cv) << std::endl;
    

    itpp::Vec<double>


    return m_test;
    }
  */

} // close radmat namespace
