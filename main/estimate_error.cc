/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

* File Name : estimate_error.cc

* Purpose :

* Creation Date : 08-10-2014

* Last Modified : Mon 13 Oct 2014 12:32:27 PM EDT

* Created By : shultz

_._._._._._._._._._._._._._._._._._._._._.*/
#include "itpp/itbase.h"
#include "ensem/ensem.h"
#include "semble/semble_semble.h"
#include "math.h"
#include <string>
#include <exception>
#include <vector>
#include <map>
#include <iostream>


struct START_PARS
{
  START_PARS()
    : fname(std::string("error_estimate"))
  {}

  int ncfg; 
  itpp::Mat<double> covariance; 
  itpp::Vec<double> mean; 
  itpp::Vec<double> x_sample; 
  std::string fname; 
  double (*F)(const double x, const itpp::Vec<double> pars); 
};

itpp::Vec<double> sample_range(double low, double high, double a)
{
  int sz = int( (high-low)/a );
  if( sz <= 0 ) 
    throw std::string("invalid range");

  itpp::Vec<double> ret( sz ); 
  for(int l = 0; l < sz; ++l)
    ret(l) = low + double(l)*a; 

  return ret; 
}


itpp::Mat<double>  construct_covariance(const itpp::Mat<double> &correlation, 
    const itpp::Vec<double> & error)
{
  itpp::Mat<double> covariance = correlation; 
  if( error.length() != correlation.cols() ) 
    throw std::string("invalid dimension 1"); 

  if( correlation.rows() != correlation.cols() ) 
    throw std::string("invalid dimension 2"); 

  for(int i = 0; i < covariance.rows(); ++i)
    for(int j = 0; j < covariance.cols(); ++j)
      covariance(i,j) *= sqrt(error(i)*error(j));

  return covariance; 
}

////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////

double PiForm(const double x, const itpp::Vec<double> p)
{
  return exp( - (x / 16. / (p(0)*p(0))) ); 
}

double RhoForm(const double x, const itpp::Vec<double> p)
{
  return p(0)*exp( - x / 16. / (p(1)*p(1)) ); 
}

double RhoPiForm(const double x, const itpp::Vec<double> p)
{
  return p(0)*exp( - (x / 16. / (p(1)*p(1)))*(1 + x*p(2) )); 
}

double PiPiStarForm(const double x, const itpp::Vec<double> p)
{
  return x*p(0)*p(0)*exp( - x*p(1)*p(1) ); 
}

double Polynomial2(const double x, const itpp::Vec<double> p)
{
  return p(0) + p(1)*x + p(2)*pow(x,2); 
}


////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////

START_PARS initialize_Gm(void)
{
  START_PARS foo;
  int ncfg = 500; 
  itpp::Vec<double> mean(2); 
  itpp::Vec<double> error(2); 
  itpp::Mat<double> correlation(2,2); 
  itpp::Vec<double> x_sample = sample_range(0,0.06,0.0001); 

  mean(0) = 2.0046; 
  error(0) = 0.01876; 
  mean(1) = 0.0588774; 
  error(1) = 0.000492; 
  correlation(0,0) = 1.; 
  correlation(1,0) = -0.903; 
  correlation(0,1) = -0.903; 
  correlation(1,1) = 1.; 

  itpp::Mat<double> covariance = construct_covariance(correlation,error); 

  foo.ncfg = ncfg; 
  foo.mean = mean; 
  foo.covariance = covariance; 
  foo.x_sample = x_sample; 
  foo.F = &RhoForm; 
  foo.fname = std::string("GmErrorEstimate.txt"); 

  return foo; 
}

START_PARS initialize_Gq(void)
{
  START_PARS foo;
  int ncfg = 500; 
  itpp::Vec<double> mean(2); 
  itpp::Vec<double> error(2); 
  itpp::Mat<double> correlation(2,2); 
  itpp::Vec<double> x_sample = sample_range(0,0.06,0.0001); 

  mean(0) = -0.381464; 
  error(0) = 0.04572; 
  mean(1) = 0.0529989; 
  error(1) = 0.005015; 
  correlation(0,0) = 1.; 
  correlation(1,0) = 0.936; 
  correlation(0,1) = 0.936; 
  correlation(1,1) = 1.; 

  itpp::Mat<double> covariance = construct_covariance(correlation,error); 

  foo.ncfg = ncfg; 
  foo.mean = mean; 
  foo.covariance = covariance; 
  foo.x_sample = x_sample; 
  foo.F = &RhoForm; 
  foo.fname = std::string("GqErrorEstimate.txt"); 

  return foo; 
}

START_PARS initialize_Pi_short(void)
{
  START_PARS foo;
  int ncfg = 1500; 
  itpp::Vec<double> mean(1); 
  itpp::Vec<double> error(1); 
  itpp::Mat<double> correlation(1,1); 
  itpp::Vec<double> x_sample = sample_range(0,0.0123,0.0001); 

  mean(0) = 0.0559166; 
  error(0) = 0.0004132; 
  correlation(0,0) = 1.; 

  itpp::Mat<double> covariance = construct_covariance(correlation,error); 

  foo.ncfg = ncfg; 
  foo.mean = mean; 
  foo.covariance = covariance; 
  foo.x_sample = x_sample; 
  foo.F = &PiForm; 
  foo.fname = std::string("ShortErrorEstimate.txt"); 

  return foo; 
}

START_PARS initialize_PiStar_short(void)
{
  START_PARS foo;
  int ncfg = 500; 
  itpp::Vec<double> mean(2); 
  itpp::Vec<double> error(2); 
  itpp::Mat<double> correlation(2,2); 
  itpp::Vec<double> x_sample = sample_range(0,0.0128281,0.0001); 

  mean(0) = 0.974356; 
  error(0) = 0.02557; 
  mean(1) = 0.0372104; 
  error(1) = 0.001093; 
  correlation(0,0) = 1.; 
  correlation(1,0) = -0.813; 
  correlation(0,1) = -0.813; 
  correlation(1,1) = 1.; 

  itpp::Mat<double> covariance = construct_covariance(correlation,error); 

  foo.ncfg = ncfg; 
  foo.mean = mean; 
  foo.covariance = covariance; 
  foo.x_sample = x_sample; 
  foo.F = &RhoForm; 
  foo.fname = std::string("ShortErrorEstimate.txt"); 

  return foo; 
}

START_PARS initialize_PiStar_long(void)
{
  START_PARS foo;
  int ncfg = 500; 
  itpp::Vec<double> mean(2); 
  itpp::Vec<double> error(2); 
  itpp::Mat<double> correlation(2,2); 
  itpp::Vec<double> x_sample = sample_range(0,0.0259126,0.0001); 

  mean(0) = 0.949467; 
  error(0) = 0.03461; 
  mean(1) = 0.0390253; 
  error(1) = 0.001041; 
  correlation(0,0) = 1.; 
  correlation(1,0) = -0.839; 
  correlation(0,1) = -0.839; 
  correlation(1,1) = 1.; 

  itpp::Mat<double> covariance = construct_covariance(correlation,error); 

  foo.ncfg = ncfg; 
  foo.mean = mean; 
  foo.covariance = covariance; 
  foo.x_sample = x_sample; 
  foo.F = &RhoForm; 
  foo.fname = std::string("LongErrorEstimate.txt"); 

  return foo; 
}

START_PARS initialize_PiPiStar(void)
{
  START_PARS foo;
  int ncfg = 500; 
  itpp::Vec<double> mean(2); 
  itpp::Vec<double> error(2); 
  itpp::Mat<double> correlation(2,2); 
  itpp::Vec<double> x_sample = sample_range(-0.015,0.03,0.0001); 

  mean(0) = 3.00576; 
  error(0) = 0.04111; 
  mean(1) = 5.77671; 
  error(1) = 0.1417; 
  correlation(0,0) = 1.; 
  correlation(1,0) = 0.904; 
  correlation(0,1) = 0.904; 
  correlation(1,1) = 1.; 

  itpp::Mat<double> covariance = construct_covariance(correlation,error); 

  foo.ncfg = ncfg; 
  foo.mean = mean; 
  foo.covariance = covariance; 
  foo.x_sample = x_sample; 
  foo.F = &PiPiStarForm; 
  foo.fname = std::string("ErrorEstimate.txt"); 

  return foo; 
}

START_PARS initialize_RhoPi(void)
{
  START_PARS foo;
  int ncfg = 1500; 
  itpp::Vec<double> mean(3); 
  itpp::Vec<double> error(3); 
  itpp::Mat<double> correlation(3,3); 
  itpp::Vec<double> x_sample = sample_range(-0.01,0.055,0.0001); 

  mean(0) = 0.492124; 
  error(0) = 0.005871; 
  mean(1) = 0.0511643; 
  error(1) = 0.001057; 
  mean(2) = -4.26168; 
  error(2) = 0.5463; 
  correlation(0,0) = 1.; 
  correlation(1,0) = -0.900; 
  correlation(2,0) = -0.757; 
  correlation(0,1) = -0.900; 
  correlation(1,1) = 1.; 
  correlation(2,1) = 0.951; 
  correlation(0,2) = -0.757; 
  correlation(1,2) = 0.951; 
  correlation(2,2) = 1.; 

  itpp::Mat<double> covariance = construct_covariance(correlation,error); 

  foo.ncfg = ncfg; 
  foo.mean = mean; 
  foo.covariance = covariance; 
  foo.x_sample = x_sample; 
  foo.F = &RhoPiForm; 
  foo.fname = std::string("ErrorEstimate.txt"); 

  return foo; 
}

START_PARS initialize_Rho1Pi(void)
{
  START_PARS foo;
  int ncfg = 500; 
  itpp::Vec<double> mean(3); 
  itpp::Vec<double> error(3); 
  itpp::Mat<double> correlation(3,3); 
  itpp::Vec<double> x_sample = sample_range(-0.01,0.04,0.0001); 

  mean(0) = -0.0518954; 
  error(0) = 0.0004874; 
  mean(1) = -1.88021; 
  error(1) = 0.1086; 
  mean(2) = 23.1301; 
  error(2) = 2.841; 
  correlation(0,0) = 1.; 
  correlation(1,0) = -0.653; 
  correlation(2,0) = 0.565; 
  correlation(0,1) = correlation(1,0); 
  correlation(1,1) = 1.; 
  correlation(2,1) = -0.988; 
  correlation(0,2) = correlation(2,0); 
  correlation(1,2) = correlation(2,1); 
  correlation(2,2) = 1.; 

  itpp::Mat<double> covariance = construct_covariance(correlation,error); 

  foo.ncfg = ncfg; 
  foo.mean = mean; 
  foo.covariance = covariance; 
  foo.x_sample = x_sample; 
  foo.F = &Polynomial2;
  foo.fname = std::string("ErrorEstimate.txt"); 

  return foo; 
}

START_PARS initialize_Rho2Pi(void)
{
  START_PARS foo;
  int ncfg = 500; 
  itpp::Vec<double> mean(3); 
  itpp::Vec<double> error(3); 
  itpp::Mat<double> correlation(3,3); 
  itpp::Vec<double> x_sample = sample_range(-0.042,0.03,0.0001); 

  mean(0) = 0.0230801; 
  error(0) = 0.002285; 
  mean(1) = -1.26713; 
  error(1) = 0.05989; 
  mean(2) = 36.4194; 
  error(2) = 3.478; 
  correlation(0,0) = 1.; 
  correlation(1,0) = 0.121; 
  correlation(2,0) = -0.814; 
  correlation(0,1) = correlation(1,0); 
  correlation(1,1) = 1.; 
  correlation(2,1) = -0.540; 
  correlation(0,2) = correlation(2,0); 
  correlation(1,2) = correlation(2,1); 
  correlation(2,2) = 1.; 

  itpp::Mat<double> covariance = construct_covariance(correlation,error); 

  foo.ncfg = ncfg; 
  foo.mean = mean; 
  foo.covariance = covariance; 
  foo.x_sample = x_sample; 
  foo.F = &Polynomial2;
  foo.fname = std::string("ErrorEstimate.txt"); 

  return foo; 
}

////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////


// void init_rho_mass()
// {
//   rho_masses.resize(3); 
//   rho_masses[0] = 0.216292;
//   rho_masses[1] = 0.396776; 
//   rho_masses[2] = 0.419859; 
//   rho_masses_init = true; 
// }

double RhoPiWidth(const double x, const itpp::Vec<double> p)
{
  double mRho = 0.216292; 
  double mpi = 0.1483; 
  double a = 1./137.; 
  double q = (pow(mRho,2.) - pow(mpi,2.))/(2.*mRho);   
  double prefactor = (a*4.*pow(q,3.))/(3*pow(mRho+mpi,2.)); 
  double conversion = 702./mpi;

  // put this into MeV
  return conversion*prefactor*pow(RhoPiForm(x,p),2.); 
}

double RhoPiWidthPhysical(const double x, const itpp::Vec<double> p)
{
  double mRho = 770; 
  double mpi = 135; 
  double a = 1./137.; 
  double q = (pow(mRho,2.) - pow(mpi,2.))/(2.*mRho);   
  double prefactor = (a*4.*pow(q,3.))/(3*pow(mRho+mpi,2.)); 
  double conversion = 1.;

  // put this into MeV
  return conversion*prefactor*pow(RhoPiForm(x,p),2.); 
}

double KstarKWidthPhysical(const double x, const itpp::Vec<double> p)
{
  double mRho = 892; 
  double mpi = 494; 
  double a = 1./137.; 
  double q = (pow(mRho,2.) - pow(mpi,2.))/(2.*mRho);   
  double prefactor = (a*4.*pow(q,3.))/(3*pow(mRho+mpi,2.)); 
  double conversion = 1.;

  // put this into MeV
  return conversion*prefactor*pow(RhoPiForm(x,p),2.); 
}

double Rho1PiWidth(const double x, const itpp::Vec<double> p)
{
  double mRho = 0.396776; 
  double mpi = 0.1483; 
  double a = 1./137.; 
  double q = (pow(mRho,2.) - pow(mpi,2.))/(2.*mRho);   
  double prefactor = (a*4.*pow(q,3.))/(3*pow(mRho+mpi,2.)); 
  double conversion = (702./mpi)*1000.;

  // put this into keV
  return conversion*prefactor*pow(Polynomial2(x,p),2.); 
}

double Rho2PiWidth(const double x, const itpp::Vec<double> p)
{
  double mRho = 0.419859; 
  double mpi = 0.1483; 
  double a = 1./137.; 
  double q = (pow(mRho,2.) - pow(mpi,2.))/(2.*mRho);   
  double prefactor = (a*4.*pow(q,3.))/(3*pow(mRho+mpi,2.)); 
  double conversion = (702./mpi)*1000.;

  // put this into keV
  return conversion*prefactor*pow(Polynomial2(x,p),2.); 
}

START_PARS rho_pi_width()
{
  START_PARS s =  initialize_RhoPi();
  itpp::Vec<double> x_sample(1);  
  x_sample.zeros(); 

  s.x_sample = x_sample; 

  s.F = &RhoPiWidth; 
  s.fname = std::string("WidthEstimateMeV.txt"); 

  return s; 
}

START_PARS rho_pi_width_phys()
{
  START_PARS s =  initialize_RhoPi();
  itpp::Vec<double> x_sample(1);  
  x_sample.zeros(); 

  s.x_sample = x_sample; 

  s.F = &RhoPiWidthPhysical; 
  s.fname = std::string("WidthEstimateMeV.phys.txt"); 

  return s; 
}

START_PARS kstar_k_width_phys()
{
  START_PARS s =  initialize_RhoPi();
  itpp::Vec<double> x_sample(1);  
  x_sample.zeros(); 

  s.x_sample = x_sample; 

  s.F = &KstarKWidthPhysical; 
  s.fname = std::string("WidthEstimateMeV.Kstarphys.txt"); 

  return s; 
}

START_PARS rho1_pi_width()
{
  START_PARS s =  initialize_Rho1Pi();
  itpp::Vec<double> x_sample(1);  
  x_sample.zeros(); 

  s.x_sample = x_sample; 

  s.F = &Rho1PiWidth; 
  s.fname = std::string("WidthEstimatekeV.txt"); 

  return s; 
}

START_PARS rho2_pi_width()
{
  START_PARS s =  initialize_Rho2Pi();
  itpp::Vec<double> x_sample(1);  
  x_sample.zeros(); 

  s.x_sample = x_sample; 

  s.F = &Rho2PiWidth; 
  s.fname = std::string("WidthEstimatekeV.txt"); 

  return s; 
}

////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////


void check_covariance(std::vector<itpp::Vec<double> > &draws, 
    const itpp::Vec<double> &mean,
    const itpp::Mat<double> cov)
{
  SEMBLE::SembleVector<double> foo(draws); 
  itpp::Vec<double> meand = foo.mean(); 
  std::cout << "mean(target) =  " << mean 
    << "\nmean(draws) = " << meand << std::endl;

  itpp::Mat<double> covd(cov); 
  covd.zeros(); 
  
  int ncfg = draws.size(); 
  int dim = mean.length(); 

  for(int i = 0; i < dim; ++i)
    for(int j = 0; j < dim; ++j)
      for(int bin = 0; bin < ncfg; ++bin)
        covd(i,j) += (foo[bin](i) - meand(i))*(foo[bin](j) -meand(i)); 

  covd /= double(ncfg); 

  std::cout << "cov(target) = " << cov
    << "\ncov(draws) = " << covd << std::endl; 

}
    

void check_generator()
{
  itpp::Normal_RNG rng(0.,1.); 
  int ndraw  = 500; 
  double mean(0.), variance(0.), work;
  std::vector<double> draws(ndraw,0); 
  std::vector<double>::const_iterator it; 

  std::string f("generator_histo"); 
  std::ofstream out(f.c_str()); 

  for(int i = 0; i < ndraw; ++i)
  {
    work = rng(); 
    draws[i] = work; 
    mean += work; 
    out << work << std::endl;
  }

  mean /= double(ndraw); 

  for(it = draws.begin(); it != draws.end(); ++it)
    variance += (*it - mean)*(*it - mean); 

  out.close(); 

  std::cout << "ndraw = " << ndraw 
    <<"\nmean = " << mean
    <<"\nvar  = " << variance/double(ndraw)
    << std::endl;
}

void run_code(const START_PARS &thingy)
{
  // check_generator(); 

  // draw an ensemble 
  int dim = thingy.mean.length(); 
  int ncfg = thingy.ncfg; 
  itpp::Mat<double> L = itpp::hermitian_transpose(itpp::chol(thingy.covariance)); 

  std::cout << "target covariance " << thingy.covariance << std::endl;
  std::cout << "chol " << L << std::endl;
  std::cout << "LLdag " << L * itpp::hermitian_transpose(L) << std::endl;

  std::vector<itpp::Vec<double> > draws(ncfg,thingy.mean); 
  std::vector<itpp::Vec<double> >::iterator it; 
  itpp::Normal_RNG rng(0.,1.); 
  itpp::Vec<double> draw(dim); 


  for(it = draws.begin(); it != draws.end(); ++it)
  {
    for(int i = 0; i < dim; ++i)
      draw(i) = rng(); 

    *it += L*draw; 
  }

  check_covariance(draws,thingy.mean, thingy.covariance);

  SEMBLE::SembleVector<double> ensemble(draws); 
  ensemble.rescaleSembleDown(); 

  // Function 
  double (*F)(const double x, const itpp::Vec<double> pars) = thingy.F;

  std::map<double,ENSEM::EnsemReal> result; 
  int nsamp = thingy.x_sample.length(); 
  for(int i = 0; i < nsamp; ++i)
  {
    ENSEM::EnsemReal ensem; 
    ensem.resize(ncfg); 
    double x = thingy.x_sample(i); 

    for(int b = 0 ; b < ncfg; ++b)
      ENSEM::pokeEnsem(ensem, ENSEM::Real( (*F)(x,ensemble[b]) ) , b); 


    result.insert(std::make_pair(x,ENSEM::rescaleEnsemUp(ensem))); 
  }

  // now dump it 

  std::ofstream out(thingy.fname.c_str()); 
  std::map<double,ENSEM::EnsemReal>::const_iterator rit; 
  for(rit = result.begin(); rit != result.end(); ++rit)
    out << rit->first << " " << ENSEM::toDouble(ENSEM::mean(rit->second))
      << " " << sqrt(ENSEM::toDouble(ENSEM::variance(rit->second))) << std::endl;

  out.close(); 
}


int main(void)
{
//  START_PARS Gm = initialize_Gm(); 
//  START_PARS Gq = initialize_Gq(); 
//  run_code(Gm); 
//  run_code(Gq); 
  
//  START_PARS s = initialize_PiStar_short(); 
//  START_PARS l = initialize_PiStar_long(); 
//  run_code(s); 
//  run_code(l); 

//  START_PARS s = initialize_Pi_short(); 
//  run_code(s); 

  START_PARS s = initialize_PiPiStar(); 
  run_code(s); 

//  START_PARS s =  initialize_RhoPi();
//  START_PARS sw = rho_pi_width();
//  START_PARS swp = rho_pi_width_phys();
//  START_PARS swk = kstar_k_width_phys();
//  run_code(s); 
//  run_code(sw); 
//  run_code(swp); 
//  run_code(swk); 

//  START_PARS s =  initialize_Rho1Pi();  
//  START_PARS sw = rho1_pi_width();
//
//  run_code(s); 
//  run_code(sw); 

//  START_PARS s =  initialize_Rho2Pi();  
//  START_PARS sw = rho2_pi_width();
//
//  run_code(s); 
//  run_code(sw); 

  return 0; 
}


