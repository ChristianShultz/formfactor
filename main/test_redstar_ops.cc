/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : test_redstar__ops.cc

 * Purpose :

 * Creation Date : 8-122013

 * Last Modified : Mon 28 Apr 2014 04:26:30 PM EDT

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/


#include "hadron/twoquark_dirac_ops.h"
#include "hadron/twoquark_displace.h"
#include "hadron/twoquark_cont_rest_ops.h"
#include "hadron/twoquark_cont_helicity_ops.h"
#include "ensem/ensem.h"
#include "io/adat_xmlio.h"
#include "semble/semble_meta.h"
#include "radmat/register_all/register_all.h"
#include "radmat/redstar_interface/redstar_interface.h"
#include "radmat/ff/ff.h"
#include "itpp/itbase.h"
#include <map>
#include <algorithm>
#include <string>
#include <sstream>
#include <iostream>
#include <utility>


// a hacky thingy
struct KeyRetriever 
{
  template<typename T>
    typename T::first_type operator() (T KeyValPair) const
    {
      return KeyValPair.first; 
    }
};

// what are we trying to do 
  std::pair<std::string,int>
usage(int argc, char *argv[] ) 
{
  if ( argc != 3 ) 
  {
    std::cerr << "error usage: print_twoquark_ops <pattern> <row> " << std::endl; 
    exit(1337); 
  }

  std::istringstream val(argv[1]), val2(argv[2]);
  std::string pattern; 
  int row; 
  val >> pattern; 
  val2 >> row; 
  std::cout << "Matching on " << pattern 
    << " row " << row << std::endl ; 


  return std::pair<std::string,int>(pattern,row);
}

  itpp::Mat<std::complex<double> > 
convertSpinMatrix(const ENSEM::SpinMatrix &m)
{
  itpp::Mat<std::complex<double> > spin_mat(4,4);
  spin_mat.zeros(); 

  for(int i =0; i < 4; ++i)
    for(int j = 0; j<4; ++j) 
      spin_mat(i,j) = SEMBLE::toScalar(ENSEM::peekSpin(m,i,j)); 

  return itpp::round_to_zero(spin_mat,1e-6); 
}


void printSpinMatrix(const ENSEM::SpinMatrix &m, std::complex<double> w=std::complex<double>(1.,0))
{
  itpp::Mat<std::complex<double> > spin_matrix = convertSpinMatrix(m); 
  std::cout << itpp::round_to_zero(spin_matrix, 1e-6)  << "\n" << std::endl; 
}
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

// do work 
void print_spin_keys(int argc, char *argv[])
{
  if( argc != 2 ) 
  {
    std::cerr << "error usage: test_redstar ops ["<< __func__ << "] " 
      << std::endl;
    exit(1337);
  }

  // the junk in adat/lib/hadron
  typedef Hadron::TheTwoQuarkDiracSpinOpsRegInfoFactory RegMap; 
  typedef Hadron::TheTwoQuarkDiracSpinOpsFactory TwoQuarkFoundry; 

  // register ops
  Hadron::TwoQuarkDiracSpinOpsEnv::registerAll(); 

  // keys and iterator
  std::vector<std::string> reg_keys; 
  std::vector<std::string>::const_iterator  reg_keys_it; 

  // pull the keys out 
  std::transform(RegMap::Instance().begin(),
      RegMap::Instance().end(),
      std::back_inserter(reg_keys),
      KeyRetriever()); 

  // print them
  for(reg_keys_it = reg_keys.begin(); reg_keys_it != reg_keys.end(); ++reg_keys_it)
    std::cout <<  *reg_keys_it << std::endl; 
}

///////////////////////////////////////////////
///////////////////////////////////////////////
///////////////////////////////////////////////


  itpp::Mat<std::complex<double> > 
gamma_matrix(const int i )
{
  Hadron::MapTwoQuarkExpr_t expr = ENSEM::GammaDP(1 << i -1) * Hadron::initOne(); 
  std::vector<Hadron::KeyTwoQuarkExpr_t> expr_keys = expr.keys(); 
  std::vector<Hadron::KeyTwoQuarkExpr_t>::const_iterator k;

  itpp::Mat<std::complex<double> > spin_mat(4,4); 
  spin_mat.zeros(); 

  for( k = expr_keys.begin(); k != expr_keys.end(); ++k)
    spin_mat += convertSpinMatrix(expr[*k]); 

  return spin_mat; 
}



///////////////////////////////////////////////
///////////////////////////////////////////////
///////////////////////////////////////////////

void print_gamma_matrices (int argc , char *argv[] )
{
  if ( argc != 2 ) 
  {
    std::cerr << "error usage: test_redstar ops ["<< __func__ << "] " 
      << std::endl; 
    exit(1337); 
  }

  // now make the dictionary 
  for(int i = 0; i < 4; ++i) 
    std::cout << "GammaDP(" << i << "):\n"
                                        << gamma_matrix(i+1) << std::endl; 
}



///////////////////////////////////////////////
///////////////////////////////////////////////
///////////////////////////////////////////////
// do work 
void print_spin_ops(int argc, char *argv[])
{
  if( argc != 4 ) 
  {
    std::cerr << "error usage: test_redstar ops ["<< __func__ << "] " 
      << " <pattern> <row> " << std::endl; 
    exit(1337);
  }

  // the junk in adat/lib/hadron
  typedef Hadron::TheTwoQuarkDiracSpinOpsRegInfoFactory RegMap; 
  typedef Hadron::TheTwoQuarkDiracSpinOpsFactory TwoQuarkFoundry; 

  // get pattern from cmd line  
  std::istringstream val(argv[2]) , val2(argv[3]); 
  std::string pattern;
  int row; 

  val >> pattern; 
  val2 >> row; 

  // register ops
  Hadron::TwoQuarkDiracSpinOpsEnv::registerAll(); 

  // keys and iterator
  std::vector<std::string> reg_keys; 
  std::vector<std::string>::const_iterator  reg_keys_it; 

  // pull the keys out 
  std::transform(RegMap::Instance().begin(),
      RegMap::Instance().end(),
      std::back_inserter(reg_keys),
      KeyRetriever()); 

  // print them
  for(reg_keys_it = reg_keys.begin(); reg_keys_it != reg_keys.end(); ++reg_keys_it)
    if( *reg_keys_it ==  pattern)
    {
      Hadron::TwoQuarkDiracSpinOpsRegInfo_t reg = RegMap::Instance().find(*reg_keys_it)->second;
      Hadron::TwoQuarkDiracSpinOps* op = TwoQuarkFoundry::Instance().createObject(reg.id); 

      Hadron::MapTwoQuarkExpr_t expr =  (*op)(Hadron::initOne(),row); 

      std::vector<Hadron::KeyTwoQuarkExpr_t> expr_keys = expr.keys(); 
      std::vector<Hadron::KeyTwoQuarkExpr_t>::const_iterator k;

      for( k = expr_keys.begin(); k != expr_keys.end(); ++k)
        printSpinMatrix(expr[*k]);
    }

}

////////////////////////////////////////////////////
////////////////////////////////////////////////////
////////////////////////////////////////////////////



void print_helicity_keys(int argc, char *argv[])
{
  if ( argc != 2 ) 
  {
    std::cerr << "error usage: test_redstar ops ["<< __func__ << "] " 
      << std::endl; 
    exit(1337); 
  }

  // the junk in adat/lib/hadron
  typedef Hadron::TheTwoQuarkContRestOpsRegInfoFactory RegMap; 


  // register ops
  Hadron::TwoQuarkContHelicityOpsEnv::registerAll(); 

  // keys and iterator
  std::vector<std::string> reg_keys; 
  std::vector<std::string>::const_iterator  reg_keys_it; 

  // pull the keys out 
  std::transform(RegMap::Instance().begin(),
      RegMap::Instance().end(),
      std::back_inserter(reg_keys),
      KeyRetriever()); 

  // print them
  for(reg_keys_it = reg_keys.begin(); reg_keys_it != reg_keys.end(); ++reg_keys_it)
  {
    std::cout << *reg_keys_it << std::endl; 
  }
}


////////////////////////////////////////////////////
////////////////////////////////////////////////////
////////////////////////////////////////////////////

  itpp::Mat<std::complex<double> > 
redstar_helicity_op(const std::string pattern, 
    const ADATXML::Array<int> &mom, 
    const int row) 
{
  // the junk in adat/lib/hadron
  typedef Hadron::TheTwoQuarkContRestOpsRegInfoFactory RegMap; 
  typedef Hadron::TheTwoQuarkContHelicityOpsFactory TwoQuarkHFoundry; 
  typedef Hadron::TheTwoQuarkContRestOpsFactory TwoQuarkRFoundry; 


  // register ops
  Hadron::TwoQuarkContRestOpsEnv::registerAll();
  Hadron::TwoQuarkContHelicityOpsEnv::registerAll(); 


  // keys and iterator
  std::vector<std::string> reg_keys; 
  std::vector<std::string>::const_iterator  reg_keys_it; 

  // pull the keys out 
  std::transform(RegMap::Instance().begin(),
      RegMap::Instance().end(),
      std::back_inserter(reg_keys),
      KeyRetriever()); 

  itpp::Mat<std::complex<double> > GAMMA(4,4); 
  GAMMA.zeros(); 

  int count(0); 

  // find the op
  for(reg_keys_it = reg_keys.begin(); reg_keys_it != reg_keys.end(); ++reg_keys_it)
  {
    if( *reg_keys_it == pattern)
    {
      Hadron::MapTwoQuarkExpr_t expr; 

      if( ( mom[0] == 0 ) && ( mom[1] == 0 ) && ( mom[2] == 0) ) 
      {
        Hadron::TwoQuarkContRestOps *op; 
        op = TwoQuarkRFoundry::Instance().createObject(pattern,pattern); 
        expr = (*op)(row); 
      }
      else
      {
        Hadron::TwoQuarkContHelicityOps *op; 
        op = TwoQuarkHFoundry::Instance().createObject(pattern,pattern); 
        expr = (*op)(mom,row); 
      }

      std::vector<Hadron::KeyTwoQuarkExpr_t> expr_keys = expr.keys(); 
      std::vector<Hadron::KeyTwoQuarkExpr_t>::const_iterator k;
      for( k = expr_keys.begin(); k != expr_keys.end(); ++k)
        GAMMA += convertSpinMatrix( expr[*k] );

      ++count; 
    }
  }

  if (count != 1) 
  {
    std::cout << __func__ << ": count = " << count 
      << " returning zero matrix" << std::endl;
    GAMMA.zeros(); 
  }

  return GAMMA; 
}

////////////////////////////////////////////////////
////////////////////////////////////////////////////
////////////////////////////////////////////////////


void print_helicity_ops(int argc, char *argv[])
{
  if ( argc != 7 ) 
  {
    std::cerr << "error usage: test_redstar ops ["<< __func__ << "] " 
      << " <op> <mom> <row> " << std::endl; 
    exit(1337); 
  }

  std::string pattern;
  int row;
  ADATXML::Array<int> mom; 

  std::istringstream val(argv[2]) , mx(argv[3]), my(argv[4]);
  std::istringstream mz(argv[5]) , r(argv[6]);

  val >> pattern;
  mom.resize(3);
  mx >> mom[0]; 
  my >> mom[1];
  my >> mom[2]; 
  r >> row; 

  std::cout << itpp::round_to_zero(
      redstar_helicity_op(pattern,mom,row) , 1e-6)  <<  std::endl; 
}


////////////////////////////////////////////////////
////////////////////////////////////////////////////
////////////////////////////////////////////////////

  itpp::Mat<std::complex<double> > 
radmat_spin_1_mat( const ADATXML::Array<int> &mom, const int row)
{
  itpp::Mat<std::complex<double> > eps = radmat::eps3d(mom,true); 
  itpp::Mat<std::complex<double> > GAMMA(4,4); 
  GAMMA.zeros(); 

  for(int i = 0; i < 3; ++i)
  {
    Hadron::MapTwoQuarkExpr_t expr = ENSEM::GammaDP(1 << i) * Hadron::initOne(); 
    std::vector<Hadron::KeyTwoQuarkExpr_t> expr_keys = expr.keys(); 
    std::vector<Hadron::KeyTwoQuarkExpr_t>::const_iterator k;
    for( k = expr_keys.begin(); k != expr_keys.end(); ++k)
      GAMMA += convertSpinMatrix( expr[*k] ) * eps(row-1,i);  
  }

  return GAMMA;
}

////////////////////////////////////////////////////
////////////////////////////////////////////////////
////////////////////////////////////////////////////

void radmat_spin_1(int argc , char *argv[] ) 
{
  std::cout << __func__ << ": THIS IS A CREATION OP IN RADMAT " << std::endl;
  if ( argc !=  6)
  {
    std::cerr << "error usage: test_redstar ops ["<< __func__ << "] " 
      << "<mom> <row> " <<  std::endl;
    exit(1337);
  }

  int row;
  ADATXML::Array<int> mom; 

  std::istringstream  mx(argv[2]), my(argv[3]);
  std::istringstream mz(argv[4]) , r(argv[5]);

  mom.resize(3);
  mx >> mom[0]; 
  my >> mom[1];
  my >> mom[2]; 
  r >> row; 

  std::cout << itpp::round_to_zero( 
      radmat_spin_1_mat(mom,row), 1e-6) << std::endl;
}

////////////////////////////////////////////////////
////////////////////////////////////////////////////
////////////////////////////////////////////////////

void radmat_cartesianize_red_rho(int argc , char *argv[] ) 
{
  std::cout << __func__ << ": THIS IS A CREATION OP IN RADMAT " << std::endl;
  if ( argc !=  5)
  {
    std::cerr << "error usage: test_redstar ops ["<< __func__ << "] " 
      << "<mom>" <<  std::endl;
    exit(1337);
  }

  int row;
  ADATXML::Array<int> mom; 

  std::istringstream  mx(argv[2]), my(argv[3]), mz(argv[4]) ;

  mom.resize(3);
  mx >> mom[0]; 
  my >> mom[1];
  my >> mom[2]; 

  itpp::Vec<itpp::Mat<std::complex<double> > > red(3);  
  std::string op = "rhoxD0_J0__J1";

  for(int i = 1; i < 4; ++i)
    red[i-1] = redstar_helicity_op(op,mom,i); 

  itpp::Mat<std::complex<double> > inverse_tform;
  inverse_tform = itpp::hermitian_transpose(radmat::eps3d(mom,true)); 

  std::cout << __func__ << ": M_{j,lambda} M_{lambda,j}\n" 
    << itpp::round_to_zero( radmat::invert2Cart(mom,true) * radmat::eps3d(mom,true) , 1e-6) 
    << std::endl; 

  itpp::Mat<std::complex<double> > zero(4,4); 
  zero.zeros(); 

  // zero it out
  itpp::Vec<itpp::Mat<std::complex<double> > > cart(3); 
  for(int i = 0; i < 3; ++i)
    cart[i] = zero;

  for(int i =0; i < 3; ++i)
    for(int j =0; j <3; ++j)
      cart[i] += inverse_tform(i,j) * red[j]; 

  for(int i = 0; i < 3; ++i)
    std::cout << __func__ << ": cartesian " << i 
      << "\n" << itpp::round_to_zero(cart[i], 1e-6) << std::endl;
}

////////////////////////////////////////////////////
////////////////////////////////////////////////////
////////////////////////////////////////////////////

void radmat_test_photon_inversion(int argc , char *argv[] ) 
{
  std::cout << __func__ << ": THIS IS A CREATION OP IN RADMAT " << std::endl;
  if ( argc !=  2)
  {
    std::cerr << "error usage: test_redstar ops ["<< __func__ << "] " 
      <<  std::endl;
    exit(1337);
  }

  ADATXML::Array<int> mom; 
  mom.resize(3);
  itpp::Vec<itpp::Mat<std::complex<double> > > red(3);  
  std::string op = "rhoxD0_J0__J1";

  int nt(0),nb(0); // test count

  // zero
  itpp::Mat<std::complex<double> > zero(4,4); 
  zero.zeros(); 

  // expected answer
  itpp::Vec<itpp::Mat<std::complex<double> > > cart(3); 
  for(int i = 0; i < 3; ++i)
    cart[i] = gamma_matrix(i+1); // ensem offset

  for(int x = -2; x <= 2; ++x)
    for(int y = -2; y <= 2; ++y)
      for(int z = -2; z <= 2; ++z)
      {
        if( x*x + y*y + z*z > 4 ) 
          continue; 

        mom[0] = x; 
        mom[1] = y; 
        mom[2] = z;

        std::cout << __func__ << ": testing mom " << x << y << z << std::endl;

        // fill out the redstar
        for(int i = 1; i < 4; ++i)
          red[i-1] = redstar_helicity_op(op,mom,i); 

        // how do we move back to cartesian 
        itpp::Mat<std::complex<double> > inverse_tform;
        inverse_tform = radmat::invert2Cart(mom,true); 

        // zero it out
        itpp::Mat<std::complex<double> >  soln; 

        // the cartesian index we are testing 
        for(int i =0; i < 3; ++i)
        {
          soln = zero; 

          // do the inversion on this index
          for(int j =0; j <3; ++j)
            soln += inverse_tform(i,j) * red[j]; 

          // test this index
          if( itpp::round_to_zero(soln - cart[i],1e-6) != zero )
          {
            std::cout << __func__ << ": inversion error, mom=" 
              << x << y << z << " cart " <<  i 
              << "\nanswer:\n" << itpp::round_to_zero(soln,1e-6)
              << "\nexpected:\n" << itpp::round_to_zero(cart[i],1e-6)
              << "\ndifference:\n" << itpp::round_to_zero(soln - cart[i], 1e-6) << std::endl;
            ++nb;
          }
          ++nt; 
        }

      } // z mom

  // result 
  std::cout << nt -nb << " of " << nt 
    << " tests passed  (any fails is really really bad!!! )" << std::endl;
}


////////////////////////////////////////////////////
////////////////////////////////////////////////////
////////////////////////////////////////////////////

itpp::Mat<std::complex<double> > 
pad_out(const radmat::Tensor<std::complex<double> , 1> &l, 
    const radmat::Tensor<std::complex<double> , 1> &r)
{
  itpp::Mat<std::complex<double> > sum(4,4); 
  sum.zeros(); 

  for(int i =0; i < 4; ++i)
    for(int j = 0; j < 4; ++j)
      sum(i,j) = l[i]*r[j]; 
  return sum; 
}

void check_ortho(const radmat::Tensor<std::complex<double>,1> &eps, 
    const radmat::Tensor<double,1> &p, 
    const std::string &proj)
{
  std::complex<double>  dc(0.,0.);
  dc += eps[0] * p[0]; 
  dc -= eps[1] * p[1]; 
  dc -= eps[2] * p[2]; 
  dc -= eps[3] * p[3]; 

  if( std::norm(dc) > 1e-6 )
  {
    std::cout << __func__ <<": projection " << proj 
      << " -- polarization not orthogonal to momentum " 
      << "\neps:  " << eps << "\nmom:  " <<  p << std::endl;
  }
}

itpp::Mat<std::complex<double> > 
sum_out_polarization_lefty(const radmat::Tensor<double,1> &p, 
    const radmat::Tensor<double,1> &pp, 
    const double mom_kick)
{
  radmat::embedHelicityPolarizationTensor<1,1> foo; 

  radmat::leftPTensor<1,1> P; 
  radmat::leftPTensor<1,0> Z; 
  radmat::leftPTensor<1,-1> M; 

  itpp::Mat<std::complex<double> > sum(4,4); 
  sum.zeros(); 
  radmat::Tensor<std::complex<double> , 1> eps; 
  radmat::Tensor<std::complex<double> , 1> epsc; 

  eps = P.left_p_tensor(p,pp,mom_kick); 
  epsc = foo.conjugate(P.left_p_tensor(p,pp,mom_kick)); 
  sum += pad_out(eps,epsc); 
  check_ortho(eps,p,"left+");

  eps = Z.left_p_tensor(p,pp,mom_kick); 
  epsc = foo.conjugate(Z.left_p_tensor(p,pp,mom_kick)); 
  sum += pad_out(eps,epsc); 
  check_ortho(eps,p,"left0");

  eps = M.left_p_tensor(p,pp,mom_kick); 
  epsc = foo.conjugate(M.left_p_tensor(p,pp,mom_kick)); 
  sum += pad_out(eps,epsc); 
  check_ortho(eps,p,"left-");


  return sum; 
}

void
print_eps_lefty(const radmat::Tensor<double,1> &p, 
    const radmat::Tensor<double,1> &pp, 
    const double mom_kick)
{
  radmat::leftPTensor<1,1> P; 
  radmat::leftPTensor<1,0> Z; 
  radmat::leftPTensor<1,-1> M; 
  radmat::Tensor<std::complex<double> , 1> eps; 

  eps = P.left_p_tensor(p,pp,mom_kick); 
  std::cout << "left+ " << eps; 

  eps = Z.left_p_tensor(p,pp,mom_kick); 
  std::cout << "left0 " << eps; 

  eps = M.left_p_tensor(p,pp,mom_kick); 
  std::cout << "left- " << eps; 
}


itpp::Mat<std::complex<double> > 
sum_out_polarization_righty(const radmat::Tensor<double,1> &p, 
    const radmat::Tensor<double,1> &pp, 
    const double mom_kick)
{
  radmat::embedHelicityPolarizationTensor<1,1> foo; 

  radmat::rightPTensor<1,1> P; 
  radmat::rightPTensor<1,0> Z; 
  radmat::rightPTensor<1,-1> M; 

  itpp::Mat<std::complex<double> > sum(4,4); 
  sum.zeros(); 
  radmat::Tensor<std::complex<double> , 1> eps; 
  radmat::Tensor<std::complex<double> , 1> epsc; 

  eps = P.right_p_tensor(p,pp,mom_kick); 
  epsc = foo.conjugate(P.right_p_tensor(p,pp,mom_kick)); 
  sum += pad_out(eps,epsc); 
  check_ortho(eps,pp,"right+");

  eps = Z.right_p_tensor(p,pp,mom_kick); 
  epsc = foo.conjugate(Z.right_p_tensor(p,pp,mom_kick)); 
  sum += pad_out(eps,epsc); 
  check_ortho(eps,pp,"right0");

  eps = M.right_p_tensor(p,pp,mom_kick); 
  epsc = foo.conjugate(M.right_p_tensor(p,pp,mom_kick)); 
  sum += pad_out(eps,epsc); 
  check_ortho(eps,pp,"right-");

  return sum; 
}

void
print_eps_righty(const radmat::Tensor<double,1> &p, 
    const radmat::Tensor<double,1> &pp, 
    const double mom_kick)
{
  radmat::rightPTensor<1,1> P; 
  radmat::rightPTensor<1,0> Z; 
  radmat::rightPTensor<1,-1> M; 
  radmat::Tensor<std::complex<double> , 1> eps; 

  eps = P.right_p_tensor(p,pp,mom_kick); 
  std::cout << "right+ " << eps; 

  eps = Z.right_p_tensor(p,pp,mom_kick); 
  std::cout << "right0 " << eps; 

  eps = M.right_p_tensor(p,pp,mom_kick); 
  std::cout << "right- " << eps; 
}


//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////



void radmat_test_polarization_ortho(int argc , char *argv[] ) 
{
  if ( argc !=  2)
  {
    std::cerr << "error usage: test_redstar ops ["<< __func__ << "] " 
      <<  std::endl;
    exit(1337);
  }

  ADATXML::Array<int> mom,mom2; 
  mom.resize(3); 
  mom2.resize(3); 

  int nt(0),nb(0); // test count


  // expected answer
  itpp::Mat<std::complex<double> >  suml,sumr, g, zero(4,4); 
  radmat::Tensor<double,1> p4((radmat::TensorShape<1>())[4]),pp4;
  pp4 = p4; 
  zero.zeros();
  suml = zero; 
  sumr = zero; 
  g = zero; 

  g(0,0) = 1.; 
  g(1,1) = -1.; 
  g(2,2) = -1.; 
  g(3,3) = -1.; 

  // rho on 743 lattice
  double mom_kick = 0.11; 
  double mass = 0.216; 


  for(int x = -2; x <= 2; ++x)
    for(int y = -2; y <= 2; ++y)
      for(int z = -2; z <= 2; ++z)
      {
        if( x*x + y*y + z*z > 4 ) 
          continue; 

        mom[0] = x; 
        mom[1] = y; 
        mom[2] = z;

        double E = sqrt(mass * mass + mom_kick*mom_kick*(x*x + y*y + z*z));  
        p4[0] = E; 
        p4[1] = mom_kick*mom[0];
        p4[2] = mom_kick*mom[1];
        p4[3] = mom_kick*mom[2];

        suml = zero; 
        for(int mu =0; mu < 4; ++mu)
          for(int nu =0; nu < 4; ++nu)
            suml(mu,nu) = -g(mu,nu) + p4[mu]*p4[nu]/mass/mass; 

        for(int i = -2; i <= 2; ++i)
          for(int j = -2; j <= 2; ++j)
            for(int k = -2; k <= 2; ++k)
            {
              if( i*i + j*j + k*k > 4 ) 
                continue; 

              if( (x-i)*(x-i) + (y-j)*(y-j) + (z-k)*(z-k) > 4)
                continue; 


              mom2[0] = i; 
              mom2[1] = j; 
              mom2[2] = k;

              double E2 = sqrt(mass * mass + mom_kick*mom_kick*(i*i+ j*j + k*k));  
              pp4[0] = E2; 
              pp4[1] = mom_kick*mom2[0];
              pp4[2] = mom_kick*mom2[1];
              pp4[3] = mom_kick*mom2[2];



              sumr = zero; 
              for(int mu =0; mu < 4; ++mu)
                for(int nu =0; nu < 4; ++nu)
                  sumr(mu,nu) = -g(mu,nu) + pp4[mu]*pp4[nu]/mass/mass; 


              std::cout << __func__ << ": testing lefty" << x << y << z 
                << "   righty " << i << j << k << std::endl;

              itpp::Mat<std::complex<double> > res_left = sum_out_polarization_lefty(p4,pp4,mom_kick);  

              // test this index
              if( itpp::round_to_zero(suml - res_left,1e-6) != zero )
              {
                std::cout << __func__ << ": lefty error, mom=" << p4 
                  << "\nanswer:\n" << itpp::round_to_zero(res_left,1e-6)
                  << "\nexpected:\n" << itpp::round_to_zero(suml,1e-6)
                  << "\ndifference:\n" << itpp::round_to_zero(res_left - suml, 1e-6) << std::endl;
                ++nb;
               print_eps_lefty(p4,pp4,mom_kick);  
              }
              ++nt; 


              itpp::Mat<std::complex<double> > res_right = sum_out_polarization_righty(p4,pp4,mom_kick);  

              // test this index
              if( itpp::round_to_zero(sumr - res_right,1e-6) != zero )
              {
                std::cout << __func__ << ": righty error, mom=" << pp4  
                  << "\nanswer:\n" << itpp::round_to_zero(res_right,1e-6)
                  << "\nexpected:\n" << itpp::round_to_zero(sumr,1e-6)
                  << "\ndifference:\n" << itpp::round_to_zero(res_right - sumr, 1e-6) << std::endl;
                ++nb;
               print_eps_righty(p4,pp4,mom_kick);  
              }
              ++nt; 

              if(nb)
                exit(1); 

            }
      }



  // result 
  std::cout << nt -nb << " of " << nt 
    << " tests passed  (any fails is really really bad!!! )" << std::endl;
}


////////////////////////////////////////////////////
////////////////////////////////////////////////////
////////////////////////////////////////////////////

// only works for the rho spin operator
void diff_rad_red_rho_hel_ops(int argc, char *argv[])
{
  std::cout << __func__ << ": THIS IS A CREATION OP IN RADMAT " << std::endl;
  if ( argc != 7 ) 
  {
    std::cerr << "error usage: test_redstar ops ["<< __func__ << "] " 
      << " <op> <mom> <row> " << std::endl; 
    exit(1337); 
  }

  std::string pattern;
  int row;
  ADATXML::Array<int> mom; 

  std::istringstream val(argv[2]) , mx(argv[3]), my(argv[4]);
  std::istringstream mz(argv[5]) , r(argv[6]);

  val >> pattern;
  mom.resize(3);
  mx >> mom[0]; 
  my >> mom[1];
  my >> mom[2]; 
  r >> row; 


  itpp::Mat<std::complex<double> > red = redstar_helicity_op(pattern,mom,row); 
  itpp::Mat<std::complex<double> > rad = radmat_spin_1_mat(mom,row); 

  std::cout << "redstar:\n" << itpp::round_to_zero(red,1e-6)
    << "\nradmat:\n" << itpp::round_to_zero(rad,1e-6)
    << "\nred-rad:\n" << itpp::round_to_zero(red - rad, 1e-6) << std::endl;
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////

void loop_diff_rad_red_rho(int argc, char *argv[])
{
  std::cout << __func__ << ": THIS IS A CREATION OP IN RADMAT " << std::endl;
  if ( argc != 2 ) 
  {
    std::cerr << "error usage: test_redstar ops ["<< __func__ << "] " 
      << std::endl;
    exit(1337); 
  }

  ADATXML::Array<int> mom; 

  itpp::Mat<std::complex<double> > zero(4,4); 
  zero.zeros();
  mom.resize(3);

  itpp::Mat<std::complex<double> > red; 
  itpp::Mat<std::complex<double> > rad; 

  std::string op = "rhoxD0_J0__J1";
  int nt(0),nb(0); 

  for(int x = -2; x <= 2; ++x)
    for(int y = -2; y <= 2; ++y)
      for(int z = -2; z <= 2; ++z)
      {
        if( x*x + y*y + z*z > 4 ) 
          continue; 

        mom[0] = x; 
        mom[1] = y; 
        mom[2] = z;

        for(int row = 1; row <4; ++row)
        {
          std::cout << "mom " << x << y << z << " row " << row << std::endl;

          red = redstar_helicity_op(op,mom,row); 
          rad = radmat_spin_1_mat(mom,row); 

          if( itpp::round_to_zero(rad-red,1e-6) != zero )
          {
            std::cout << __func__ << ": rotation error, mom=" 
              << x << y << z << " row " << row
              << "\nredstar:\n" << itpp::round_to_zero(red,1e-6)
              << "\nradmat:\n" << itpp::round_to_zero(rad,1e-6)
              << "\nred-rad:\n" << itpp::round_to_zero(red - rad, 1e-6) << std::endl;
            ++nb;
          }
          ++nt; 
        } // row

      } // z mom

  std::cout << nt -nb << " of " << nt 
    << " tests passed  (any fails is really really bad!!! )" << std::endl;
}


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////


void unique_frames(int argc, char *argv[])
{
  if ( argc != 2 ) 
  {
    std::cerr << "error usage: test_redstar ops ["<< __func__ << "] " 
      << std::endl;
    exit(1337); 
  }


  std::vector<radmat::mom_pair_key> frames; 
  std::vector<radmat::mom_pair_key>::const_iterator it; 
  frames = radmat::LatticeRotationEnv::TheRotationGroupGenerator::Instance().canonical_frames;

  for(it = frames.begin(); it != frames.end(); ++it)
    std::cout << *it << std::endl;

}



//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////




void related_by_rotation(int argc, char *argv[])
{
  if ( argc != 8 ) 
  {
    std::cerr << "error usage: test_redstar ops ["<< __func__ << "] " 
      << "<mom1>  <mom2> "<< std::endl;
    exit(1337); 
  }

  std::istringstream v(argv[2]),vv(argv[3]),vvv(argv[4]);
  std::istringstream w(argv[5]),ww(argv[6]),www(argv[7]);

  ADATXML::Array<int> p1,p2; 
  p1.resize(3); 
  v >> p1[0]; 
  vv >> p1[1]; 
  vvv >> p1[2]; 
  p2.resize(3); 
  w >> p2[0];
  ww >> p2[1];
  www >> p2[2]; 


  std::vector<radmat::mom_pair_key> frames; 
  std::vector<radmat::mom_pair_key>::const_iterator frame_it; 
  std::map<radmat::mom_pair_key,radmat::mom_pair_key,radmat::mom_key_comp>::const_iterator it; 
  std::pair<radmat::mom_t,radmat::mom_t> canonical = radmat::LatticeRotationEnv::rotation_group_can_mom(p1,p2); 
  radmat::mom_pair_key canonical_key( radmat::mom_key(canonical.first) , radmat::mom_key(canonical.second)); 

  typedef radmat::LatticeRotationEnv::TheRotationGroupGenerator RG; 
  radmat::mom_key_comp comparator; 

  for( it = RG::Instance().can_frame_map.begin() ; it != RG::Instance().can_frame_map.end(); ++it)
    if( comparator.exact_equivalence( canonical_key, it->second) )
      frames.push_back( it->first ); 

  for(frame_it = frames.begin(); frame_it != frames.end(); ++frame_it)
    std::cout << *frame_it << std::endl;

}



//////////////////////////////////////////////////
//////////////////////////////////////////////////
//////////////////////////////////////////////////

typedef void (*fptr)(int argc , char *argv[]); 
std::map<std::string , fptr > options;

void insert_op(const std::string &s, const fptr &f)
{
  options.insert(std::pair<std::string,fptr>(s,f)); 
}

void init_options(void)
{
  insert_op("unique_frames",&unique_frames); 
  insert_op("related_by_rotation",&related_by_rotation); 
  insert_op("radmat_test_pol_ortho",&radmat_test_polarization_ortho); 
  insert_op("radmat_test_photon_inversion",&radmat_test_photon_inversion); 
  insert_op("cartesianize_red_rho",&radmat_cartesianize_red_rho); 
  insert_op("loop_diff_rad_red_rho",&loop_diff_rad_red_rho);
  insert_op("diff_rad_red_rho",&diff_rad_red_rho_hel_ops);
  insert_op("radmat_spin_1",&radmat_spin_1); 
  insert_op("print_gamma_matrices",&print_gamma_matrices); 
  insert_op("print_helicity_keys",&print_helicity_keys); 
  insert_op("print_helicity_ops",&print_helicity_ops); 
  insert_op("print_spin_keys",&print_spin_keys); 
  insert_op("print_spin_ops",&print_spin_ops); 
}


// pick appropriate function and pass on command line inputs 
void do_work(std::string &op, int argc,char *argv[])
{
  init_options(); 

  if(options.find(op) == options.end())
  {
    std::cerr << " unrecognized op " << op 
      << " more intelligent choices are " << std::endl; 
    std::map<std::string , fptr>::const_iterator it; 
    for(it = options.begin(); it != options.end(); ++it)
      std::cerr << it->first << std::endl; 
    exit(1); 
  }

  fptr foo = options[op];

  foo(argc,argv); 
}


// main program wrapper
int main(int argc, char *argv[])
{
  radmat::AllFactoryEnv::registerAll(); 

  // we will always have at least 2 , radmat_util operation_with_no_inputs
  if(argc < 2)
  {
    std::cerr << "usage: test_redstar_ops : <operation> <op inputs ...> " << std::endl;
    exit(1); 
  }

  std::string op;
  std::istringstream opi(argv[1]); 
  opi >> op; 

  do_work(op,argc,argv); 

  return 0;
}









