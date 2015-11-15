/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : grab_pol.cc

 * Purpose :

 * Creation Date : 22-04-2013

 * Last Modified : Thu 12 Dec 2013 10:29:40 PM EST

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/



#include "radmat/register_all/register_all.h"
#include "radmat/ff/lorentzff_polarization_embedding.h"
#include "radmat/ff/lorentzff_polarization_embedding_old.h"
#include "radmat/ff/lorentzff_formfac_utils.h"
#include <iostream>
#include <sstream>

using namespace radmat;




int local_flag = 0; 


template<int J, int H>
void print_tensor(const Tensor<double,1> &l, 
    const Tensor<double,1> &r, 
    const double &mom_kick)
{

  if(local_flag == 1)
  {
    std::cout << "USING NEW PTENS" << std::endl;

    radmat::embedHelicityPolarizationTensor<J,H> z;
    radmat::leftPTensor<J,H> lefty;
    radmat::rightPTensor<J,H> righty;

    std::cout << "lefty: " << l << "\nrighty:" << r << std::endl;
    std::cout << "\nz-axis left: " << z.z_axis_helicity_tensor(l,mom_kick);
    std::cout << "\nz-axis right: " << z.z_axis_helicity_tensor(r,mom_kick);
    std::cout << "\nlefty: " << lefty.left_p_tensor(l,r,mom_kick); 
    std::cout << "\nrighty: " << righty.right_p_tensor(l,r,mom_kick); 


    std::cout << "\ncanonical_frame: "
      << radmat::LatticeRotationEnv::rotation_group_label(
          radmat::get_space_mom(l,mom_kick),radmat::get_space_mom(r,mom_kick)) 
      << std::endl;
  }
  else
  {
    std::cout << "USING OLD PTENS" << std::endl;

    radmat::embedHelicityPolarizationTensor_old<J,H> z;
    radmat::leftPTensor_old<J,H> lefty;
    radmat::rightPTensor_old<J,H> righty;

    std::cout << "lefty: " << l << "\nrighty:" << r << std::endl;
    std::cout << "\nz-axis left: " << z.z_axis_helicity_tensor(l,mom_kick);
    std::cout << "\nz-axis right: " << z.z_axis_helicity_tensor(r,mom_kick);
    std::cout << "\nlefty: " << lefty.left_p_tensor(l,r,mom_kick); 
    std::cout << "\nrighty: " << righty.right_p_tensor(l,r,mom_kick); 


    std::cout << "\ncanonical_frame: "
      << radmat::LatticeRotationEnv::rotation_group_label(
          radmat::get_space_mom(l,mom_kick),radmat::get_space_mom(r,mom_kick)) 
      << std::endl;
  }
}



void print_tensor_1(const Tensor<double,1> &l, 
    const Tensor<double,1> &r, 
    const double &mom_kick,
    const short hel)
{
#define J 1
  switch(hel)
  {
    case 1:
      print_tensor<J,1>(l,r,mom_kick);
      break;
    case 0:
      print_tensor<J,0>(l,r,mom_kick);
      break;
    case -1:
      print_tensor<J,-1>(l,r,mom_kick);
      break;
    default:
      std::cout << "Hel " << hel << "is ill-defined" << std::endl;
  }

}





void print_tensor_2(const Tensor<double,1> &l, 
    const Tensor<double,1> &r, 
    const double &mom_kick,
    const short hel)
{
#undef J
#define J 2

  switch(hel)
  {
    case -2:
      print_tensor<J,2>(l,r,mom_kick);
      break;
    case 2:
      print_tensor<J,2>(l,r,mom_kick);
      break;
    case 1:
      print_tensor<J,1>(l,r,mom_kick);
      break;
    case 0:
      print_tensor<J,0>(l,r,mom_kick);
      break;
    case -1:
      print_tensor<J,-1>(l,r,mom_kick);
      break;
    default:
      std::cout << "Hel " << hel << "is ill-defined" << std::endl;
  }

}


void print_tensor_3(const Tensor<double,1> &l, 
    const Tensor<double,1> &r, 
    const double &mom_kick,
    const short hel)
{
#undef J
#define J 3
  switch(hel)
  {
    case -3:
      print_tensor<J,3>(l,r,mom_kick);
      break;
    case 3:
      print_tensor<J,3>(l,r,mom_kick);
      break;
    case -2:
      print_tensor<J,2>(l,r,mom_kick);
      break;
    case 2:
      print_tensor<J,2>(l,r,mom_kick);
      break;
    case 1:
      print_tensor<J,1>(l,r,mom_kick);
      break;
    case 0:
      print_tensor<J,0>(l,r,mom_kick);
      break;
    case -1:
      print_tensor<J,-1>(l,r,mom_kick);
      break;
    default:
      std::cout << "Hel " << hel << "is ill-defined" << std::endl;
  }

}


#undef J
void print_tensor(const Tensor<double,1> &l, 
    const Tensor<double,1> &r, 
    const double &mom_kick, 
    const short J,
    const short hel)
{

  switch(J)
  {
    case 1:
      print_tensor_1(l,r,mom_kick,hel);
      break;
    case 2:
      print_tensor_2(l,r,mom_kick,hel); 
      break;
    case 3:
      print_tensor_3(l,r,mom_kick,hel); 
      break;
    default: 
      std::cout << "J = " << J << "is not supported " << std::endl;
  }

}


  Tensor<double,1>
gen_t(const double m, const int x, const int y, const int z, const double kick)
{
  Tensor<double,1> r( (TensorShape<1>())[4] ,0.);
  r[0] = sqrt( m*m + kick*kick*( x*x + y*y + z*z) );
  r[1] = kick *x;
  r[2] = kick *y;
  r[3] = kick *z;
  return r; 
}


template<typename T>
  void
read_par(const int i, T &t, char *argv[])
{
  std::istringstream val(argv[i]);
  val >> t; 
}


int main(int argc, char *argv[])
{
  AllFactoryEnv::registerAll(); 

  if( argc != 13 )
  {
    std::cout << "error: usage: <grab_pol> " 
      << "m1 m2 p1 p2 p3 pp1 pp2 pp3 mom_kick J H (1->new)"
      << std::endl;
    exit(1);
  }

  double m1,m2; 
  int p1,p2,p3;
  int p4,p5,p6;  
  double kick; 
  short J,H;

  read_par(1,m1,argv);
  read_par(2,m2,argv);
  read_par(3,p1,argv);
  read_par(4,p2,argv);
  read_par(5,p3,argv);
  read_par(6,p4,argv);
  read_par(7,p5,argv);
  read_par(8,p6,argv);
  read_par(9,kick,argv);
  read_par(10,J,argv);
  read_par(11,H,argv);
  read_par(12,local_flag,argv);

  print_tensor(
      gen_t(m1,p1,p2,p3,kick),
      gen_t(m2,p4,p5,p6,kick),
      kick,J,H);

  return 0; 
}


