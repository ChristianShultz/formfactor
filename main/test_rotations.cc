/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

* File Name : test_rotations.cc

* Purpose :

* Creation Date : 11-12-2013

* Last Modified : Mon 28 Apr 2014 06:17:52 PM EDT

* Created By : shultz

_._._._._._._._._._._._._._._._._._._._._.*/

#include "radmat/register_all/register_all.h"
#include <map>
#include <string>
#include <sstream>
#include "radmat/redstar_interface/redstar_canonical_rotations.h"
#include "radmat/redstar_interface/redstar_canonical_lattice_rotations.h"
#include "radmat/ff/ff.h"
#include "radmat/utils/levi_civita.h"
#include "hadron/irrep_util.h"
#include "formfac/formfac_qsq.h"

typedef radmat::mom_t mom_t; 
typedef radmat::RotationMatrix_t RotationMatrix_t; 

// a hacky thingy
struct KeyRetriever 
{
  template<typename T>
    typename T::first_type operator() (T KeyValPair) const
    {
      return KeyValPair.first; 
    }
};

/////////////////////////////////////////////
radmat::Tensor<double,1> init_4tens(const double m, 
    int x, int y, int z, double kick)
{
  radmat::Tensor<double,1> foo( (radmat::TensorShape<1>())[4] , 0.); 

  foo[0] = sqrt( m*m + kick*kick*(x*x + y*y + z*z) ); 
  foo[1] = kick * x; 
  foo[2] = kick * y; 
  foo[3] = kick * z; 

  return foo; 
}

/////////////////////////////////////////////
std::string string_mom(const mom_t &p)
{
  std::stringstream ss; 
  ss << p[0] << p[1] << p[2];
  return ss.str(); 
}

/////////////////////////////////////////////
mom_t ref_mom(const mom_t &p)
{
  mom_t foo = FF::canonicalOrder(p); 

  if( (foo[0] == 1) && (foo[1] == 1) && (foo[2] == 1) )
  {
    return radmat::gen_mom<1,1,1>();
  }
  else if( (foo[0] == 1) && (foo[1] == 1) && (foo[2] == 0) )
  {
    return radmat::gen_mom<0,1,1>();
  }
  else if( (foo[0] == 1) && (foo[1] == 0) && (foo[2] == 0) )
  {
    return radmat::gen_mom<0,0,1>();
  }
  else if( (foo[0] == 2) && (foo[1] == 0) && (foo[2] == 0) )
  {
    return radmat::gen_mom<0,0,2>();
  }
  else if( (foo[0] == 0) && (foo[1] == 0) && (foo[2] == 0) )
  {
    return radmat::gen_mom<0,0,0>();
  }
  else
  {
    std::cout << __func__ << ": unknown momentum, can " 
      << string_mom(foo) << "  p " << string_mom(p) << std::endl;
  }

}


/////////////////////////////////////////////
  mom_t 
rotate_int_mom(const RotationMatrix_t *R,
    const mom_t &l)
{
  mom_t chk = radmat::gen_mom<0,0,0>();

  for(int i = 0; i < 3; ++i)
  {
    double res(0.); 
    for(int j = 0; j < 3; ++j)
      if ( fabs( (*R)[i+1][j+1] ) > 1e-6 )
        res += (*R)[i+1][j+1] * double(l[j]); 

    chk[i] = (res > 0.) ? int( fabs(res) + 0.5 ) : - int( fabs(res) + 0.5 ); 
  }
  return chk;
} 


/////////////////////////////////////////////
  mom_t 
rotate_z_mom(const RotationMatrix_t *R,
    const double mod_p)
{
  mom_t chk = radmat::gen_mom<0,0,0>();

  for(int i = 0; i < 3; ++i)
  {
    double res(0.); 

    if ( fabs( (*R)[i+1][3] ) > 1e-6 )
      res = (*R)[i+1][3] * mod_p; 

    chk[i] = (res > 0.) ? int( fabs(res) + 0.5 ) : - int( fabs(res) + 0.5 ); 
  }
  return chk;
} 




//////////////////////////////////////////////////
//////////////////////////////////////////////////
//////////////////////////////////////////////////

void check_rotations_mom(int argc, char *argv[])
{
  if( argc != 2 )
  {
    std::cout << "error usage: test_rotations [" << __func__ << "] "
      << std::endl; 
  }

  mom_t chk = radmat::gen_mom<0,0,0>();
  int nb(0), ng(0);

  for(int x = -2; x <=2; ++x)
    for(int y = -2; y <=2; ++y)
      for(int z = -2; z <=2; ++z)
      {
        int psq = x*x + y*y + z*z;

        if( psq > 4)
          continue; 

        chk[0] = x;
        chk[1] = y; 
        chk[2] = z; 

        std::string chk_str = string_mom(chk); 
        std::cout << "testing " << chk_str << std::endl;

        RotationMatrix_t *R = radmat::CanonicalRotationEnv::call_factory(chk); 
        mom_t pr = rotate_z_mom(R,sqrt(psq)); 

        if( chk_str == string_mom(pr)  )
        {
          ++ng; 
        }
        else
        {
          std::cout << __func__ << ": Rotation error for " 
            << chk_str << " got " << string_mom(pr) << std::endl;
          ++nb;
          ++ng; 
        }
      } 

  std::cout << ng - nb << " of " 
    << ng << " tests passed " << std::endl;

}


//////////////////////////////////////////////////
//////////////////////////////////////////////////
//////////////////////////////////////////////////

void check_lattice_rotations_mom(int argc, char *argv[])
{
  if( argc != 2 )
  {
    std::cout << "error usage: test_rotations [" << __func__ << "] "
      << std::endl; 
  }

  mom_t chk = radmat::gen_mom<0,0,0>();
  int nb(0), ng(0);

  for(int x = -2; x <=2; ++x)
    for(int y = -2; y <=2; ++y)
      for(int z = -2; z <=2; ++z)
      {
        int psq = x*x + y*y + z*z;

        if( psq > 4)
          continue; 

        chk[0] = x;
        chk[1] = y; 
        chk[2] = z; 

        std::string chk_str = string_mom(chk); 
        std::cout << "testing " << chk_str << std::endl;

        mom_t ref = ref_mom(FF::canonicalOrder(chk)); 
        std::string LG = Hadron::generateLittleGroup(ref); 
        RotationMatrix_t *Rref = radmat::CanonicalRotationEnv::call_factory(LG);  
        RotationMatrix_t *R = radmat::CanonicalLatticeRotationEnv::call_factory(chk); 
        mom_t pr = rotate_z_mom(Rref,sqrt(psq)); 

        if( string_mom(ref) == string_mom(pr) )
        {
          ++ng; 
        }
        else
        {
          std::cout << __func__ << ": Rotation error for " 
            << string_mom(ref) << ": got " 
           << string_mom(pr) << "\nR:\n" << *Rref << std::endl;
          ++nb;
          ++ng; 
        }

        mom_t p = rotate_int_mom(R,pr); 

        if( string_mom(p) == chk_str)
        {
          ++ng;
        }
        else
        {
          std::cout << __func__ << ": rotation error, tried to rotate " 
            << string_mom(pr) << " to " << chk_str << " got " 
            << string_mom(p) << "\nR:\n" << *R << std::endl;
          ++nb;
          ++ng; 
        }

        delete R;
        delete Rref; 
      } 

  std::cout << ng - nb << " of " 
    << ng << " tests passed " << std::endl;

}


//////////////////////////////////////////////////
//////////////////////////////////////////////////
//////////////////////////////////////////////////

void read_momentum(const int st, mom_t &p, char *argv[])
{
  for(int i = st; i < st+3; ++i)
  {
    std::istringstream val(argv[i]);
    val >> p[i-st];
  }
}



void get_lattice_rotation(int argc, char *argv[])
{
  if( argc != 8 )
  {
    std::cout << "error usage: test_rotations [" << __func__ << "] "
      << "<mom1> <mom2> " << std::endl; 
    exit(1);
  }
  
  mom_t m1 = radmat::gen_mom<0,0,0>(); 
  mom_t m2 = radmat::gen_mom<0,0,0>(); 

  read_momentum(2,m1,argv);
  read_momentum(5,m2,argv);

  std::cout << "read m1 " << string_mom(m1) 
    << " m2 " << string_mom(m2) << std::endl;
  std::cout << " returning R_lat * m1 = m2 " << std::endl;

   radmat::rCompEulAngles eul = radmat::generate_frame_transformation(m2,m1);
  RotationMatrix_t *r = new RotationMatrix_t( genRotationMatrix(eul) ); 

  if( !!! radmat::check_frame_transformation(r, m1 , m2 ) )
  {
    std::cout << __func__ 
      << ": WARNING -- lattice transformation failed " << std::endl;
  }

  std::cout << "R:" << *r << std::endl;

  std::cout << "det(R) = " << determinant(r) << std::endl;

  delete r; 
}


//////////////////////////////////////////////////
//////////////////////////////////////////////////
//////////////////////////////////////////////////


void get_frame_rotation(int argc, char *argv[])
{
  if( argc != 14 )
  {
    std::cout << "error usage: test_rotations [" << __func__ << "] "
      << "<mom1> <mom2> <mom_prime1> <mom_prime2>" << std::endl; 
    exit(1);
  }
  
  mom_t m1 = radmat::gen_mom<0,0,0>(); 
  mom_t m2 = radmat::gen_mom<0,0,0>(); 
  mom_t ma = radmat::gen_mom<0,0,0>(); 
  mom_t mb = radmat::gen_mom<0,0,0>(); 

  read_momentum(2,m1,argv);
  read_momentum(5,m2,argv);
  read_momentum(8,ma,argv);
  read_momentum(11,mb,argv);

  std::cout << "read m1 " << string_mom(m1) << " m2 " << string_mom(m2) 
    << " m_prime1 " << string_mom(ma) << " m_prime2 " << string_mom(mb) << std::endl;

   radmat::rCompEulAngles eul_l = radmat::generate_frame_transformation(m1,ma);
   radmat::rCompEulAngles eul_r = radmat::generate_frame_transformation(m2,mb);

  RotationMatrix_t rl = genRotationMatrix(eul_l); 
  RotationMatrix_t rr = genRotationMatrix(eul_r); 
  std::cout << "Rleft:" << rl << std::endl;
  std::cout << "Rright:" << rr << std::endl;

}


//////////////////////////////////////////////////
//////////////////////////////////////////////////
//////////////////////////////////////////////////




//////////////////////////////////////////////////
//////////////////////////////////////////////////
//////////////////////////////////////////////////

void get_left_triad_wigner(int argc, char *argv[])
{
  if( argc != 9 )
  {
    std::cout << "error usage: test_rotations [" << __func__ << "] "
      << "<mom1> <mom2> <J> " << std::endl; 
    exit(1);
  }

  mom_t l = radmat::gen_mom<0,0,0>(); 
  mom_t r = radmat::gen_mom<0,0,0>(); 
  int J; 

  read_momentum(2,l,argv);
  read_momentum(5,r,argv);
  std::istringstream val(argv[8]);
  val >> J; 

  radmat::DMatrixManager Wig; 
  radmat::rCompEulAngles eul = Wig.rotation_matrix(l,r); 
  radmat::RotationMatrix_t * Rtriad = new RotationMatrix_t( genRotationMatrix(eul) ); 
  radmat::WignerMatrix_t *D,*Dtriad; 
  D = Wig.left_wigner_matrix(eul,l,r,J,false,1); 

  clean_up_rot_mat(Rtriad);
  Wig.clean(D); 

  std::pair<radmat::mom_t,radmat::mom_t> can_frame = Wig.get_frame(l,r).moms(); 
  std::cout << "canonical frame: l " << string_mom(l) << " r " << string_mom(r) << std::endl;

  std::cout << "l " << string_mom(l) << " r " << string_mom(r)
    << "\nD_left: " << *D << "\nRtriad:" << *Rtriad << std::endl;  

  delete D; 
  delete Rtriad; 
}



//////////////////////////////////////////////////
//////////////////////////////////////////////////
//////////////////////////////////////////////////

void get_right_triad_wigner(int argc, char *argv[])
{
  if( argc != 9 )
  {
    std::cout << "error usage: test_rotations [" << __func__ << "] "
      << "<mom1> <mom2> <J> " << std::endl; 
    exit(1);
  }

  mom_t l = radmat::gen_mom<0,0,0>(); 
  mom_t r = radmat::gen_mom<0,0,0>(); 
  int J; 

  read_momentum(2,l,argv);
  read_momentum(5,r,argv);
  std::istringstream val(argv[8]);
  val >> J; 


  radmat::DMatrixManager Wig; 
  radmat::rCompEulAngles eul = Wig.rotation_matrix(l,r);
  radmat::RotationMatrix_t * Rtriad(new RotationMatrix_t( genRotationMatrix(eul) )); 

  radmat::WignerMatrix_t *D,*Dtriad; 
  D = Wig.right_wigner_matrix(eul,l,r,J,false,1); 

  clean_up_rot_mat(Rtriad);
  Wig.clean(D); 

  std::pair<radmat::mom_t,radmat::mom_t> can_frame = Wig.get_frame(l,r).moms(); 
  std::cout << "canonical frame: l " << string_mom(l) << " r " << string_mom(r) << std::endl;

  std::cout << "l " << string_mom(l) << " r " << string_mom(r)
    << "\nD_right: " << *D << "\nRtriad:" << *Rtriad << std::endl;  

  delete D; 
  delete Rtriad; 
}

///////////////////////////////////////////////
///////////////////////////////////////////////
///////////////////////////////////////////////
void get_triad_rotation(int argc, char *argv[])
{
  if( argc != 8 )
  {
    std::cout << "error usage: test_rotations [" << __func__ << "] "
      << "<mom1> <mom2> " << std::endl; 
    exit(1);
  }

  mom_t l = radmat::gen_mom<0,0,0>(); 
  mom_t r = radmat::gen_mom<0,0,0>(); 

  read_momentum(2,l,argv);
  read_momentum(5,r,argv);

  radmat::DMatrixManager Wig; 
  radmat::rCompEulAngles eul = Wig.rotation_matrix(l,r);
  radmat::RotationMatrix_t * Rtriad( new RotationMatrix_t( genRotationMatrix(eul) ) ); 

  std::pair<radmat::mom_t,radmat::mom_t> frame; 
  frame = Wig.get_frame(l,r).moms();
  std::cout << "canonical frame " << string_mom(frame.first) << " " 
    << string_mom(frame.second) << std::endl;

  clean_up_rot_mat(Rtriad);

  std::cout << "l " << string_mom(l) << " r " << string_mom(r)
    << "\nRtriad:" << *Rtriad << std::endl;  
  std::cout << "R*cl " << string_mom(rotate_int_mom(Rtriad,frame.first)) 
    << "\nR*cr " << string_mom(rotate_int_mom(Rtriad,frame.second)) << std::endl;
  std::cout << "det(R) = " << determinant(Rtriad) << std::endl;

  delete Rtriad; 
}





//////////////////////////////////////////////////
//////////////////////////////////////////////////
//////////////////////////////////////////////////

void get_wigner_matrix(int argc, char *argv[])
{
  if( argc != 6 )
  {
    std::cout << "error usage: test_rotations [" << __func__ << "] "
      << "<mom> <J> " << std::endl; 
    exit(1);
  }

  mom_t m1 = radmat::gen_mom<0,0,0>(); 
  int J; 

  std::istringstream val(argv[5]);
  val >> J; 

  read_momentum(2,m1,argv);

  std::cout << "read p " << string_mom(m1) 
    << " J " << J << std::endl;
  std::cout << " Wigner_D(p,J) " << std::endl;

  radmat::WignerMatrix_t *W = radmat::WignerDMatrixEnv::call_factory(m1,J);

  std::cout << "W:" << *W << std::endl;

  delete W; 
}

//////////////////////////////////////////////////
//////////////////////////////////////////////////
//////////////////////////////////////////////////

void get_prod_wigner_matrix(int argc, char *argv[])
{
  if( argc != 9 )
  {
    std::cout << "error usage: test_rotations [" << __func__ << "] "
      << "<mom> <mom> <J> " << std::endl; 
    exit(1);
  }

  mom_t m1 = radmat::gen_mom<0,0,0>(); 
  mom_t m2 = radmat::gen_mom<0,0,0>(); 
  int J; 

  std::istringstream val(argv[8]);
  val >> J; 

  read_momentum(2,m1,argv);
  read_momentum(5,m2,argv);

  std::cout << "read p " << string_mom(m1) 
    << " J " << J << std::endl;
  std::cout << " Wigner_D(p,J) * Wigner_D(p2,J) " << std::endl;

  radmat::WignerMatrix_t *W = radmat::WignerDMatrixEnv::call_factory(m1,J);
  radmat::WignerMatrix_t *W2 = radmat::WignerDMatrixEnv::call_factory(m2,J);

  W->lower_index(1); 
  radmat::WignerMatrix_t WW = radmat::contract(*W,*W2,1,0); 

  std::cout << "W:" << WW << std::endl;

  delete W; 
  delete W2; 
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
  insert_op("get_triad_rotation",&get_triad_rotation);
  insert_op("get_left_triad_wigner",&get_left_triad_wigner);
  insert_op("get_right_triad_wigner",&get_right_triad_wigner);
  insert_op("check_rotations_mom",&check_rotations_mom); 
  insert_op("check_lattice_rotations_mom",&check_lattice_rotations_mom); 
  insert_op("get_lattice_rotation",&get_lattice_rotation); 
  insert_op("get_frame_rotation",&get_frame_rotation); 
  insert_op("get_wigner_matrix",&get_wigner_matrix); 
  insert_op("get_prod_wigner_matrix",&get_prod_wigner_matrix); 
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
    std::cerr << "usage: test_rotations : <operation> <op inputs ...> " << std::endl;
    exit(1); 
  }

  std::string op;
  std::istringstream opi(argv[1]); 
  opi >> op; 

  do_work(op,argc,argv); 

  return 0;
}










