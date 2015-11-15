/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : grab_ff.cc

 * Purpose :

 * Creation Date : 22-04-2013

 * Last Modified : Fri 25 Jul 2014 02:04:47 PM EDT

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/


#include <string>
#include <sstream>

#include "radmat/register_all/register_all.h"
#include "radmat/ff_interface/ff_interface.h"
#include "radmat/rotation_interface/rotation_interface.h"
#include "io/adat_xmlio.h"


struct chunk
{
  std::string origin, data;
  double E; 
  int row; 
  ADATXML::Array<int> mom; 
};

struct input
{
  chunk left, gamma, right; 
  std::string tran; 
  double mom_fac; 
};


  template<typename T>
void doXMLReadT(ADATXML::XMLReader &ptop, const std::string &path, T &place, const char * f)
{
  if(ptop.count(path) > 0)
    read(ptop,path,place);
  else
  {
    std::cerr << __PRETTY_FUNCTION__ << ": Error, called by " << f << " trying to read path, " << path
      << ", path was empty, exiting" << std::endl;
    exit(1);
  }
}


void read( ADATXML::XMLReader &xml, const std::string &pth, chunk &c)
{
  ADATXML::XMLReader ptop(xml,pth); 
  doXMLReadT(ptop,"origin",c.origin,__PRETTY_FUNCTION__);
  doXMLReadT(ptop,"data",c.data,__PRETTY_FUNCTION__);
  doXMLReadT(ptop,"mass",c.E,__PRETTY_FUNCTION__);
  doXMLReadT(ptop,"row",c.row,__PRETTY_FUNCTION__);
  doXMLReadT(ptop,"mom",c.mom,__PRETTY_FUNCTION__);
}; 


int mom_sq(const ADATXML::Array<int> &p)
{
  return p[0]*p[0] + p[1]*p[1] + p[2]*p[2]; 
}

void read( ADATXML::XMLReader &xml, const std::string &pth, input &i)
{
  ADATXML::XMLReader ptop(xml,pth); 
  doXMLReadT(ptop,"left",i.left,__PRETTY_FUNCTION__);
  doXMLReadT(ptop,"gamma",i.gamma,__PRETTY_FUNCTION__);
  doXMLReadT(ptop,"right",i.right,__PRETTY_FUNCTION__);
  doXMLReadT(ptop,"tran",i.tran,__PRETTY_FUNCTION__);
  doXMLReadT(ptop,"mom_fac",i.mom_fac,__PRETTY_FUNCTION__);


  i.left.E = sqrt( i.left.E * i.left.E + i.mom_fac*i.mom_fac*double(mom_sq(i.left.mom)));
  i.gamma.E = sqrt( i.gamma.E * i.gamma.E + i.mom_fac*i.mom_fac*double(mom_sq(i.gamma.mom)));
  i.right.E = sqrt( i.right.E * i.right.E + i.mom_fac*i.mom_fac*double(mom_sq(i.right.mom)));
}


  ENSEM::EnsemReal 
ensem_energy( const double d, const int sz)
{
  ENSEM::EnsemReal ret; 
  ret.resize(sz); 
  ret = ENSEM::Real(d); 
  return ret; 
}



  radmat::ThreePointDataTag
generate_tag( const chunk &l , 
    const chunk &g, 
    const chunk &r, 
    const std::string &tran, 
    const double mom_fac)
{
  radmat::ThreePointDataTag foo; 
  foo.origin_rep = radmat::DataRep3pt( l.origin, g.origin, r.origin); 
  foo.data_rep = radmat::DataRep3pt( l.data, g.data, r.data); 
  foo.mat_elem_id = l.origin + r.origin + "_" + tran; 

  if( (l.data != l.origin ) || ( r.data != r.origin) )
    foo.mat_elem_id += "__" + l.data + "," + r.data; 

  // make it consistently 1 based 
  foo.left_row = l.row; 
  foo.gamma_row = g.row; 
  foo.right_row = r.row; 

  foo.left_mom = l.mom;
  foo.q = g.mom;
  foo.right_mom = r.mom; 

  foo.mom_fac = mom_fac; 

  foo.left_E = ensem_energy( l.E , 2 ); 
  foo.right_E = ensem_energy( r.E , 2 ); 

  return foo;
}


int main(int argc, char *argv[])
{
  radmat::AllFactoryEnv::registerAll(); 

  chunk left, right , gamma ;

  if( argc != 2 )
  {
    std::cout << "error: usage grab_ff <xml_ini>" << std::endl; 
    exit(1); 
  }

  input inp; 
  std::istringstream val(argv[1]); 
  std::string ini; 
  val >> ini; 

  try 
  {
    ADATXML::XMLReader xml(ini); 
    doXMLReadT(xml,"/Props",inp,__PRETTY_FUNCTION__); 
  }
  catch(...)
  {
    std::cout << "something went wrong" << std::endl;
    exit(1); 
  }


  radmat::ThreePointDataTag foo = generate_tag( inp.left , inp.gamma , inp.right , inp.tran , inp.mom_fac ); 

  radmat::FFKinematicFactors_t KK; 
  radmat::FFKinematicFactors_t::KinematicFactorMatrix rm = KK.genFactorsMat( &foo ); 
  radmat::FFKinematicFactors_t::KinematicFactorRow r = KK.genFactors( &foo ); 

  std::cout << "mat:\n" << rm.mean() << std::endl;
 
  std::cout << "pulled:" << std::endl;
  std::cout << r.mean() << std::endl; 

  radmat::DMatrixManager Wig; 
  std::pair<radmat::mom_t,radmat::mom_t> can = Wig.get_frame(foo.left_mom,foo.right_mom).moms(); 

  std::cout << " canonical frame was " 
    << " lefty p" << can.first[0] << can.first[1] << can.first[2] 
    << " righty p" << can.second[0] << can.second[1] << can.second[2] << std::endl;

  return 0;
}
