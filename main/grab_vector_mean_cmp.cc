/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

* File Name : grab_mean.cc

* Purpose :

* Creation Date : 02-05-2013

* Last Modified : Thu May  2 13:58:50 2013

* Created By : shultz

_._._._._._._._._._._._._._._._._._._._._.*/


#include <string>
#include <iostream>
#include "ensem/ensem.h"
#include "semble/semble_meta.h"
#include "io/adat_io.h"
#include "io/adat_xmlio.h"





int main(int argc , char *argv[])
{

  if(argc != 4)
  {
    std::cerr << "usage: " << argv[0] << " : <file.jack>  <tlow> <thigh> " << std::endl;
    exit(1);
  }

  std::string file;
  int tlow; 
  int thigh;

  {std::stringstream val(argv[1]); val >> file;}
  {std::stringstream val(argv[2]); val >> tlow;}
  {std::stringstream val(argv[3]); val >> thigh;}


  ENSEM::EnsemVectorComplex tmp;
  ENSEM::EnsemVectorReal rl, cmp;
  double mr(0), mc(0); 

  ENSEM::read(file,tmp);
  rl = ENSEM::real(tmp);
  cmp = ENSEM::imag(tmp);


  for(int i = tlow; i < thigh; ++i)
  {
    mr += SEMBLE::toScalar(ENSEM::mean(ENSEM::peekObs(rl,i))); 
    mc += SEMBLE::toScalar(ENSEM::mean(ENSEM::peekObs(cmp,i))); 
  }


  mr /= double(thigh - tlow);
  mc /= double(thigh - tlow);


  std::cout << mr << " " << mc << std::endl;


  return 0;
}



