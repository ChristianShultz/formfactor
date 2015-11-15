/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

* File Name : simple_test.cc

* Purpose :

* Creation Date : 23-02-2014

* Last Modified : Sun 23 Feb 2014 01:21:10 PM EST

* Created By : shultz

_._._._._._._._._._._._._._._._._._._._._.*/
#include <string>
#include <sstream>
#include <iostream>
#include "radmat/register_all/register_all.h"
#include "radmat/driver/radmat_driver.h"
#include "radmat/llsq/llsq_formfactor_data.h"

using namespace radmat; 

int main(int argc, char *argv[])
{
  try
  {
  radmat::AllFactoryEnv::registerAll(); 
  }
  catch(...)
  {
    std::cout << __func__ << ": reg err " << std::endl;
    exit(1); 
  }


  LLSQRealFormFactorData_t bar; 

  std::cout << bar.esize() << std::endl;
  return 0; 
}

