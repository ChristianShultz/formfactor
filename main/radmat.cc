/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

* File Name : radmat.cc

* Purpose :

* Creation Date : 25-02-2013

* Last Modified : Tue 24 Dec 2013 03:23:04 PM EST

* Created By : shultz

_._._._._._._._._._._._._._._._._._._._._.*/


#include <string>
#include <sstream>
#include <iostream>
#include "radmat/register_all/register_all.h"
#include "radmat/driver/radmat_driver.h"

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

  if(argc != 2)
  {
    std::cerr << "usage: radmat : <xmlinigile> " << std::endl;
    exit(1); 
  }

  std::string xmlini;
  std::istringstream val(argv[1]);
  val >> xmlini;

  radmat::RadmatDriver my_driver;
  my_driver.run_program(xmlini);

// grep on this to test completion
  std::cout << "ran successfully" << std::endl; 

  return 0;
}
