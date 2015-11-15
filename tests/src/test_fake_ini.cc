// test_fake_ini.cc -
//
// Thursday, June 21 2012
//

#include "radmat/fake_data/fake_data_ini.h"
#include "radmat/utils/splash.h"
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <exception>
#include "io/adat_xmlio.h"

using namespace radmat;

int 
main(int argc, char *argv[])
{
  if(argc != 2)
    {
      SPLASH("usage: test_fake_ini : <xmlinifile> ");
      exit(1);
    }

  std::string xmlinifile;
  std::istringstream val(argv[1]);
  val >> xmlinifile;
  std::cout << "Loading xmlinifile: " << xmlinifile << std::endl;

  FakeDataIni_t fakeinikeys;

  try
    {
      XMLReader xml(xmlinifile);
      read(xml,"/FakeDataIni",fakeinikeys);
    }
  catch(std::exception &e)
    {
      std::cout << "std exception: " << e.what();
    }
  catch(std::string &e)
    {
      std::cout << __func__ << ": ERROR: can't read xmlinifile (" << xmlinifile << "): " << e << std::endl;
      exit(1);
    }
  catch(...)
    {
      SPLASH("An error occured while loading the fake data inifile");
      exit(1);
    }

  std::cout << fakeinikeys << std::endl;

  return 0;
}
