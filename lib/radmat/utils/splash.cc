// err_mac.cc -
//
// Wednesday, June 20 2012
//

#include"splash.h"
#include <cstdio>
#include <cstdarg>


void splashMessage::splash::splash(const char * prettyFunction, const char * file, const int line, const char *msg)
{
  std::printf("SPLASH : %s(%d) :" , file, line);
  std::printf("\n    In %s :\n", prettyFunction);
  if(msg)
    std::printf("    %s \n",msg);
}
