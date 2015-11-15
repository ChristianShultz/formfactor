#ifndef DEBUG_PROPS_H
#define DEBUG_PROPS_H 

#include "radmat/utils/pow2assert.h"
#include <iostream>
#include "radmat/utils/handle.h"


// #undef DEBUG_MSG_ON
// #undef DEBUG_HANDLE_ON

#ifdef DEBUG_MSG_ON
#define DEBUG_MSG(X) \
  std::cout << __PRETTY_FUNCTION__ << " " << #X << std::endl
#else
#define DEBUG_MSG(X) ""
#endif

#ifdef DEBUG_HANDLE_ON
#define DEBUG_HANDLE(X) check_handle( X , __PRETTY_FUNCTION__ , __LINE__  ) 
#else
#define DEBUG_HANDLE(X) ""
#endif

namespace
{
    template<typename T> 
      void check_handle(const radmat::rHandle<T> &h, const char *func, const int line )
      {
        if ( ! &*h ) 
          std::cerr << "failing, called by " << func << line << std::endl; 
        POW2_ASSERT( &*h ); 
      }
}

#endif /* DEBUG_PROPS_H */
