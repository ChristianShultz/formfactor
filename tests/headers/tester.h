#ifndef TESTER_H_H_GUARD
#define TESTER_H_H_GUARD

#include <iostream>

struct tester
{
  tester(void) ;

  tester(const char * _test_unit)
  : test_unit(_test_unit) , ct_err(0) , ct_test(0)
  { }

  tester(const tester &o)
  : test_unit(o.test_unit) , ct_err(o.ct_err) , ct_test(o.ct_test)
  {  }

  tester& operator=(const tester &o)
  {
    if(this != &o)
      {
	test_unit = o.test_unit;
	ct_err = o.ct_err;
	ct_test = o.ct_test;
      }
    return *this;
  }

  void operator()(bool success , const char * file, unsigned int line, const char * msg)
  {
    ++ct_test;
    if(!!!success)
      {
	++ ct_err;
	std::cout << file << " : line " << line << " : " << msg << std::endl;
      }
  }
  
  const char * test_unit;
  unsigned short ct_err;
  unsigned short ct_test;
};

#define TESTER_TEST(X,Y,Z) X(Y,__FILE__,__LINE__,Z) 

#endif
