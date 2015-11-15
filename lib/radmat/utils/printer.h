#ifndef PRINTER_H
#define PRINTER_H 

#include <iostream>
#include <string> 

namespace radmat
{

  // empty turns it off, compiler should also 
  // optimize the function call away, yay
  // comment this away to turn off 
  struct debug_print_0
  {
    static void print(const std::string &msg) 
    { std::cout << msg << std::endl; }
  };

  // comment this away to turn off 
  struct debug_print_1
  {
    static void print(const std::string &msg) 
    { std::cout << msg << std::endl; }
  };

  // comment this away to turn off 
  struct debug_print_2
  {
    static void print(const std::string &msg) 
    { std::cout << msg << std::endl; }
  };

  // comment this away to turn off 
  struct debug_print_3
  {
    static void print(const std::string &msg) 
    { std::cout << msg << std::endl; }
  };

  // comment this away to turn off 
  struct debug_print_4
  {
    static void print(const std::string &msg) 
    { std::cout << msg << std::endl; }
  };

  // comment this away to turn off 
  struct debug_print_5
  {
    static void print(const std::string &msg) 
    { std::cout << msg << std::endl; }
  };

  // comment this away to turn off 
  struct debug_print_reg_all
  {
    static void print(const std::string &msg) 
     {}
    //{ std::cout << msg << std::endl; }
  };

  struct console_print
  {
    static void print(const std::string &msg)
    { std::cout << msg << std::endl; }
  };

  template<class printer>
    void printer_function( const std::string &msg ) 
    {
      printer::print( msg ); 
    } 

}

#endif /* PRINTER_H */
