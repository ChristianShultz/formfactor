#ifndef STRINGIFY_H_H_GUARD
#define STRINGIFY_H_H_GUARD

#include <complex>
#include <string>
#include <algorithm>
#include <functional>

/**
  @file stringify.h

  @brief this template is used for testing and prints the type that was tested
  @details ie cout << Stringify<double>() << endl;
  will print double to the terminal, its stupid but its nice for the testing
  */

namespace radmat
{

  struct StringifyChopper
  {
    std::string& l_trim( std::string &s )
    {
      s.erase(
          s.begin(), 
          std::find_if(
            s.begin(),
            s.end(),
            std::not1(std::ptr_fun<int,int>(std::isspace))
            )
          );
      return s; 
    }
  
    std::string& r_trim( std::string &s )
    {
      s.erase( 
          std::find_if(
            s.rbegin(),
            s.rend(),
            std::not1(std::ptr_fun<int,int>(std::isspace))).base(),
          s.end()
          );
      return s; 
    }

    std::string trim( std::string s  )
    {
      return l_trim(r_trim(s)); 
    }

  };


  struct StringifyBase
  {
    StringifyBase() {}
    virtual ~StringifyBase() {}
    virtual std::string name() const = 0; 
  }; 


  template<class T>
    struct StringifyType : public StringifyBase
  { 
    ~StringifyType() {}
  };


  // only specializations may be instatiated
  //   guard against preprocessor whitespace via 
  //   this stupid chopper thing
  template<typename T>
    std::string Stringify(void)
    {
      StringifyType<T> f;
      StringifyChopper chop; 
      return chop.trim( f.name() ); 
    }


#define REGISTER_STRINGIFY_TYPE(X)                \
  template<>                                      \
  struct StringifyType<X> : public StringifyBase  \
  {                                               \
    std::string name() const                      \
    { return std::string( #X  );}                 \
  };                                              \


#define REGISTER_STRINGIFY_TYPE2(X,Y)                 \
  template<>                                          \
  struct StringifyType<X,Y> : public StringifyBase    \
  {                                                   \
    std::string name() const                          \
    { return std::string( #X )                        \
      + "," + std::string( #Y ) ;}                    \
  };                                                  \


#define REGISTER_STRINGIFY_TYPE3(X,Y,Z)                 \
  template<>                                            \
  struct StringifyType<X,Y,Z> : public StringifyBase    \
  {                                                     \
    std::string name() const                            \
    { return std::string( #X ) + "," +                  \
      std::string( #Y ) + "," + std::string( #Z ) ;}    \
  };                                                    \

  // and on and on and on -- this works 

}
#endif
