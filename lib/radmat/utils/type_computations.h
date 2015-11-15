#ifndef TYPE_COMPUTATIONS_H_H_GUARD
#define TYPE_COMPUTATIONS_H_H_GUARD

#include <complex>

namespace radmat
{

template<typename T, typename U> 
struct Promote
{ 
typedef T Type_t; 
}; 

 

// bool 
template<> 
struct Promote<bool , bool > 
{ 
typedef bool Type_t; 
}; 

template<> 
struct Promote<bool , char > 
{ 
typedef char Type_t; 
}; 

template<> 
struct Promote<bool , short > 
{ 
typedef short Type_t; 
}; 

template<> 
struct Promote<bool , int > 
{ 
typedef int Type_t; 
}; 

template<> 
struct Promote<bool , long > 
{ 
typedef long Type_t; 
}; 

template<> 
struct Promote<bool , float > 
{ 
typedef float Type_t; 
}; 

template<> 
struct Promote<bool , double > 
{ 
typedef double Type_t; 
}; 

template<> 
struct Promote<bool , long double > 
{ 
typedef long double Type_t; 
}; 

template<> 
struct Promote<bool , std::complex<double> > 
{ 
typedef std::complex<double> Type_t; 
}; 

 

// char 
template<> 
struct Promote<char , bool > 
{ 
typedef char Type_t; 
}; 

template<> 
struct Promote<char , char > 
{ 
typedef char Type_t; 
}; 

template<> 
struct Promote<char , short > 
{ 
typedef short Type_t; 
}; 

template<> 
struct Promote<char , int > 
{ 
typedef int Type_t; 
}; 

template<> 
struct Promote<char , long > 
{ 
typedef long Type_t; 
}; 

template<> 
struct Promote<char , float > 
{ 
typedef float Type_t; 
}; 

template<> 
struct Promote<char , double > 
{ 
typedef double Type_t; 
}; 

template<> 
struct Promote<char , long double > 
{ 
typedef long double Type_t; 
}; 

template<> 
struct Promote<char , std::complex<double> > 
{ 
typedef std::complex<double> Type_t; 
}; 

 

// short 
template<> 
struct Promote<short , bool > 
{ 
typedef short Type_t; 
}; 

template<> 
struct Promote<short , char > 
{ 
typedef short Type_t; 
}; 

template<> 
struct Promote<short , short > 
{ 
typedef short Type_t; 
}; 

template<> 
struct Promote<short , int > 
{ 
typedef int Type_t; 
}; 

template<> 
struct Promote<short , long > 
{ 
typedef long Type_t; 
}; 

template<> 
struct Promote<short , float > 
{ 
typedef float Type_t; 
}; 

template<> 
struct Promote<short , double > 
{ 
typedef double Type_t; 
}; 

template<> 
struct Promote<short , long double > 
{ 
typedef long double Type_t; 
}; 

template<> 
struct Promote<short , std::complex<double> > 
{ 
typedef std::complex<double> Type_t; 
}; 

 

// int 
template<> 
struct Promote<int , bool > 
{ 
typedef int Type_t; 
}; 

template<> 
struct Promote<int , char > 
{ 
typedef int Type_t; 
}; 

template<> 
struct Promote<int , short > 
{ 
typedef int Type_t; 
}; 

template<> 
struct Promote<int , int > 
{ 
typedef int Type_t; 
}; 

template<> 
struct Promote<int , long > 
{ 
typedef long Type_t; 
}; 

template<> 
struct Promote<int , float > 
{ 
typedef float Type_t; 
}; 

template<> 
struct Promote<int , double > 
{ 
typedef double Type_t; 
}; 

template<> 
struct Promote<int , long double > 
{ 
typedef long double Type_t; 
}; 

template<> 
struct Promote<int , std::complex<double> > 
{ 
typedef std::complex<double> Type_t; 
}; 

 

// long 
template<> 
struct Promote<long , bool > 
{ 
typedef long Type_t; 
}; 

template<> 
struct Promote<long , char > 
{ 
typedef long Type_t; 
}; 

template<> 
struct Promote<long , short > 
{ 
typedef long Type_t; 
}; 

template<> 
struct Promote<long , int > 
{ 
typedef long Type_t; 
}; 

template<> 
struct Promote<long , long > 
{ 
typedef long Type_t; 
}; 

template<> 
struct Promote<long , float > 
{ 
typedef float Type_t; 
}; 

template<> 
struct Promote<long , double > 
{ 
typedef double Type_t; 
}; 

template<> 
struct Promote<long , long double > 
{ 
typedef long double Type_t; 
}; 

template<> 
struct Promote<long , std::complex<double> > 
{ 
typedef std::complex<double> Type_t; 
}; 

 

// float 
template<> 
struct Promote<float , bool > 
{ 
typedef float Type_t; 
}; 

template<> 
struct Promote<float , char > 
{ 
typedef float Type_t; 
}; 

template<> 
struct Promote<float , short > 
{ 
typedef float Type_t; 
}; 

template<> 
struct Promote<float , int > 
{ 
typedef float Type_t; 
}; 

template<> 
struct Promote<float , long > 
{ 
typedef float Type_t; 
}; 

template<> 
struct Promote<float , float > 
{ 
typedef float Type_t; 
}; 

template<> 
struct Promote<float , double > 
{ 
typedef double Type_t; 
}; 

template<> 
struct Promote<float , long double > 
{ 
typedef long double Type_t; 
}; 

template<> 
struct Promote<float , std::complex<double> > 
{ 
typedef std::complex<double> Type_t; 
}; 

 

// double 
template<> 
struct Promote<double , bool > 
{ 
typedef double Type_t; 
}; 

template<> 
struct Promote<double , char > 
{ 
typedef double Type_t; 
}; 

template<> 
struct Promote<double , short > 
{ 
typedef double Type_t; 
}; 

template<> 
struct Promote<double , int > 
{ 
typedef double Type_t; 
}; 

template<> 
struct Promote<double , long > 
{ 
typedef double Type_t; 
}; 

template<> 
struct Promote<double , float > 
{ 
typedef double Type_t; 
}; 

template<> 
struct Promote<double , double > 
{ 
typedef double Type_t; 
}; 

template<> 
struct Promote<double , long double > 
{ 
typedef long double Type_t; 
}; 

template<> 
struct Promote<double , std::complex<double> > 
{ 
typedef std::complex<double> Type_t; 
}; 

 

// long double 
template<> 
struct Promote<long double , bool > 
{ 
typedef long double Type_t; 
}; 

template<> 
struct Promote<long double , char > 
{ 
typedef long double Type_t; 
}; 

template<> 
struct Promote<long double , short > 
{ 
typedef long double Type_t; 
}; 

template<> 
struct Promote<long double , int > 
{ 
typedef long double Type_t; 
}; 

template<> 
struct Promote<long double , long > 
{ 
typedef long double Type_t; 
}; 

template<> 
struct Promote<long double , float > 
{ 
typedef long double Type_t; 
}; 

template<> 
struct Promote<long double , double > 
{ 
typedef long double Type_t; 
}; 

template<> 
struct Promote<long double , long double > 
{ 
typedef long double Type_t; 
}; 

template<> 
struct Promote<long double , std::complex<double> > 
{ 
typedef std::complex<double> Type_t; 
}; 

 

// std::complex<double> 
template<> 
struct Promote<std::complex<double> , bool > 
{ 
typedef std::complex<double> Type_t; 
}; 

template<> 
struct Promote<std::complex<double> , char > 
{ 
typedef std::complex<double> Type_t; 
}; 

template<> 
struct Promote<std::complex<double> , short > 
{ 
typedef std::complex<double> Type_t; 
}; 

template<> 
struct Promote<std::complex<double> , int > 
{ 
typedef std::complex<double> Type_t; 
}; 

template<> 
struct Promote<std::complex<double> , long > 
{ 
typedef std::complex<double> Type_t; 
}; 

template<> 
struct Promote<std::complex<double> , float > 
{ 
typedef std::complex<double> Type_t; 
}; 

template<> 
struct Promote<std::complex<double> , double > 
{ 
typedef std::complex<double> Type_t; 
}; 

template<> 
struct Promote<std::complex<double> , long double > 
{ 
typedef std::complex<double> Type_t; 
}; 

template<> 
struct Promote<std::complex<double> , std::complex<double> > 
{ 
typedef std::complex<double> Type_t; 
}; 

}

#endif
