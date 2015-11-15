#ifndef LORENTZFF_CANONICAL_JJLIST_FF_H
#define LORENTZFF_CANONICAL_JJLIST_FF_H 


#include "lorentzff_canonical_frame_formfacs_rotation_manager.h"
#include "radmat/utils/obj_expr_t.h"
#include "radmat/utils/stringify.h"
#include <complex>
#include <iostream>

namespace radmat
{
 
  template< typename T>  
  struct cartesianIndex 
  { 
    cartesianIndex() : idx(T(-1)) , h_left(-100) , h_right(-100) {}
    cartesianIndex(const T &i, int hl, int hr) : idx(i) , h_left(hl) , h_right(hr) {}
    T idx;
    int h_left;
    int h_right; 
  };

 template<typename T> 
   std::ostream& operator<<(std::ostream &o, const cartesianIndex<T> &i)
   {
     o << "[ " <<  i.idx  << " ] hl" << i.h_left << " hr " << i.h_right ; 
     return o; 
   }

  typedef ListObjExpr_t< std::complex<double> , cartesianIndex<int> > canIdx_t;  

  template<typename T>
    ObjExpr_t<std::complex<double> , cartesianIndex<T> >
    operator*(const double d, const ObjExpr_t<std::complex<double>, cartesianIndex<T> >  &i)
    {
      return ObjExpr_t<std::complex<double> , cartesianIndex<T> >(  d* i.m_coeff, i.m_obj) ; 
    }

  template<typename T>
    ListObjExpr_t<std::complex<double> , cartesianIndex<T> >
    operator*(const double d, const ListObjExpr_t<std::complex<double>, cartesianIndex<T> >  &i)
    {
      ListObjExpr_t<std::complex<double> , cartesianIndex<T> > ret; 
      if ( fabs( d ) > 1e-6 ) 
      {
        typename ListObjExpr_t<std::complex<double> , cartesianIndex<T> >::const_iterator it; 
        for( it = i.begin(); it != i.end(); ++it)
          ret.push_back(  d * *it ) ; 
      }
      return ret; 
    }

  template<typename T>
    ListObjExpr_t<std::complex<double> , cartesianIndex<T> >
    operator*(const ListObjExpr_t<std::complex<double>, cartesianIndex<T> >  &fop,
        const ListObjExpr_t<std::complex<double>, cartesianIndex<T> >  &)
    {
      std::cout << __PRETTY_FUNCTION__ << ": YOU ARE USING A STUB " << std::endl;
      return fop; 
    }


  template<int J_l, int j_r>
    struct JJFFimpl; 

  REGISTER_STRINGIFY_TYPE2( JJFFimpl<0,0> );
  REGISTER_STRINGIFY_TYPE2( JJFFimpl<1,0> );
  REGISTER_STRINGIFY_TYPE2( JJFFimpl<2,0> );
  REGISTER_STRINGIFY_TYPE2( JJFFimpl<3,0> );
  REGISTER_STRINGIFY_TYPE2( JJFFimpl<0,1> );
  REGISTER_STRINGIFY_TYPE2( JJFFimpl<1,1> );
  REGISTER_STRINGIFY_TYPE2( JJFFimpl<2,1> );
  REGISTER_STRINGIFY_TYPE2( JJFFimpl<3,1> );
  REGISTER_STRINGIFY_TYPE2( JJFFimpl<0,2> );
  REGISTER_STRINGIFY_TYPE2( JJFFimpl<1,2> );
  REGISTER_STRINGIFY_TYPE2( JJFFimpl<2,2> );
  REGISTER_STRINGIFY_TYPE2( JJFFimpl<3,2> );
  REGISTER_STRINGIFY_TYPE2( JJFFimpl<0,3> );
  REGISTER_STRINGIFY_TYPE2( JJFFimpl<1,3> );
  REGISTER_STRINGIFY_TYPE2( JJFFimpl<2,3> );
  REGISTER_STRINGIFY_TYPE2( JJFFimpl<3,3> );



  template<int j_l, int j_r> 
  struct JJFFimpl
    : public FormFacRotationManager<JJFFimpl<j_l,j_r>, canIdx_t, j_l, j_r, DONT_AVERAGE_MASSES>
  {
    typedef canIdx_t Data_t; 
    virtual ~JJFFimpl() {}

    virtual std::string ff_impl() const {return std::string("fooboo");}

    virtual Tensor< Data_t , 1 >
      impl(const Tensor<double,1> &pf,
          const Tensor<double,1> &pi,
          const double momfac,
          const int hl, 
          const int hr) const
      {
        Tensor< Data_t , 1 > ret( (TensorShape<1>())[4] , Data_t()); 
        for(int i = 0; i < 4; ++i)
          ret[i] = Data_t( Data_t::ListObj_t(std::complex<double>(1.,0.) , cartesianIndex<int>(i,hl,hr))); 
        return ret; 
      }

  };


}


#endif /* LORENTZFF_CANONICAL_JJLIST_FF_H */
