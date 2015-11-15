#ifndef AUX_H_H_GUARD
#define AUX_H_H_GUARD

#include "pow2assert.h"
#include "type_computations.h"
#include <string>

namespace radmat
{

  template<class T, class U >
    T downcastAndDeref(const U *ptr)
    {
      const T *p = dynamic_cast<const T*>(ptr);  // cast
      POW2_ASSERT(p);                            // check cast
      T ret = *p;                                // make a T
      return ret; 
    }

  // stl complex is a pain.. , use this to get some asymetric complex operatons done

  template<class T,class U>
    inline typename Promote<T, U>::Type_t convert_stl_type(const T &t)
    {
      return typename Promote<T, U>::Type_t(t);
    }

  template<typename T, typename U, class BinaryOperator>
    inline typename Promote<T, U>::Type_t binary_op(const T &t, const U &u, BinaryOperator b_op)
    {
      return b_op(convert_stl_type<T,U>(t),convert_stl_type<U,T>(u));
    }

  // this is a useful functor for fake data that gets stuck here since it doesn't really have a home
  template<typename T, typename U, typename V, T(*ptr)(const U, const V) >
    struct bind1st_2ParFunction_cc
    {

      bind1st_2ParFunction_cc(void)
        : bound(false)
      {  }

      void bind1st(const U &par1)
      {
        m_par = par1;
        bound = true;
      }

      T operator()(const V &par2) const
      {
        POW2_ASSERT(bound);
        return (*ptr)(m_par,par2);
      }

 /*     bind1st_2ParFunction_cc(const bind1st_2ParFunction_cc<T,U,V,T(*)(const U, const V) > &o)
        :m_par(o.m_par) , bound(o.bound)
      {  }

      bind1st_2ParFunction_cc<T,U,V,T(*)(const U, const V) >& operator=(const bind1st_2ParFunction_cc<T,U,V,T(*)(const U, const V) > &o)
      {
        if(this != &o)
        {
          m_par = o.m_par;
          bound = o.bound;
        }
        return *this;
      }

*/
      private:
      U m_par;
      bool bound;
    };



}

#endif
