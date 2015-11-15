#ifndef TENSOR_H_H_GUARD
#define TENSOR_H_H_GUARD

#include <exception>
#include <algorithm>
#include <functional>
#include <numeric>
#include <limits>
#include <iostream>
#include <complex>
#include <vector>
#include <cmath>
#include "hadron/irrep_util.h"
#include "pow2assert.h"
#include "type_computations.h"
#include "xml_array.h"
#include "splash.h"
#include "aux.h"

// disclaimer -- there is at least one bug I haven't found..
// stole idea for PrimitiveTensors from http://www.drdobbs.com/184401319?pgno=1



namespace radmat
{

  // declare the index type
  typedef short idx_t;

  // fwd declare some structs
  template<typename T, idx_t N>
    struct PrimitiveTensor;

  template<typename T, idx_t N>
    struct Tensor;

  // fwd declare a unary minus for tensors
  template<typename T, idx_t N>
    Tensor<T, N> operator-(const Tensor<T, N> &);		     

  //! overload multiplication to contract the last index of the first arg with the first index of the second arg
  template<typename T, typename U, idx_t N, idx_t M>
    Tensor < typename Promote<T,U>::Type_t, N + M - 2 > operator*(const Tensor<T, N> &, const Tensor<U, M> &);

  /*
  //! overload multiplication to inner product for vectors
  template<typename T, typename U>
  typename Promote<T,U>::Type_t operator*(const Tensor<T, 1> &, const Tensor<U, 1> &);
   */

  //! an approximately equal object, compare by element 
  template<typename T, idx_t N>
    bool equals(const Tensor<T, N> &, const Tensor<T, N> &, T precision = std::numeric_limits<T>::epsilon());

  //! approximately not equal
  template<typename T, idx_t N>
    bool not_equals(const Tensor<T, N> &a, const Tensor<T, N> &b, T precision = std::numeric_limits<T>::epsilon())
    {
      return !equals(a, b, precision);
    }

  //! actual equivalence 
  template<typename T, idx_t N>
    bool operator==(const Tensor<T, N> &a, const Tensor<T, N> &b)
    {
      return equals(a, b, T(0));
    }

  //! not actual equivalence 
  template<typename T, idx_t N>
    bool operator!=(const Tensor<T, N> &a, const Tensor<T, N> &b)
    {
      return not_equals(a, b, T(0));
    }

  template<typename T>
    T value(const Tensor<T,0> &a)
    {
      return a.value();
    }


  // some useful things like levi civita symbols and g_munu
  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////

  // my gcc doesn't like short * complex or int * complex ...
  //! diag(+,---) -- indicies lowered
  Tensor<double,2> g_dd(void);   

  //! diag(+,---) -- indicies raised
  Tensor<double,2> g_uu(void);


  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////


  /** @ brief a struct generating a tensor shape for use in constructors

    @ details     example usage:

    PrimitiveTensor<double,5> foobar(TensorShape<5>()[4][4][4][4][630] , 7);

    foobar is now a 5th rank tensor, the first 4 indicies run [0,3] and the last runs [0,629], all elements
    are initialized to 7.  This allows us to load the dimensions or shape in the constructor which is nice.
   */


  template<idx_t N>
    struct TensorShape
    {
      TensorShape(void)
      {  }

      TensorShape < N - 1 > operator[](const idx_t dimN) // this basically uses some recursion to
      {
        // map the sequence of elements into a vector
        dim.push_back(dimN);                          // which is then used in the constructor of the
        return TensorShape < N - 1 > (dim);           // PrimitiveTensor class to specify the "shape"
      }

      std::vector<idx_t> vec(void) const {return dim;}
      private:
      TensorShape(const std::vector<idx_t> &dims)
        : dim(dims)
      {  }

      friend struct TensorShape < N + 1 >;

      std::vector<idx_t> dim;
    };

  /*
     @ breif a neat sub array type to deal with indexing multi dimensional tensors

     \tparam T the underlying data type
     \tparam N the sub_t 'rank'

     @details a recursive definition for sub arrays slicing with [] operator on each call
     this will allow access as array[d1][d2]...[dn] w/o needing to specialize for
     each possible value of n, basically a generic tensor

     roughly we are mapping a generic n-tensor to linear storage. each application of []
     can be thought of as providing a new view of the tensor as a slice along that index
     the recursive definition of the template operator [] allows us to keep slicing
     untill we hit a single element, the stop point

     this chains together for access, basically in the code you would call something like
     my_array[a][b][c] where my_array is a narray_t<T,3>, the first application of the
     operator [] will leave you with the slice b,c along a and a sub_t<T,2> object.
     my_array<T,3>[a][b][c] -> my_array<T,2>[b][c] -> my_array<T,1>[c] -> tensor element (T^{abc})
   */

  // a sub array type to deal with indexing of tensors via template recursion -- yay!
  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  template<typename T, idx_t N>
    struct sub_t
    {
      // stl like iterators
      typedef T *iterator;
      typedef const T *const_iterator;

      private:
      T *const pElems;                    //! pointer to elements we have access to at this level
      const idx_t *const pDimensions;     //! pointer to dimensions of this and lower levels
      const idx_t *const pSubArrLen;      //! pointer to sub array length to calculate next pos

      //! private constructor
      sub_t<T, N>(T *_pElems, const idx_t *_pDimensions, const idx_t *_pSubArrLen)
        : pElems(_pElems) , pDimensions(_pDimensions) , pSubArrLen(_pSubArrLen)
      {    }

      public:
      //! same as stl
      iterator begin(void)
      {
        return pElems;
      }
      //! same as stl
      const_iterator begin(void) const
      {
        return pElems;
      }
      //! same as stl
      iterator end(void)
      {
        return pElems + pDimensions[0] * pSubArrLen[0];
      }
      //! same as stl
      const_iterator end(void) const
      {
        return pElems + pDimensions[0] * pSubArrLen[0];
      }
      //! set all elements in this view to ini
      void init(const T &ini = T())
      {
        std::fill(begin(), end(), ini);
      }
      //! recursive copy pattern
      void copy(const sub_t<T, N> &o, const T &ini = T())
      {
        idx_t can_copy = std::min(pDimensions[0], o.pDimensions[0]); //! things we can copy
        idx_t cant_copy = pDimensions[0];                            //! things we cant

        for(idx_t i = 0; i < can_copy; ++i)  // recurse down on the overlapping parts
          (*this)[i].copy(o[i], ini);

        for(idx_t i = can_copy; i < cant_copy; ++i)  // set the non-overlapping bits to ini
          (*this)[i].init(ini);
      }

      //! recursive accessors
      sub_t < T, N - 1 > operator [](idx_t Index)
      {
        POW2_ASSERT(Index < pDimensions[0]);
        return sub_t < T, N - 1 > (&pElems[Index * pSubArrLen[0]], pDimensions + 1, pSubArrLen + 1);
      }

      //! recursive accessors
      const sub_t < T, N - 1 > operator [](idx_t Index) const
      {
        POW2_ASSERT(Index < pDimensions[0]);
        return sub_t < T, N - 1 > (&pElems[Index * pSubArrLen[0]], pDimensions + 1, pSubArrLen + 1);
      }

      //! declare friends
      friend struct PrimitiveTensor < T, N + 1 >;
      friend struct Tensor < T, N + 1 >;
      friend struct sub_t < T, N + 1 >;

      //! recursive stream operator -- printing is a bit funny right now
      friend std::ostream &operator<<(std::ostream &o, const sub_t<T, N> &s)
      {
        o << "[\n" << s[0];

        for(idx_t i = 1; i < s.pDimensions[0]; ++i)
          o << s[i] ;

        o << "]\n";

        return o;
      }
    };

  // sub_t is recursive, we need to specialize the behavior for some base case (N=1)
  /////////////////////////////////////////////////////////////////////////////////////////////////////


  // a sub array type to deal with indexing -- base case  -- this does the real work
  /////////////////////////////////////////////////////////////////////////////////////////////////////////


  /**
    @brief a template partial specialization to specify behavior
    at the recursion stop point, N =1

    \tparam T the underlying data type
    \tparam N the sub_t 'rank'

    @details a recursive definition for sub arrays slicing with [] operator on each call
    this will allow access as array[d1][d2]...[dn] w/o needing to specialize for
    each possible value of n, basically a generic tensor

    roughly we are mapping a generic n-tensor to linear storage. each application of []
    can be thought of as providing a new view of the tensor as a slice along that index
    the recursive definition of the template operator [] allows us to keep slicing
    untill we hit a single element, the stop point

    this chains together for access, basically in the code you would call something like
    my_array[a][b][c] where my_array is a narray_t<T,3>, the first application of the
    operator [] will leave you with the slice b,c along a and a sub_t<T,2> object.
    my_array<T,3>[a][b][c] -> my_array<T,2>[b][c] -> my_array<T,1>[c] -> tensor element (T^{abc})
   */

  template<typename T>
    struct sub_t<T, 1>
    {
      // stl like iterators
      typedef T *iterator;
      typedef const T *const_iterator;

      private:
      T *const pElems;
      const idx_t *const pDimensions;
      const idx_t *const pSubArrLen;

      sub_t<T, 1>(T *_pElems, const idx_t *_pDimensions, const idx_t *_pSubArrLen)
        : pElems(_pElems) , pDimensions(_pDimensions) , pSubArrLen(_pSubArrLen)
      {    }

      public:
      iterator begin(void)
      {
        return pElems;
      }
      const_iterator begin(void) const
      {
        return pElems;
      }
      iterator end(void)
      {
        return pElems + pDimensions[0] * pSubArrLen[0];
      }
      const_iterator end(void) const
      {
        return pElems + pDimensions[0] * pSubArrLen[0];
      }

      void init(const T &ini = T())
      {
        std::fill(begin(), end(), ini);
      }

      void copy(const sub_t<T, 1> &o, const T &ini = T())
      {
        idx_t can_copy = std::min(pDimensions[0], o.pDimensions[0]);
        idx_t cant_copy = pDimensions[0];

        // copy what we can
        std::copy(o.begin(), o.begin() + can_copy, begin());

        // fill anything else
        if(cant_copy > can_copy)
          std::fill(begin() + can_copy, end() , ini);
      }

      // accessors
      T &operator [](idx_t Index)
      {
        POW2_ASSERT((Index < pDimensions[0]) && (Index > -1));
        return pElems[Index];
      }

      const T &operator [](idx_t Index) const
      {
        POW2_ASSERT((Index < pDimensions[0]) && (Index > -1));
        return pElems[Index];
      }

      friend struct PrimitiveTensor<T, 2>;
      friend struct Tensor<T, 2>;
      friend struct sub_t<T, 2>;

      friend std::ostream &operator<<(std::ostream &o, const sub_t<T, 1> &s)
      {
        o << "[" << s[0];

        for(idx_t i = 1; i < s.pDimensions[0]; ++i)
          o << "," << s[i] ;

        o << "]\n";

        return o;
      }
    };



  // a pure virtual ABC for polymorphism
  /////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////

  /**
    @brief A pure virtual base class to derive from
   */
  struct TensorBase
  {
    virtual TensorBase *clone(void) const = 0;
    virtual ~TensorBase(void) {}  
  };



  // a primitive tensor, basically an N-D array type thats resizable and templated
  /////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////

  /**
    @brief the n-dimensional array template

    @details the n-dimensional array type, Tensor will derive from this to implement added tensor functions.
    This is basically on the top of the sub_t stack and is the container for all elements. Note that sub_t
    never actually took control of anything, it was just a convenient way to define the storage.
   */
  template<typename T, idx_t N>
    struct PrimitiveTensor : public TensorBase
  {
    typedef typename sub_t<T, 1>::iterator iterator;             //! stl like iterator to loop elements
    typedef typename sub_t<T, 1>::const_iterator const_iterator; //! stl like const_iterator to loop elements
    typedef T value_type;                                        //! stl like value_type

    //! default constructor
    PrimitiveTensor(void)
      : pElems(NULL) , numElems(0)
    {
      fill_zero();
    }

    //! construct with a set of user provided dimensions
    PrimitiveTensor(const idx_t *_Dimensions, const T &ini = T())
      : pElems(NULL), numElems(0)
    {
      fill_zero();
      resize(_Dimensions, ini);
    }


    //! construct with a tensor shape of user provided dimensions, safer, used in the TensorShape construction method
    PrimitiveTensor(const TensorShape<0> &shape, const T &ini = T())
      : pElems(NULL) , numElems(0)
    {
      fill_zero();

      std::vector<idx_t> _Dimensions = shape.vec();

      POW2_ASSERT(_Dimensions.size() == N);

      for(idx_t i = 0; i < N; ++i)
        Dimensions[i] = _Dimensions[i];

      resize(Dimensions, ini);
    }

    //! construct with a vector of user provided dimensions
    PrimitiveTensor(const std::vector<idx_t> &_Dimensions, const T &ini = T())
      : pElems(NULL) , numElems(0)
    {
      fill_zero();
      POW2_ASSERT(_Dimensions.size() == N);

      for(idx_t i = 0; i < N; ++i)
        Dimensions[i] = _Dimensions[i];

      resize(Dimensions, ini);
    }

    //! copy constructor
    PrimitiveTensor(const PrimitiveTensor<T, N> &o)
      : pElems(NULL) , numElems(0)
    {
      fill_zero();
      PrimitiveTensor<T, N> oo;

      if(!!!o.empty() && oo.resize(o.Dimensions))
        std::copy(o.begin(), o.end(), oo.begin());

      swap(oo);
    }


    //! assignment operator
    PrimitiveTensor<T, N> &operator=(const PrimitiveTensor<T, N> &o)
    {
      if(this != &o)
      {
        PrimitiveTensor<T, N> oo(o);
        swap(oo);
      }

      return *this;
    }


    //! asymmetric assignment
    template<typename T2>
      PrimitiveTensor<T,N>& operator=(const PrimitiveTensor<T2,N> &o)
      {
        if(this == (void*) &o)
          return *this;
        PrimitiveTensor<T,N> oo(o.getDim(),(T)0);
        std::copy(o.begin(),o.end(),oo.begin());
        swap(oo);
        return *this;
      }

    //! virtual destructor cleaning up pointers
    virtual ~PrimitiveTensor<T, N>(void)
    {
      if(!!!empty())
        delete [] pElems;
    }

    //! a clone method
    virtual PrimitiveTensor<T, N> *clone(void) const
    {
      return new PrimitiveTensor<T, N>(*this);
    }

    //! accessor
    sub_t < T, N - 1 > operator [](idx_t Index)
    {
      POW2_ASSERT((Index < Dimensions[0]) && (Index > -1));
      return sub_t < T, N - 1 > (&pElems[Index * SubArrLen[0]], Dimensions + 1, SubArrLen + 1);
    }

    //! accessor
    const sub_t < T, N - 1 > operator [](idx_t Index) const
    {
      POW2_ASSERT((Index < Dimensions[0]) && (Index > -1));
      return sub_t < T, N - 1 > (&pElems[Index * SubArrLen[0]], Dimensions + 1, SubArrLen + 1);
    }


    //! resize a tensor -- can't change rank
    virtual bool resize(const TensorShape<0> &shape, const T &ini = T(), const bool keepElems=false)
    {
      std::vector<idx_t> _Dimensions = shape.vec();
      POW2_ASSERT(_Dimensions.size() == N);

      for(idx_t i = 0; i < N; ++i)
        Dimensions[i] = _Dimensions[i];

      return resize(Dimensions, ini, keepElems);
    }

    //! resize a tensor -- can't change rank
    virtual bool resize(const std::vector<idx_t> _Dimensions, const T &ini = T(), const bool keepElems = false)
    {
      POW2_ASSERT(_Dimensions.size() == N);

      for(idx_t i = 0; i < N; ++i)
        Dimensions[i] = _Dimensions[i];

      return resize(Dimensions, ini, keepElems);
    }

    //! resize a tensor -- can't change rank
    virtual bool resize(const idx_t *_Dimensions, const T &ini = T(), const bool keepElems = false)
    {
      POW2_ASSERT(_Dimensions);
      PrimitiveTensor<T, N> tmp;

      tmp.numElems = 1;

      for(idx_t i = 0; i < N; ++i)
      {
        // if(_Dimensions[i] <= 1)    // 1 is possible in our framework but doesn't make practical sense
        //   POW2_ASSERT(false);

        tmp.numElems *= _Dimensions[i];
        tmp.Dimensions[i] = _Dimensions[i];
        tmp.SubArrLen[i] = 1;

        for(idx_t j = N - 1; j > i; --j)
          tmp.SubArrLen[i] *= _Dimensions[j];
      }

      tmp.pElems = new T[tmp.numElems];

      if(!!!tmp.pElems)
        POW2_ASSERT(false);

      if(!!!empty() && keepElems)
        tmp.copy(*this, ini);
      else
        tmp.initialize(ini);

      swap(tmp);

      return true;
    }

    //! check if we allocated or if this is empty
    virtual bool empty(void) const
    {
      return !pElems;
    }

    //! return the dimension of idx
    virtual idx_t getDim(const idx_t idx) const
    {
      POW2_ASSERT(idx < N);
      return Dimensions[idx];
    }

    //! iterator for looping
    virtual iterator begin(void)
    {
      return pElems;
    }
    //! iterator for looping
    virtual const_iterator begin(void) const
    {
      return pElems;
    }
    //! iterator for looping
    virtual iterator end(void)
    {
      return pElems + numElems;
    }
    //! iterator for looping
    virtual const_iterator end(void) const
    {
      return pElems + numElems;
    }

    //! delete all elements and set all dimensiosn to zero, must be resized to work again
    virtual void clear(void) 
    {
      if(!!!empty())
      {
        delete [] pElems;
        pElems = NULL;
        fill_zero();
      }
    }

    virtual std::vector<idx_t> getDim(void) const
    {
      std::vector<idx_t> foo(N,0);
      for(idx_t i = 0; i < N; ++i)
        foo[i] = Dimensions[i];
      return foo;
    }

    // use everything below this at your own risk..


    //! get elements via pointers, useful for contracting things
    T &getElem(idx_t const *const _pos)
    {
      idx_t pos = 0;

      for(idx_t i = 0; i < N; ++i)
        pos += _pos[i] * SubArrLen[i];

      return pElems[pos];
    }

    //! get elements via pointers, useful for contracting things
    const T &getElem(idx_t const *const _pos) const
    {
      idx_t pos = 0;

      for(idx_t i = 0; i < N; ++i)
        pos += _pos[i] * SubArrLen[i];

      return pElems[pos];
    }

    //! get elements via pointers, useful for contracting things
    T &getElem(idx_t const *const *const _pos)
    {
      idx_t pos = 0;

      for(idx_t i = 0; i < N; ++i)
        pos += *_pos[i] * SubArrLen[i];

      return pElems[pos];
    }

    //! get elements via pointers, useful for contracting things
    const T &getElem(idx_t const *const *const _pos) const
    {
      idx_t pos = 0;

      for(idx_t i = 0; i < N; ++i)
        pos += *_pos[i] * SubArrLen[i];

      return pElems[pos];
    }

    //! set all elements to ini
    virtual void initialize(const T &ini = T())
    {
      std::fill(begin(), end(), ini);
    }

    //! recursively copy the views of the tensor
    virtual void copy(const PrimitiveTensor<T, N> &o, const T &ini = T())
    {
      idx_t can_copy = std::min(Dimensions[0], o.Dimensions[0]);
      idx_t cant_copy = Dimensions[0];

      for(idx_t i = 0; i < can_copy; ++i)  // recurse down on the overlapping parts
        (*this)[i].copy(o[i], ini);

      for(idx_t i = can_copy; i < cant_copy; ++i)  // set the non-overlapping bits to ini
        (*this)[i].init(ini);
    }

    //! swap around the pointers, allow quick copies of tmp data in mem fcns via swap
    virtual void swap(PrimitiveTensor<T, N> &o)
    {
      std::swap(pElems, o.pElems);
      std::swap(numElems, o.numElems);

      std::swap_ranges(Dimensions, Dimensions + N, o.Dimensions);
      std::swap_ranges(SubArrLen, SubArrLen + N, o.SubArrLen);
    }

    //! initialize the Dimensions and SubArrLen arrays to zero
    virtual void fill_zero(void)
    {
      std::fill(Dimensions, Dimensions + N, 0);
      std::fill(SubArrLen, SubArrLen + N, 0);
    }

    T *pElems;              //! pointer to the elements
    idx_t numElems;         //! number of elements in the tensor
    idx_t Dimensions[N];    //! Dimension of each index
    idx_t SubArrLen[N];     //! sub array lengths
  };



  // a primitive tensor, basically an N-D array type -- specialized for N = 1
  /////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////

  /**
    @brief specialization for the 1-dimensional array template

    @details the 1-dimensional array type, Tensor will derive from this to implement added tensor functions.
    This is basically on the top of the sub_t stack and is the container for all elements. Note that sub_t
    doesn't work in this case.
    */


  template<typename T>
    struct PrimitiveTensor<T, 1> : public TensorBase
    {
      typedef typename sub_t<T, 1>::iterator iterator;              //! stl like iterator to loop elements
      typedef typename sub_t<T, 1>::const_iterator const_iterator;  //! stl like const_iterator to loop elements
      typedef T value_type;                                         //! stl like value_type

      //! default constructor
      PrimitiveTensor(void)
        : pElems(NULL) , numElems(0)
      {
        Dimensions[0] = 0;
      }

      //! construct with a set of user provided dimensions
      PrimitiveTensor(const idx_t *_Dimensions, const T &ini = T())
        : pElems(NULL), numElems(0)
      {
        Dimensions[0] = 0;
        resize(_Dimensions, ini);
      }

      //! construct with a vector of user provided dimensions
      PrimitiveTensor(const TensorShape<0> &shape, const T &ini=T())
        : pElems(NULL) , numElems(0)
      {
        std::vector<idx_t> _Dimensions = shape.vec();
        POW2_ASSERT(_Dimensions.size() == 1);
        Dimensions[0] = _Dimensions.at(0);
        resize(Dimensions, ini);
      }

      PrimitiveTensor(const std::vector<idx_t> &_Dimensions, const T &ini =T())
        : pElems(NULL) , numElems(0)
      {
        POW2_ASSERT(_Dimensions.size() == 1);
        Dimensions[0] = _Dimensions.at(0);
        resize(Dimensions, ini);
      }

      //! copy constructor
      PrimitiveTensor(const PrimitiveTensor<T, 1> &o)
        : pElems(NULL) , numElems(0)
      {
        PrimitiveTensor<T, 1> oo;

        if(!!!o.empty() && oo.resize(o.Dimensions))
          std::copy(o.begin(), o.end(), oo.begin());

        swap(oo);
      }


      //! assignment operator
      PrimitiveTensor<T, 1> &operator=(const PrimitiveTensor<T, 1> &o)
      {
        if(this != &o)
        {
          PrimitiveTensor<T, 1> oo(o);
          swap(oo);
        }

        return *this;
      }


      //! asymmetric assignment
      template<typename T2>
        PrimitiveTensor<T,1>& operator=(const PrimitiveTensor<T2,1> &o)
        {
          if(this == (void*) &o)
            return *this;

          PrimitiveTensor<T,1> oo(o.getDim(),(T)0);
          std::copy(o.begin(),o.end(),oo.begin());
          swap(oo);

          return *this;
        }

      //! virtual destructor cleaning up pointers
      virtual ~PrimitiveTensor<T, 1>(void)
      {
        if(!!!empty())
          delete [] pElems;
      }

      //! a clone method
      virtual PrimitiveTensor<T, 1> *clone(void) const
      {
        return new PrimitiveTensor<T, 1>(*this);
      }

      //! accessor
      T &operator [](idx_t Index)
      {
        POW2_ASSERT((Index < Dimensions[0]) && (Index > -1));
        return pElems[Index];
      }

      //! accessor
      const T &operator [](idx_t Index) const
      {
        POW2_ASSERT((Index < Dimensions[0]) && (Index > -1));
        return pElems[Index];
      }

      //! resize a tensor -- can't change rank
      virtual bool resize(const TensorShape<0> &shape, const T &ini = T(), const bool keepElems = false)
      {
        std::vector<idx_t> _Dimensions = shape.vec();
        POW2_ASSERT(_Dimensions.size() == 1);

        Dimensions[0] = _Dimensions.at(0);

        return resize(Dimensions, ini, keepElems);
      }

      //! resize a tensor -- can't change rank
      virtual bool resize(const std::vector<idx_t> &_Dimensions, const T &ini = T(), const bool keepElems = false)
      {
        POW2_ASSERT(_Dimensions.size() == 1);

        Dimensions[0] = _Dimensions.at(0);

        return resize(Dimensions, ini, keepElems);
      }

      //! resize a tensor -- can't change rank
      virtual bool resize(const idx_t *_Dimensions, const T &ini = T(), const bool keepElems = false)
      {
        POW2_ASSERT(_Dimensions);
        PrimitiveTensor<T, 1> tmp;

        tmp.Dimensions[0] = _Dimensions[0];
        tmp.numElems = _Dimensions[0];
        POW2_ASSERT(tmp.numElems > 1);
        tmp.pElems = new T[tmp.numElems];

        if(!!!tmp.pElems)
          POW2_ASSERT(false);

        if(!!!empty() && keepElems)
          tmp.copy(*this, ini);
        else
          tmp.initialize(ini);

        swap(tmp);

        return true;
      }

      //! check if we allocated or if this is empty
      virtual bool empty(void) const
      {
        return !pElems;
      }

      //! return the dimension of idx
      virtual idx_t getDim(const idx_t idx) const
      {
        POW2_ASSERT(idx == 0);
        return Dimensions[0];
      }

      //! iterator for looping
      virtual iterator begin(void)
      {
        return pElems;
      }
      //! iterator for looping
      virtual const_iterator begin(void) const
      {
        return pElems;
      }
      //! iterator for looping
      virtual iterator end(void)
      {
        return pElems + numElems;
      }
      //! iterator for looping
      virtual const_iterator end(void) const
      {
        return pElems + numElems;
      }

      //! delete all elements and set all dimensiosn to zero, must be resized to work again
      virtual void clear(void) 
      {
        if(!!!empty())
        {
          delete [] pElems;
          pElems = NULL;
          Dimensions[0] = 0;
        }
      }

      virtual std::vector<idx_t> getDim(void) const
      {
        std::vector<idx_t> foo(1,0);
        for(idx_t i = 0; i < 1; ++i)
          foo[i] = Dimensions[i];
        return foo;
      }


      // allow for inheritance
      protected:

      //! get elements via pointers, useful for contracting things
      T &getElem(idx_t const   *const _pos)
      {
        return pElems[*_pos];
      }

      //! get elements via pointers, useful for contracting things
      const T &getElem(idx_t const *const _pos) const
      {
        return pElems[*_pos];
      }

      //! get elements via pointers, useful for contracting things
      T &getElem(idx_t const *const *const _pos)
      {
        return pElems[ **_pos];
      }

      //! get elements via pointers, useful for contracting things
      const T &getElem(idx_t const *const *const _pos) const
      {
        return pElems[ **_pos];
      }

      //! set all elements to ini
      virtual void initialize(const T &ini = T())
      {
        std::fill(begin(), end(), ini);
      }

      //! copy the elements
      virtual void copy(const PrimitiveTensor<T, 1> &o, const T &ini = T())
      {
        idx_t can_copy = std::min(Dimensions[0], o.Dimensions[0]);
        idx_t cant_copy = Dimensions[0];

        // copy what we can
        std::copy(o.begin(), o.begin() + can_copy, begin());

        // fill anything else
        if(cant_copy > can_copy)
          std::fill(begin() + can_copy, end() , ini);
      }

      //! swap around the pointers, allow quick copies of tmp data in mem fcns via swap
      virtual void swap(PrimitiveTensor<T, 1> &o)
      {
        std::swap(pElems, o.pElems);
        std::swap(numElems, o.numElems);
        std::swap(Dimensions[0], o.Dimensions[0]);
      }

      T *pElems;             //! pointer to the elements
      idx_t numElems;        //! number of elements
      idx_t Dimensions[1];   //! size of the dimension

    };

  // a primitive tensor, basically an N-D array type -- specialized for N = 0 -- scalar
  /////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////

  /**
    @brief specialization for the 0-dimensional array template

    @details the 0-dimensional array type, Tensor will derive from this to implement added tensor functions.
    This is basically on the top of the sub_t stack and is the container for all elements. Note that sub_t
    doesn't work in this case.
    */

  template<typename T>
    struct PrimitiveTensor<T,0> : public TensorBase
    {
      typedef typename sub_t<T, 1>::iterator iterator;              //! stl like iterator to loop elements
      typedef typename sub_t<T, 1>::const_iterator const_iterator;  //! stl like const_iterator to loop elements
      typedef T value_type;   


      PrimitiveTensor(void)
        : pElem(NULL) , numElem(0) 
      {  }

      PrimitiveTensor(const T &ini)
        : pElem(NULL), numElem(0)
      {
        resize(ini);
      }

      PrimitiveTensor(const TensorShape<0> &shape, const T &ini = T())
        : pElem(NULL) , numElem(0)
      {
        POW2_ASSERT(shape.vec().size() == 0);
        resize(ini);
      }

      PrimitiveTensor(const std::vector<idx_t> &dim, const T &ini = T())
        : pElem(NULL) , numElem(0)
      {
        POW2_ASSERT(dim.size() == 0);
        resize(ini);
      }

      PrimitiveTensor(const PrimitiveTensor<T,0> &o)
        : pElem(NULL) , numElem(0)
      {
        PrimitiveTensor<T,0> oo;

        if(!!!o.empty() && oo.resize(o.pElem) )
          std::copy(o.begin(),o.end(),oo.begin());

        swap(oo);
      }


      PrimitiveTensor<T,0> & operator=(const PrimitiveTensor<T,0> &o)
      {
        if(this != &o)
        {
          PrimitiveTensor<T,0> oo(o);
          swap(oo);
        }
        return *this;
      }

      PrimitiveTensor<T,0> & operator=(const T &ini)
      {
        resize(ini);
      }

      template<typename U>
        PrimitiveTensor<T,0> & operator=(const PrimitiveTensor<U,0> &o)
        {
          if(this == (void*) &o)
            return *this;

          PrimitiveTensor<T,0> oo(T(0));
          std::copy(o.begin(),o.end(),oo.begin());
          swap(oo);

          return *this;
        }

      virtual ~PrimitiveTensor<T,0>(void)
      {
        if(!!!empty())
          delete pElem;
        pElem = NULL;
        numElem = 0;
      }

      virtual PrimitiveTensor<T,0> *clone(void) const
      {
        return new PrimitiveTensor<T,0>(*this);
      }

      T& value(void)
      {
        POW2_ASSERT(!!!empty());
        return *pElem;
      }

      const T& value(void) const
      {
        POW2_ASSERT(!!!empty());
        return *pElem;
      }

      virtual bool resize(const T &ini = T())
      {
        if(empty())
          pElem = new T;
        *pElem = ini;
        numElem = 1;
        return true;
      }

      virtual bool empty(void) const
      {
        return !pElem;
      }

      virtual iterator begin(void)
      {
        return pElem;
      }

      virtual const_iterator begin(void) const
      {
        return pElem;
      }

      virtual iterator end(void)
      {
        return pElem;
      }

      virtual const_iterator end(void) const
      {
        return pElem;
      }

      virtual void clear(void)
      {
        if(!!!empty())
        {
          delete pElem;
          pElem = NULL;
          numElem = 0;
        }
      }

      virtual std::vector<idx_t> getDim(void) const
      {
        return std::vector<idx_t>();
      }

      protected:

      virtual bool resize(const T * ptr)
      {
        if(ptr)
        {
          if(!!!empty())
            *pElem = *ptr;
          else
          {
            pElem = new T();
            POW2_ASSERT(pElem);
            *pElem = *ptr;
            numElem = 1;
          }
        }
        return true;
      }

      virtual void initialize(const T &ini = T())
      {
        resize(ini);
      }

      virtual void copy(const PrimitiveTensor<T,0> &o, const T &ini = T())
      {
        if(o.empty())
          resize(ini);
        else
        {
          PrimitiveTensor<T,0> oo(o);
          swap(oo);
        }
      }

      virtual void swap(PrimitiveTensor<T,0> &o)
      {
        std::swap(pElem,o.pElem);
        std::swap(numElem,o.numElem);
      }

      T *pElem;
      idx_t numElem;
    };



  // the tensor class -- inherit then add in the contractions deal with index positions
  /////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////


  // NB need to use this-> or some type of scoping for methods begin(), end() etc
  // that are defined in the primitive class due to name lookup rules.. poo

  template<typename T, idx_t N>
    struct Tensor : public PrimitiveTensor<T, N>
  {
    // stl junk for inherited begin()/end()
    typedef typename PrimitiveTensor<T, N>::iterator iterator;             //! stl like iterator to loop elements
    typedef typename PrimitiveTensor<T, N>::const_iterator const_iterator; //! stl like const_iterator to loop elements
    typedef typename PrimitiveTensor<T, N>::value_type value_type;         //! stl like value type


    // nested class
    struct equals
    {
      equals(void)
        : eps(std::numeric_limits<T>::epsilon())
      {  }

      equals(T _eps)
        : eps(_eps)
      {  }

      bool operator()(const T &x, const T &y)     // std::abs(eps) for complex types
      {
        // check that the norm of the difference is
        return (std::abs(x - y) <= std::abs(eps)); // <= eps
      }

      T eps;
    };

    // constructors, destructor, assignment

    //! default constructor, defaults to all indicies raised
    Tensor(void)
      : PrimitiveTensor<T, N>()
    {
      is_raised.resize(N, true);
    }

    //! construct with a set of user provided dimensions, defaults to all indicies raised
    Tensor(const idx_t *_Dimensions, const T &ini = T())
      : PrimitiveTensor<T, N>(_Dimensions, ini)
    {
      is_raised.resize(N, true);
    }

    //! construct using a c++ vector of user provided dimensions, default to all indicies raised
    Tensor(const std::vector<idx_t> &_Dimensions, const T &ini=T())
      : PrimitiveTensor<T,N>(_Dimensions, ini)
    {
      is_raised.resize(N,true);
    }

    //! construct with a set of user provided dimensions, defaults to all indicies raised    
    Tensor(const TensorShape<0> &shape, const T &ini = T())
      : PrimitiveTensor<T, N>(shape, ini)
    {
      is_raised.resize(N, true);
    }

    //! copy constuctor
    Tensor(const Tensor<T, N> &o)
      : PrimitiveTensor<T, N>(o) , is_raised(o.is_raised)
    {  }

    //! this needs to be here since compilers are stupid 
    Tensor<T, N> &operator=(const Tensor<T, N> &o)
    {
      if(this != &o)
      {
        try
        {
          PrimitiveTensor<T, N>::operator=(o); // involves a new and if we change
        }                                      // underlying implmentation we don't
        catch(...)                               // have to do anything here
        {
          POW2_ASSERT(false);
        }

        is_raised = o.getRaised();
      }

      return *this;
    }

    //! asymmetric assignment
    template<typename T2>
      Tensor<T,N> & operator=(const Tensor<T2,N> &o)
      {
        if(this == (void*) &o)
          return *this;

        try 
        {
          PrimitiveTensor<T,N>::operator=(o);
        }
        catch(...)
        {
          POW2_ASSERT(false);
        }

        is_raised = o.getRaised();

        return *this;
      }

    //! do nothing destructor for derived bit
    virtual ~Tensor(void) {}

    //! a clone method
    virtual Tensor<T, N> *clone(void) const
    {
      return new Tensor<T, N>(*this);
    }

    //
    // raising/lowering methods

    //! check if an index is raised
    const bool isRaised(const idx_t idx) const
    {
      return is_raised.at(idx);
    }

    //! set all the indicies to the passed vector
    void setRaised(const std::vector<bool> &idxs)
    {
      POW2_ASSERT(idxs.size() == N);
      is_raised = idxs;
    }

    //! raise an index
    void raise_index(const idx_t idx)
    {
      is_raised[idx] = true;
    }

    //! lower an index
    void lower_index(const idx_t idx)
    {
      is_raised[idx] = false;
    }

    //! raise by application of a metric
    virtual void raise(const Tensor<T, 2> &metric, const idx_t idx);

    //! lower by application of a metric
    virtual void lower(const Tensor<T, 2> &metric, const idx_t idx);

    std::vector<bool> getRaised(void) const
    {
      return is_raised;
    }

    //
    // fill all elems to ini
    void fill(const T & ini = T())
    {
      std::fill(this->begin(),this->end(),ini);
    }

    //
    // some basic algebra

    //! multiply by a factor
    Tensor<T,N>& operator*=(const T &factor);
    //! divide by a factor
    Tensor<T,N>& operator/=(const T &factor);
    //! add a constant to all elements
    Tensor<T,N>&  operator+=(const T &factor);
    //! subtract a constant from all elements
    Tensor<T,N>&  operator-=(const T &factor);
    //! subtract another tensor, dimensions and index positions must be the same
    Tensor<T,N>&  operator-=(const Tensor<T, N> &o);
    //! add another tensor, dimensions and index positions must be the same
    Tensor<T,N>&  operator+=(const Tensor<T, N> &o);

    // friends
    //

    //! unary minus
    friend Tensor<T, N> operator-<>(const Tensor<T, N> &plus);

    //! do a dimensional check
    template<typename U, idx_t M>
      friend bool chckDim(const Tensor<U, M> &a, const Tensor<U, M> &b);

    /**
      @brief primitive contraction that ignores indicies -- do not use
      @details useful to piggyback and make other work shorter, specialized
      for M = MM = 1
      */
    template<typename U, typename UU, idx_t M, idx_t MM>
      friend
      Tensor < typename Promote<U, UU>::Type_t, M + MM - 2 >
      contract_p(const Tensor<U, M> &A,
          const Tensor<UU, MM> &B,
          const idx_t A_idx,
          const idx_t B_idx);
    /**
      @brief raise (lower) an index by application of a metric

      @ details ignores indicies on the metric and raises (lowers)
      the index of the return tensor after contracting
      */
    template<typename U, typename UU, idx_t M>
      friend
      Tensor<typename Promote<U, UU>::Type_t, M>
      applyMetric(const Tensor<U, M> &_tensor,
          const Tensor<UU, 2> &metric,
          const idx_t idx);

    //! contract without a metric, raised/lowered checking, can not be used in place of applyMetric
    template<typename U, typename UU, idx_t M, idx_t MM>
      friend
      Tensor < typename Promote<U, UU>::Type_t, M + MM - 2 >
      contract(const Tensor<U, M> &A,
          const Tensor<UU, MM> &B,
          const idx_t A_idx = (M - 1),
          const idx_t B_idx = 0);

    //! contract to inner_product for two rank1 tensors, checks raised/lowered
    template<typename U, typename UU>
      friend
      typename Promote<U, UU>::Type_t
      contract_b(const Tensor<U, 1> &a,
          const Tensor<UU, 1> &b);


    //! contract with a metric, apply the metric if needed
    template<typename U, typename UU, typename UUU, idx_t M, idx_t MM>
      friend
      Tensor < typename Promote<U, typename Promote<UU,UUU>::Type_t >::Type_t, M + MM - 2 >
      contract(const Tensor<U, M> &A,
          const Tensor<UU, MM> &B,
          const Tensor<UUU, 2> &metric,
          const idx_t A_idx = (M - 1),
          const idx_t B_idx = 0);

    //! tensor product
    template<typename U, typename UU, idx_t M, idx_t MM>
      friend 
      Tensor<typename Promote<U,UU>::Type_t, M+MM>
      tensorProduct(const Tensor<U,M> &a, const Tensor<UU,MM> &b);

    //! basic asymmetric algebraic operation
    template<typename U, typename UU, idx_t M>
      friend
      Tensor<typename Promote<U, UU>::Type_t, M>
      operator*(const Tensor<U, M> &lhs, const UU &rhs);


    //! basic asymmetric algebraic operation
    template<typename U, typename UU, idx_t M>
      friend
      Tensor<typename Promote<U, UU>::Type_t, M>
      operator/(const Tensor<U, M> &lhs, const UU rhs);

    /**
      @brief basic asymmetric algebraic operation
      @details indicies of lhs and rhs must be in same positions for this operation to make sense
      */
    template<typename U, typename UU, idx_t M>
      friend
      Tensor<typename Promote<U, UU>::Type_t, M>
      operator+(const Tensor<U, M> &lhs, const Tensor<UU, M> &rhs);

    /**
      @brief basic asymmetric algebraic operation
      @details indicies of lhs and rhs must be in same positions for this operation to make sense
      */
    template<typename U, typename UU, idx_t M>
      friend
      Tensor<typename Promote<U, UU>::Type_t, M>
      operator-(const Tensor<U, M> &lhs, const Tensor<UU, M> &rhs);

    //! stream operator
    template<typename U, idx_t M>
      friend
      std::ostream &operator<<(std::ostream &o, const Tensor<U, M> &t);

    protected:

    //! redefine base class swap
    virtual void swap(Tensor<T, N> &o)
    {
      PrimitiveTensor<T, N>::swap(o);
      std::swap(is_raised, o.is_raised);
    }

    const T& getByIndex(const idx_t i) const 
    {
      return this->pElems[i];
    }

    idx_t size(void) const
    {
      return this->numElems;
    }

    std::vector<bool> is_raised;  //! the index positions, true is raised
  };


  // the tensor class -- specialization for scalars 
  /////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////


  // NB need to use this-> or some type of scoping for methods begin(), end() etc
  // that are defined in the primitive class due to name lookup rules.. poo

  template<typename T>
    struct Tensor<T,0> : public PrimitiveTensor<T,0>
    {

      typedef typename PrimitiveTensor<T,0>::iterator iterator;
      typedef typename PrimitiveTensor<T,0>::const_iterator const_iterator;
      typedef typename PrimitiveTensor<T,0>::value_type value_type;


      // nested class
      struct equals
      {
        equals(void)
          : eps(std::numeric_limits<T>::epsilon())
        {  }

        equals(T _eps)
          : eps(_eps)
        {  }

        bool operator()(const T &x, const T &y)     // std::abs(eps) for complex types
        {
          // check that the norm of the difference is
          return (std::abs(x - y) <= std::abs(eps)); // <= eps
        }

        T eps;
      };


      Tensor(void)
        : PrimitiveTensor<T,0>()
      { }

      Tensor(const T &ini)
        : PrimitiveTensor<T,0>(ini)
      { }

      Tensor(const TensorShape<0> &shape, const T& ini = T())
        : PrimitiveTensor<T,0>(shape,ini)
      { }

      Tensor(const std::vector<idx_t> &dim, const T &ini = T())
        : PrimitiveTensor<T,0>(dim,ini)
      { }

      Tensor(const Tensor<T,0> & o)
        : PrimitiveTensor<T,0>(o)
      {  }

      Tensor<T,0> & operator=(const Tensor<T,0> &o)
      {
        if(this != &o)
        {
          try
          {
            PrimitiveTensor<T,0>::operator=(o);
          }
          catch(...)
          {
            POW2_ASSERT(false);
          }
        }
        return *this;  
      }

      template<typename U>
        Tensor<T,0>& operator=(const Tensor<U,0> &o)
        {
          if(this == (void*) &o)
            return *this;

          try
          {
            PrimitiveTensor<T,0>::operator=(o);
          }
          catch(...)
          {
            POW2_ASSERT(false);
          }

          return *this;

        }

      virtual ~Tensor(void) {}

      void fill(const T& ini=T())
      {
        std::fill(this->begin(),this->end(),ini);
      }


      //
      // some basic algebra

      //! multiply by a factor
      Tensor<T,0>& operator*=(const T &factor);
      //! divide by a factor
      Tensor<T,0>& operator/=(const T &factor);
      //! add a constant to all elements
      Tensor<T,0>&  operator+=(const T &factor);
      //! subtract a constant from all elements
      Tensor<T,0>&  operator-=(const T &factor);
      //! subtract another tensor, dimensions and index positions must be the same
      Tensor<T,0>&  operator-=(const Tensor<T, 0> &o);
      //! add another tensor, dimensions and index positions must be the same
      Tensor<T,0>&  operator+=(const Tensor<T, 0> &o);

      friend Tensor<T,0> operator-<>(const Tensor<T,0> &plus);

      template<typename U, idx_t M>
        friend bool chckDim(const Tensor<U,M> &, const Tensor<U,M> &);

      //! tensor product
      template<typename U, typename UU, idx_t M, idx_t MM>
        friend 
        Tensor<typename Promote<U,UU>::Type_t, M+MM>
        tensorProduct(const Tensor<U,M> &a, const Tensor<UU,MM> &b);

      //! basic asymmetric algebraic operation
      template<typename U, typename UU, idx_t M>
        friend
        Tensor<typename Promote<U, UU>::Type_t, M>
        operator*(const Tensor<U, M> &lhs, const UU &rhs);




      //! basic asymmetric algebraic operation
      template<typename U, typename UU, idx_t M>
        friend
        Tensor<typename Promote<U, UU>::Type_t, M>
        operator/(const Tensor<U, M> &lhs, const UU rhs);

      /**
        @brief basic asymmetric algebraic operation
        @details indicies of lhs and rhs must be in same positions for this operation to make sense
        */
      template<typename U, typename UU, idx_t M>
        friend
        Tensor<typename Promote<U, UU>::Type_t, M>
        operator+(const Tensor<U, M> &lhs, const Tensor<UU, M> &rhs);

      /**
        @brief basic asymmetric algebraic operation
        @details indicies of lhs and rhs must be in same positions for this operation to make sense
        */
      template<typename U, typename UU, idx_t M>
        friend
        Tensor<typename Promote<U, UU>::Type_t, M>
        operator-(const Tensor<U, M> &lhs, const Tensor<UU, M> &rhs);

      //! stream operator
      template<typename U, idx_t M>
        friend
        std::ostream &operator<<(std::ostream &o, const Tensor<U, M> &t);

      T& value(void)
      {
        return *(this->pElem);
      }

      const T& value(void) const
      {
        return *(this->pElem);
      }

      protected:

      idx_t size(void) const
      {
        return this->numElem;
      }


    };



  // implementation of the interesting bits
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////


  // class methods
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  template<typename T, idx_t N>
    void Tensor<T, N>::raise(const Tensor<T, 2> &metric, const idx_t idx)
    {
      Tensor<T,N> tmp = applyMetric(*this, metric, idx);
      swap(tmp);
    }

  template<typename T, idx_t N>
    void Tensor<T, N>::lower(const Tensor<T, 2> &metric, const idx_t idx)
    {
      Tensor<T,N> tmp = applyMetric(*this, metric, idx);
      swap(tmp);
    }

  // the stl is magical
  template<typename T, idx_t N>
    Tensor<T,N>& Tensor<T, N>::operator*=(const T &factor)
    {
      std::transform(this->begin(), this->end(), this->begin(), std::bind1st(std::multiplies<T>(), factor));
      return *this;
    }

  template<typename T, idx_t N>
    Tensor<T,N>&  Tensor<T, N>::operator/=(const T &factor)
    {
      std::transform(this->begin(), this->end(), this->begin(), std::bind1st(std::divides<T>(), factor));
      return *this;

    }

  template<typename T, idx_t N>
    Tensor<T,N>&  Tensor<T, N>::operator+=(const T &factor)
    {
      std::transform(this->begin(), this->end(), this->begin(), std::bind1st(std::plus<T>(), factor));
      return *this;

    }

  template<typename T, idx_t N>
    Tensor<T,N>&  Tensor<T, N>::operator-=(const T &factor)
    {
      std::transform(this->begin(), this->end(), this->begin(), std::bind2nd(std::minus<T>(), factor));
      return *this;

    }

  template<typename T, idx_t N>
    Tensor<T,N>&  Tensor<T, N>::operator-=(const Tensor<T, N> &o)
    {
      POW2_ASSERT((is_raised == o.is_raised)
          && (std::equal(this->Dimensions, this->Dimensions + N, o.Dimensions)));
      std::transform(this->begin(), this->end(), o.begin(), this->begin(), std::minus<T>());
      return *this;

    }

  template<typename T, idx_t N>
    Tensor<T,N>&  Tensor<T, N>::operator+=(const Tensor<T, N> &o)
    {
      POW2_ASSERT((is_raised == o.is_raised)
          && (std::equal(this->Dimensions, this->Dimensions + N, o.Dimensions)));
      std::transform(this->begin(), this->end(), o.begin(), this->begin(), std::plus<T>());
      return *this;

    }

  // friends
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  template<typename T, idx_t N>
    Tensor<T, N> operator-(const Tensor<T, N> &plus)
    {
      Tensor<T, N> minus;
      minus.resize(plus.Dimensions);
      std::transform(plus.begin(), plus.end(), minus.begin(), std::negate<T>());
      return minus;
    }


  //! do a dimensional check
  template<typename T, idx_t N>
    bool chckDim(const Tensor<T, N> &a, const Tensor<T, N> &b)
    {
      if(a.empty())
      {
        if(b.empty()) 
          return true;
        else          
          return false;
      }
      return std::equal(a.Dimensions,a.Dimensions+N,b.Dimensions);
    }


  //! specialize scalars to just check if they are both instantiated/not instantiated
  template<typename T>
    bool chckDim(const Tensor<T,0> &a, const Tensor<T,0> &b)
    {
      if(a.empty() == b.empty())
        return true;
      return false;
    }







  /////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////


  // primitive contractions



  // this is a bit long but its commented at least, it handels the contraction for all other contractions
  template<typename T, typename U, idx_t N, idx_t M>
    Tensor < typename Promote<T,U>::Type_t , N + M - 2 >
    contract_p(const Tensor<T, N> &A, const Tensor<U, M> &B, const idx_t A_idx, const idx_t B_idx)
    {
      POW2_ASSERT(A.Dimensions[A_idx] == B.Dimensions[B_idx]);


      const idx_t c_dim = A.Dimensions[A_idx];    // dimension of the index we are contracting over
      idx_t contract_idx;                         // the contraction index

      // allocate some storage
      idx_t **a, **b, *c;
      a = new idx_t*[N];                          // specifies the index in the A tensor for getElem
      b = new idx_t*[M];                          // same for B
      c = new idx_t[N + M - 2];                   // same for the return tensor

      // set up the pointer arrays for A
      bool minus = false;                         // c will specify the index in the new array

      for(idx_t i = 0; i < N; ++i)           // since it is a contraction we can use the
        if(i != A_idx)                            // bits of c to tell us what elements we needed
          if(!minus)                              // from a and b, using pointers to the c array
            a[i] = &c[i];                         // means that we don't need to worry about
          else                                    // keeping the indicies together.
            a[i] = &c[i - 1];
        else                                      // the set of lines following 'allocate some storage'
        {
          // declare the storage for the c array which is the
          a[i] = &contract_idx;                 // index in the new tensor and then tie the pointers
          minus = true;                         // for a and b to their corresponding index in c
        }

      // set up the pointer arrays for B
      minus = false;

      for(idx_t i = 0; i < M; ++i)
        if(i != B_idx)
          if(!minus)
            b[i] = &c[i + N - 1];
          else
            b[i] = &c[i + N - 2];
        else
        {
          b[i] = &contract_idx;
          minus = true;
        }

      // set up the array to be returned
      Tensor < typename Promote<T,U>::Type_t, N + M - 2 > ret;
      idx_t * shape = new idx_t[N+M-2];
      std::vector<idx_t> dim;

      // get the shape
      for(idx_t i = 0; i < N; ++i)
        if(i != A_idx)
          dim.push_back(A.Dimensions[i]);

      for(idx_t i = 0; i < M; ++i)
        if(i != B_idx)
          dim.push_back(B.Dimensions[i]);


      const idx_t sz =N+M-2;
      // Now fill the return array

      // zero the index explicitly
      for(idx_t i = 0; i < sz; ++i)
      {
        c[i] = 0;
        shape[i] = dim[i];
      }

      // allocate the storage
      POW2_ASSERT(ret.resize(shape,0.));

      // set up an n-dimensional for loop to loop over an fill in all the slots
      while(true)
      {
        // set the element to zero
        ret.getElem(c) = 0;

        // a,b contain pointers to contract_idx so looping over this is like saying
        // T^{abcd} = \sum_e T^{aeb}*T^{ecd}
        for(contract_idx = 0; contract_idx < c_dim; ++contract_idx)
          ret.getElem(c) += 
            binary_op(A.getElem(a),B.getElem(b),std::multiplies<typename Promote<T,U>::Type_t >());

        /*
           ret.getElem(c) += convert_stl_type<T,U>(A.getElem(a)) * convert_stl_type<T,U>

           ((typename Promote<T,U>::Type)A.getElem(a) 
         * (typename Promote<T,U>::Type)B.getElem(b));
         */
        // update idx
        idx_t j;

        for(j = 0; j < sz; ++j)        // the leftmost index cycles fastest
        {
          // perhaps it should be changed to be
          ++c[j];                    // the rightmost based on how the arrays

          if(c[j] < dim[j])          // are stored internally?
            break;

          c[j] = 0;
        }

        // break while loop
        if(j == dim.size())
          break;
      }

      // deallocate
      delete [] shape; 
      delete [] c;     // this was an array of idx_t
      delete [] a;     // this was an array of pointers to idx_t type
      delete [] b;     // same as a

      return ret;
    }




  // inner product to builtin/underlying type
  template<typename T, typename U>
    typename Promote<T,U>::Type_t contract_b(const Tensor<T, 1> &a, const Tensor<U, 1> &b)
    {
      idx_t sz = a.size();
      POW2_ASSERT((a.isRaised(0) != b.isRaised(0))
          && (sz == b.size()));
      typename Promote<T,U>::Type_t sum(0);

      for(idx_t i = 0; i < sz; ++i)
        sum += a.getByIndex(i)*b.getByIndex(i);

      return sum;
    }


  // specializations

  template<typename T, typename U>
    Tensor<typename Promote<T,U>::Type_t,0>
    contract_p(const Tensor<T,1> &A, const Tensor<U,1> &B, const idx_t a_idx, const idx_t b_idx)
    {
      POW2_ASSERT(a_idx == 0);
      POW2_ASSERT(b_idx == 0);
      return Tensor<typename Promote<T,U>::Type_t, 0>(contract_b(A,B));
    }

  template<typename T, typename U>
    Tensor<typename Promote<T,U>::Type_t , 0>
    contract_p(const Tensor<T,0> &a, const Tensor<U,0> &b, const idx_t , const idx_t)
    {
      return Tensor<typename Promote<T,U>::Type_t, 0>(value(a) * value(b));
    }

  template<typename T, typename U, idx_t N>
    Tensor<typename Promote<T,U>::Type_t, N>
    contract_p(const Tensor<T,0> &a, const Tensor<U,N> &b, const idx_t, const idx_t)
    {
      Tensor<typename Promote<T,U>::Type_t, N> foo = b;
      return value(a)*b;
    }


  template<typename T, typename U, idx_t N>
    Tensor<typename Promote<T,U>::Type_t, N>
    contract_p(const Tensor<T,N> &a, const Tensor<U,0> &b, const idx_t, const idx_t)
    {
      Tensor<typename Promote<T,U>::Type_t, N> foo = a;
      return value(b)*a;
    }





  /////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////

  // same structure as above
  template<typename T, typename U, idx_t N>
    Tensor<typename Promote<T,U>::Type_t , N>
    applyMetric(const Tensor<T, N> &tensor, const Tensor<U, 2> &metric, const idx_t idx)
    {

      //    std::cout << __func__ << " T in \n" << tensor << std::endl;
      //    std::cout << __func__ << " metric \n" << metric << std::endl;

      POW2_ASSERT((tensor.Dimensions[idx] == metric.Dimensions[0])
          && (metric.Dimensions[0] == metric.Dimensions[1])
          && (idx < N));

      idx_t **a, **b, *cycle, contract;
      std::vector<idx_t> dim;


      a = new idx_t*[N];
      b = new idx_t*[2];
      cycle = new idx_t[N];

      for(idx_t i = 0; i < N; ++i)
      {
        if(i != idx)
          a[i] = &cycle[i];
        else
          a[i] = &contract;
      }

      b[0] = &contract;
      b[1] = &cycle[idx];

      Tensor<typename Promote<T,U>::Type_t, N> ret;
      POW2_ASSERT(ret.resize(tensor.Dimensions));

      const idx_t cdim = tensor.Dimensions[idx];

      std::fill(cycle, cycle + N, 0);

      while(true)
      {
        ret.getElem(cycle) = 0.;

        for(contract = 0; contract < cdim; ++contract)
          ret.getElem(cycle) = ret.getElem(cycle) +  tensor.getElem(a) * metric.getElem(b);

        idx_t j;

        for(j = 0; j < N; ++j)
        {
          ++cycle[j];

          if(cycle[j] < tensor.Dimensions[j])
            break;

          cycle[j] = 0;
        }

        if(j == N)
          break;
      }

      delete [] cycle;
      delete [] a;
      delete [] b;

      // change the index position 
      ret.setRaised(tensor.is_raised);
      if(tensor.isRaised(idx))
        ret.lower_index(idx);
      else
        ret.raise_index(idx);

      //   std::cout << "T out \n" << ret << std::endl;

      return ret;
    }

  // cause a crash on metrics applied to scalars
  template<typename T, typename U>
    Tensor<typename Promote<T,U>::Type_t ,0>
    applyMetric(const Tensor<T, 0> &tensor, const Tensor<U, 2> &metric, const idx_t idx)
    {
      SPLASH("you were caught trying to apply a metric to a scalar, this isn't supported");
      std::cout << "it is also likely that what you're" 
        <<" doing is mathematically ill-defined" << std::endl;
      POW2_ASSERT(false);

      return Tensor<typename Promote<T,U>::Type_t , 0>(); // gcc warning 
    }


  /////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////


  // real world contractions


  // general contraction -- checks index positions
  template<typename T, typename U, idx_t N, idx_t M>
    Tensor < typename Promote<T,U>::Type_t, N + M - 2 >
    contract(const Tensor<T, N> &A, const Tensor<U, M> &B, const idx_t A_idx, const idx_t B_idx)
    {
      POW2_ASSERT(A.isRaised(A_idx) != B.isRaised(B_idx));

      Tensor < typename Promote<T,U>::Type_t, N + M - 2 > ret = contract_p(A, B, A_idx, B_idx);

      std::vector<bool> raised;

      for(idx_t a = 0; a < N; ++a)
        if(a != A_idx)
          raised.push_back(A.isRaised(a));

      for(idx_t b = 0; b < M; ++b)
        if(b != B_idx)
          raised.push_back(B.isRaised(b));

      ret.setRaised(raised);

      return ret;
    }

  template<typename T, typename U, idx_t N>
    Tensor<typename Promote<T,U>::Type_t, N>  
    contract(const Tensor<T,N> &a, const Tensor<U,0> &b, const idx_t, const idx_t)
    {
      Tensor<typename Promote<T,U>::Type_t , N> foo = a;
      return foo*value(b);
    }

  template<typename T, typename U, idx_t N>
    Tensor<typename Promote<T,U>::Type_t, N>  
    contract(const Tensor<T,0> &a, const Tensor<U,0> &b, const idx_t, const idx_t)
    {
      Tensor<typename Promote<T,U>::Type_t , N> foo = b;
      return foo*value(a);
    }

  template<typename T, typename U>
    Tensor<typename Promote<T,U>::Type_t, 0>
    contract(const Tensor<T,1> &a, const Tensor<U,1> &b, const idx_t aa, const idx_t bb)
    {
      return contract_p(a,b,aa,bb);
    }

  template<typename T, typename U>
    Tensor<typename Promote<T,U>::Type_t , 0>
    contract(const Tensor<T,0> &a, const Tensor<U,0> &b, const idx_t, const idx_t)
    {
      return contract_p(a,b,0,0);
    }

  /////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////


  // real world contractions with metrics

  // apply a metric if we need it, lazy
  template<typename T, typename U, typename V, idx_t N, idx_t M>
    Tensor < typename Promote<T, typename Promote<U, V>::Type_t >::Type_t , N + M - 2 >
    contract(const Tensor<T, N> &A,
        const Tensor<U, M> &B,
        const Tensor<V, 2> &metric,
        const idx_t A_idx = (M -1),
        const idx_t B_idx = 0)
    {
      if(A.isRaised(A_idx) == B.isRaised(B_idx))
      {
        Tensor<typename Promote<T, V>::Type_t , N> AA = applyMetric(A, metric, A_idx);
        return contract(AA, B, A_idx, B_idx);
      }

      return contract(A, B, A_idx, B_idx);
    }


  // the next set all cause crashes since the operations don't make sense to me
  template<typename T, typename U, typename V, idx_t N>
    Tensor < typename Promote<T, typename Promote<U, V>::Type_t >::Type_t , N>
    contract(const Tensor<T, N> &A,
        const Tensor<U, 0> &B,
        const Tensor<V, 2> &metric,
        const idx_t A_idx,
        const idx_t B_idx)
    {
      SPLASH("you were caught trying to apply a metric to a scalar, this isn't supported");
      std::cout << "it is also likely that what you're" 
        <<" doing is mathematically ill-defined" << std::endl;
      POW2_ASSERT(false);

      // gcc warning 
      return Tensor < typename Promote<T, typename Promote<U, V>::Type_t >::Type_t , N>(); 
    }


  template<typename T, typename U, typename V, idx_t N>
    Tensor < typename Promote<T, typename Promote<U, V>::Type_t >::Type_t , N>
    contract(const Tensor<T, 0> &A,
        const Tensor<U, N> &B,
        const Tensor<V, 2> &metric,
        const idx_t A_idx,
        const idx_t B_idx)
    {
      SPLASH("you were caught trying to apply a metric to a scalar, this isn't supported");
      std::cout << "it is also likely that what you're" 
        <<" doing is mathematically ill-defined" << std::endl;
      POW2_ASSERT(false);

      // gcc warning 
      return Tensor < typename Promote<T, typename Promote<U, V>::Type_t >::Type_t , N>(); 
    }


  template<typename T, typename U, typename V>
    Tensor < typename Promote<T, typename Promote<U, V>::Type_t >::Type_t , 0>
    contract(const Tensor<T, 0> &A,
        const Tensor<U, 0> &B,
        const Tensor<V, 2> &metric,
        const idx_t A_idx,
        const idx_t B_idx)
    {
      SPLASH("you were caught trying to apply a metric to a scalar, this isn't supported");
      std::cout << "it is also likely that what you're" 
        <<" doing is mathematically ill-defined" << std::endl;
      POW2_ASSERT(false);

      // gcc warning 
      return Tensor < typename Promote<T, typename Promote<U, V>::Type_t >::Type_t , 0>(); 
    }


  /////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////


  // tensor products



  // a very stupid way to store tensor products..
  template<typename T, typename U, idx_t N, idx_t M>
    Tensor < typename Promote<T,U>::Type_t , N + M >
    tensorProduct(const Tensor<T, N> &a, const Tensor<U, M> &b)
    {
      idx_t **aa, **bb, *cycle, j;

      aa = new idx_t*[N];
      bb = new idx_t*[M];
      cycle = new idx_t[N + M];

      std::vector<idx_t> dim_c(N + M, 0);
      std::vector<bool> raised(N + M, false);

      for(idx_t i = 0; i < N; ++i)
      {
        raised[i] = a.isRaised(i);
        cycle[i] = 0;
        aa[i] = &cycle[i];
        dim_c[i] = a.getDim(i);
      }

      for(idx_t i = 0; i < M; ++i)
      {
        raised[i + N] = b.isRaised(i);
        cycle[i + N] = 0;
        bb[i] = &cycle[i + N];
        dim_c[i + N] = b.getDim(i);
      }

      Tensor < typename Promote<T,U>::Type_t , N + M > ret;
      POW2_ASSERT(ret.resize(dim_c));

      while(true)
      {
        ret.getElem(cycle) = a.getElem(aa) * b.getElem(bb);

        for(j = 0; j < N + M; ++j)
        {
          ++cycle[j];

          if(cycle[j] < dim_c[j])
            break;

          cycle[j] = 0;
        }

        if(j == N + M)
          break;
      }

      delete [] aa;
      delete [] bb;
      delete [] cycle;

      ret.setRaised(raised);

      return ret;
    }


  // tensor producting with a scalar is multiplication by a scalar


  template<typename T, typename U,  idx_t N>
    Tensor < typename Promote<T,U>::Type_t , N  >
    tensorProduct(const Tensor<T, N> &a, const Tensor<U, 0> &b)
    {
      Tensor < typename Promote<T,U>::Type_t , N  > foo = a;
      return foo * value(b);
    }

  template<typename T, typename U,  idx_t N>
    Tensor < typename Promote<T,U>::Type_t , N  >
    tensorProduct(const Tensor<T, 0> &a, const Tensor<U, N> &b)
    {
      Tensor < typename Promote<T,U>::Type_t , N  > foo = b;
      return foo * value(a);
    }





  /////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////



  // multiply by a factor
  template<typename T, typename U, idx_t N>
    Tensor<typename Promote<T,U>::Type_t, N> operator*(const Tensor<T, N> &lhs, const U &rhs)
    {
      Tensor<typename Promote<T,U>::Type_t , N> ret;
      POW2_ASSERT(ret.resize(lhs.Dimensions));

      // if they arent the same type then
      // either T must be convertible to U or U must be convertible to T for this to work
      std::transform(lhs.begin(), lhs.end(), ret.begin(),
          std::bind1st(std::multiplies<typename Promote<T,U>::Type_t>(), rhs));

      ret.setRaised(lhs.is_raised);

      return ret;
    }

  // specialize
  template<typename T, typename U>
    Tensor<typename Promote<T,U>::Type_t, 0> operator*(const Tensor<T, 0> &lhs, const U &rhs)
    {
      Tensor<typename Promote<T,U>::Type_t , 0> ret = lhs;
      ret *= rhs;
      return ret;  
    }


  template<typename T, typename U, idx_t N>
    Tensor<typename Promote<T,U>::Type_t , N> operator*(const U &lhs, const Tensor<T,N> &rhs)
    {
      return rhs * lhs;
    }

  // divide by a factor
  template<typename T, typename U, idx_t N>
    Tensor<typename Promote<T,U>::Type_t, N> operator/(const Tensor<T, N> &lhs, const U &rhs)
    {
      Tensor<typename Promote<T,U>::Type_t , N> ret;
      POW2_ASSERT(ret.resize(lhs.Dimensions));

      // see note above
      std::transform(lhs.begin(), lhs.end(), ret.begin(),
          std::bind1st(std::divides<typename Promote<T,U>::Type() >(), rhs));

      ret.setRaised(lhs.is_raised);

      return ret;
    }


  // specialize
  template<typename T, typename U>
    Tensor<typename Promote<T,U>::Type_t, 0> operator/(const Tensor<T, 0> &lhs, const U &rhs)
    {
      Tensor<typename Promote<T,U>::Type_t , 0> ret = lhs;
      ret /= rhs;
      return ret;  
    }




  // add a tensor
  template<typename T, typename U, idx_t N>
    Tensor<typename Promote<T,U>::Type_t , N> 
    operator+(const Tensor<T, N> &lhs, const Tensor<U, N> &rhs)
    {
      POW2_ASSERT((lhs.is_raised == rhs.is_raised)
          && (std::equal(lhs.Dimensions, lhs.Dimensions + N, rhs.Dimensions)));

      Tensor<typename Promote<T,U>::Type_t , N> ret;
      POW2_ASSERT(ret.resize(lhs.Dimensions));

      // see note above
      std::transform(lhs.begin(), lhs.end(), rhs.begin(), ret.begin(), std::plus<typename Promote<T,U>::Type_t>());
      ret.setRaised(lhs.is_raised);

      return ret;
    }


  // specialize
  template<typename T, typename U>
    Tensor<typename Promote<T,U>::Type_t, 0> 
    operator+(const Tensor<T, 0> &lhs, const Tensor<U,0> &rhs)
    {
      return Tensor<typename Promote<T,U>::Type_t , 0>(value(lhs) + value(rhs));
    }




  // subtract a tensor
  template<typename T, typename U, idx_t N>
    Tensor<typename Promote<T,U>::Type_t, N> operator-(const Tensor<T, N> &lhs, const Tensor<U, N> &rhs)
    {
      POW2_ASSERT((lhs.is_raised == rhs.is_raised)
          && (std::equal(lhs.Dimensions, lhs.Dimensions + N, rhs.Dimensions)));

      Tensor<typename Promote<T,U>::Type_t , N> ret;
      POW2_ASSERT(ret.resize(lhs.Dimensions));

      // see note above
      std::transform(lhs.begin(), lhs.end(), rhs.begin(), ret.begin(), std::minus<typename Promote<T,U>::Type_t>());
      ret.setRaised(lhs.is_raised);

      return ret;
    }



  // specialize
  template<typename T, typename U>
    Tensor<typename Promote<T,U>::Type_t, 0> 
    operator-(const Tensor<T, 0> &lhs, const Tensor<U,0> &rhs)
    {
      return Tensor<typename Promote<T,U>::Type_t , 0>(value(lhs) - value(rhs));
    }






  // stream 
  template<typename T, idx_t N>
    std::ostream &operator<<(std::ostream &o, const Tensor<T, N> &t)
    {
      o << "<";
      std::vector<bool> ipos = t.getRaised();  
      std::vector<bool>::const_iterator it; 
      for(it = ipos.begin(); it != ipos.end(); ++it)
        o << *it; 

      o << ">\n[" << t[0];

      for(idx_t i = 1; i < t.Dimensions[0]; ++i)
        o << t[i];

      o << "]";

      return o;
    }

  // stream specialization for N = 1
  template<typename T>
    std::ostream &operator<<(std::ostream &o, const Tensor<T, 1> &t)
    {
      o << "<";
      std::vector<bool> ipos = t.getRaised();  
      std::vector<bool>::const_iterator it; 
      for(it = ipos.begin(); it != ipos.end(); ++it)
        o << *it; 
      o << ">\n[" << t[0];

      for(idx_t i = 1; i < t.getDim(0); ++i)
        o << "," << t[i];

      o << "]\n";

      return o;
    }

  // specialize N = 0
  template<typename T>
    std::ostream& operator<<(std::ostream &o , const Tensor<T,0> &t)
    {
      o << value(t) << "\n";
      return o;
    }

  // an approximately equals operator
  template<typename T, idx_t N>
    bool equals(const Tensor<T, N> &a, const Tensor<T, N> &b, T precision)
    {
      POW2_ASSERT(chckDim(a, b));

      if(a.getRaised() != b.getRaised())
        return false;
      return std::equal(a.begin(), a.end(), b.begin(),
          typename Tensor<T, N>::equals(precision));
    }

  // specialize
  template<typename T>
    bool equals(const Tensor<T,0> &a, const Tensor<T,0> &b, T precision)
    {
      return std::equal(a.begin(), a.end(), b.begin(),
          typename Tensor<T, 0>::equals(precision));
    }


  // non-friend algebraic stuff
  /////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////



  // overload mul operator to contract
  template<typename T, typename U, idx_t N, idx_t M>
    Tensor < typename Promote<T,U>::Type_t, N + M - 2 >
    operator*(const Tensor<T, N> &a, const Tensor<U, M> &b)
    {
      return contract(a, b, N - 1, 0);
    }

  /*  // overload mul operator to contract
      template<typename T, typename U>
      typename Promote<T,U>::Type_t operator*(const Tensor<T, 1> &a, const Tensor<U, 1> &b)
      {
      return contract(a, b, 0, 0);
      }
      */

  //! overload wedge to tensor product
  template<typename T, typename U, idx_t N, idx_t M>
    Tensor<typename Promote<T,U>::Type_t, N+M>
    operator^(const Tensor<T,N> &a, const Tensor<U,M> &b)
    {
      return tensorProduct(a,b);
    }


  // a conversion routine
  /////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////

  // this should work as an asymmetric type conversion and fail to compile if there is 
  // no implicit conversion available within the code base / compiler
  template<typename T, typename U, idx_t rank>
    Tensor<T,rank> convertTensorUnderlyingType(const Tensor<U,rank> &u)
    {
      Tensor<T,rank> ret;
      ret.setRaised(u.getRaised()); 
      if(u.empty())
        return ret;

      ret.resize(u.getDim());

      typename Tensor<T,1>::iterator result;
      typename Tensor<U,1>::const_iterator first, last;

      first = u.begin();
      last = u.end();
      result = ret.begin();

      while(first != last)
        *result++ = *first++;

      return ret;
    }

  template<typename T, typename U>
    Tensor<T,0> convertTensorUnderlyingType(const Tensor<U,0> &u)
    {
      return Tensor<T,0> (T(value(u)));
    }


} // close radmat namespace
#endif
