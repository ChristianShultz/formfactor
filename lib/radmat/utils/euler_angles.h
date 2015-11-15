#ifndef EULER_ANGLES_H
#define EULER_ANGLES_H 

#include "hadron/irrep_util.h"
#include "io/adat_xmlio.h"
#include "tensor.h"
#include <vector>

namespace radmat
{

  // struct to hold euler angles 
  struct rEulAngles
  {
    rEulAngles() 
      : alpha(0.) , beta(0.) , gamma(0.)
    {}
    rEulAngles( const double a, const double b , const double g)
      : alpha(a) , beta(b) , gamma(g)
    {}
    rEulAngles( const Hadron::CubicCanonicalRotation_t &eul)
      : alpha(eul.alpha) , beta(eul.beta) , gamma(eul.gamma)
    {}
    rEulAngles( const rEulAngles &o)
      : alpha(o.alpha) , beta(o.beta) , gamma(o.gamma) 
    {}

    rEulAngles& operator=(const rEulAngles &o)
    {
      if( this != &o)
      {
        alpha = o.alpha; 
        beta = o.beta; 
        gamma = o.gamma; 
      }

      return *this; 
    }

    double alpha, beta, gamma;
  };


  // composite transformation 
  struct rCompEulAngles
  {
    typedef std::vector<rEulAngles> tform_t; 
    typedef tform_t::const_iterator const_iterator; 

    rCompEulAngles() {}

    rCompEulAngles( const rEulAngles &r)
    { tform.push_back(r); }

    rCompEulAngles( const rCompEulAngles &o)
      : tform(o.tform)
    { }

    rCompEulAngles& operator=(const rCompEulAngles &o)
    {
      if( this != &o)
      {
        tform = o.tform; 
      }
      return *this; 
    }

    const_iterator begin() const { return tform.begin(); }
    const_iterator end() const { return tform.end(); }

    void push_back(const rEulAngles &r) { tform.push_back(r); }

    tform_t tform; 
  };

  rEulAngles get_euler_angles(const ADATXML::Array<int> &); 

  // dagger it 
  rEulAngles inverse(const rEulAngles &); 
  rCompEulAngles inverse(const rCompEulAngles &); 


  Tensor<double,2> genRotationMatrix(const rEulAngles &); 
  Tensor<double,2> genRotationMatrix(const rCompEulAngles &); 

}

#endif /* EULER_ANGLES_H */
