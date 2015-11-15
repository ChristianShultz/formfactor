/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

* File Name : euler_angles.cc

* Purpose :

* Creation Date : 28-04-2014

* Last Modified : Mon 28 Apr 2014 11:33:54 AM EDT

* Created By : shultz

_._._._._._._._._._._._._._._._._._._._._.*/

#include "euler_angles.h"


namespace radmat
{

  rEulAngles inverse( const rEulAngles &eul)
  {
    return rEulAngles( -eul.gamma , -eul.beta, -eul.alpha); 
  }

  rCompEulAngles inverse(const rCompEulAngles &eul)
  {
    rCompEulAngles ret; 
    rCompEulAngles::const_iterator e; 

    for( e = eul.begin(); e != eul.end(); ++e)
      ret.push_back( inverse( *e ) ); 

    return ret; 
  }

  namespace
  {
    Tensor<double,2> I4(void)
    {
      Tensor<double,2> I((TensorShape<2>())[4][4],0.);
      for(int i =0; i < 4; ++i)
        I[i][i] = 1.;
      I.lower_index(1);
      return I;
    }

    bool isRest(const XMLArray::Array<int> &mom)
    {
      return ((mom[0] == 0) && (mom[1] == 0) && (mom[2] == 0));
    }

  } // anonomyous 


  // pull from adat 
  rEulAngles get_euler_angles(const ADATXML::Array<int> &p)
  {
    Hadron::CubicCanonicalRotation_t eul;
    if( (p[0] == 0) && (p[1] == 0) && (p[2] == 0) )
    {
      eul.alpha = 0.;
      eul.beta = 0.;
      eul.gamma = 0.; 
    } 
    else
      eul = Hadron::cubicCanonicalRotation(p); 

    return rEulAngles(eul);
  }



  // Z-Y-Z
  Tensor<double,2> genRotationMatrix(const rEulAngles &eulerangles)
  {

    Tensor<double,2> A((TensorShape<2>())[4][4],0.),B((TensorShape<2>())[4][4],0.),C((TensorShape<2>())[4][4],0.);


    double a,b,g;
    a = eulerangles.alpha;
    b = eulerangles.beta;
    g = eulerangles.gamma;


    A[0][0] = 1.;
    B[0][0] = 1.;
    C[0][0] = 1.;


    A[1][1] = cos(a)   ;    A[1][2] = -sin(a)   ;    A[1][3] =  0.      ; 
    A[2][1] = sin(a)   ;    A[2][2] = cos(a)    ;    A[2][3] =  0.      ; 
    A[3][1] = 0.       ;    A[3][2] = 0.        ;    A[3][3] =  1.      ; 


#if 0
    // B_x
    B[1][1] = 1.       ;    B[1][2] = 0.        ;    B[1][3] =  0.      ; 
    B[2][1] = 0.       ;    B[2][2] = cos(b)    ;    B[2][3] = sin(b)   ; 
    B[3][1] = 0.       ;    B[3][2] = -sin(b)   ;    B[3][3] = cos(b)   ; 
#endif 

    // B_y 
    B[1][1] = cos(b)   ;    B[1][2] = 0.        ;    B[1][3] = sin(b)   ; 
    B[2][1] = 0.       ;    B[2][2] = 1.        ;    B[2][3] =  0.      ; 
    B[3][1] = -sin(b)  ;    B[3][2] = 0.        ;    B[3][3] = cos(b)   ; 



    C[1][1] = cos(g)   ;    C[1][2] = -sin(g)   ;    C[1][3] =  0.      ; 
    C[2][1] = sin(g)   ;    C[2][2] = cos(g)    ;    C[2][3] =  0.      ; 
    C[3][1] = 0.       ;    C[3][2] = 0.        ;    C[3][3] =  1.      ; 



    A.lower_index(1);
    B.lower_index(1);
    C.lower_index(1);

    // return C * B * A;
    return A*B*C;
  }

  Tensor<double,2> genRotationMatrix(const rCompEulAngles &eul)
  {
    Tensor<double,2> ret,tmp;
    ret = I4(); 

    rCompEulAngles::const_iterator e; 

    for( e = eul.begin(); e != eul.end(); ++e)
    {
      tmp = ret * genRotationMatrix(*e); 
      ret = tmp; 
    }

    return ret; 
  }


} // radmat 

