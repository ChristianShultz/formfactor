/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

* File Name : Wigner_D_matrix_factory.cc

* Purpose :

* Creation Date : 14-04-2014

* Last Modified : Mon 28 Apr 2014 10:49:40 AM EDT

* Created By : shultz

_._._._._._._._._._._._._._._._._._._._._.*/

#include "Wigner_D_matrix_factory.h"
#include "rotation_utils.h"
#include "hadron/irrep_util.h"
#include "hadron/clebsch.h"
#include "semble/semble_meta.h"
#include <sstream>
#include <exception> 
#include "radmat/utils/printer.h"
#include "radmat/utils/euler_angles.h"


namespace radmat
{

  namespace
  {
    std::complex<double> complex_zero(0.,0.); 

    std::complex<double> round_to_zero(const std::complex<double> &cd, const double thresh=1e-6)
    {return ( std::norm(cd) < thresh ) ? complex_zero : cd ; }
    
    WignerMatrix_t gen_wigner_matrix(const mom_t &p, const int J)
    {
      int bound = 2*J + 1; 

      WignerMatrix_t W( (TensorShape<2>())[bound][bound] , std::complex<double>(0.,0.) ); 

      rEulAngles eul = get_euler_angles(p); 

      for(int m1 = -J; m1 <= J; ++m1)
        for(int m2 = -J; m2 <= J; ++m2)
        {
          std::complex<double> cd =  SEMBLE::toScalar( 
              Hadron::Wigner_D( 2*J , 2*m1 , 2*m2 , eul.alpha, eul.beta, eul.gamma) );  
            W[J-m1][J-m2] = round_to_zero(cd,1e-6); 
        }

      return W; 
    }

    WignerKey gen_id(const mom_t &p, const int J)
    {
      return WignerKey(p,J); 
    }

    typedef radmat::WignerDMatrixEnv::TheWignerDMatrixFactory WDFac; 

    bool reg_wigner_matrix(const WignerMatrix_t &W, const WignerKey &id)
    {
      if( WDFac::Instance().find(id) != WDFac::Instance().end() )
        throw std::string("WignerDMatrixEnv double reg error"); 

      WDFac::Instance().insert(std::make_pair(id,W)); 
      return true; 
    }

    std::ostream& operator<<(std::ostream &o, const WignerKey &k)
    {
      return o << "J_" << k.J << "__p" << k.px << k.py << k.pz; 
    }

    WignerMatrix_t* query_factory(const mom_t &p, const int J)
    {
      WignerKey id = gen_id(p,J); 
      std::map<WignerKey,WignerMatrix_t,WignerKeyClassComp>::const_iterator it; 
      it = WDFac::Instance().find(id); 
      if( it == WDFac::Instance().end())
      {
        std::cout << __func__ << ": WignerDMatrixEnv, missing id " << id; 

#pragma omp critical
        {
          std::cout << "avail keys " << std::endl;
          for(it = WDFac::Instance().begin(); it != WDFac::Instance().end(); ++it)
            std::cout << it->first << std::endl;
        }
        throw std::string("WignerDMatrixEnv missed key"); 
      }

      return it->second.clone(); 
    }

    template<int X, int Y, int Z> 
      bool do_reg(const int J)
      {
        mom_t mom = gen_mom<X,Y,Z>(); 
        WignerKey id = gen_id(mom,J); 
        WignerMatrix_t W = gen_wigner_matrix(mom,J);  
        return reg_wigner_matrix(W,id); 
      }


    bool do_mom_reg(const int J)
    {
      bool success = true; 

      success &= do_reg<0,0,0>(J); 

      // D4 
      success &= do_reg<1,0,0>(J); 
      success &= do_reg<0,1,0>(J); 
      success &= do_reg<0,0,1>(J); 
      success &= do_reg<-1,0,0>(J); 
      success &= do_reg<0,-1,0>(J); 
      success &= do_reg<0,0,-1>(J); 

      success &= do_reg<2,0,0>(J); 
      success &= do_reg<0,2,0>(J); 
      success &= do_reg<0,0,2>(J); 
      success &= do_reg<-2,0,0>(J); 
      success &= do_reg<0,-2,0>(J); 
      success &= do_reg<0,0,-2>(J); 

      // D2      
      success &= do_reg<1,1,0>(J); 
      success &= do_reg<0,1,1>(J); 
      success &= do_reg<1,0,1>(J); 
      success &= do_reg<1,-1,0>(J); 
      success &= do_reg<0,1,-1>(J); 
      success &= do_reg<-1,0,1>(J); 
      success &= do_reg<-1,1,0>(J); 
      success &= do_reg<0,-1,1>(J); 
      success &= do_reg<1,0,-1>(J); 
      success &= do_reg<-1,-1,0>(J); 
      success &= do_reg<0,-1,-1>(J); 
      success &= do_reg<-1,0,-1>(J); 


      // D3
      success &= do_reg<1,1,1>(J); 
      success &= do_reg<-1,1,1>(J); 
      success &= do_reg<1,-1,1>(J); 
      success &= do_reg<1,1,-1>(J); 
      success &= do_reg<-1,-1,1>(J); 
      success &= do_reg<1,-1,-1>(J); 
      success &= do_reg<-1,1,-1>(J); 
      success &= do_reg<-1,-1,-1>(J); 

      return success; 
    }

  }  // anonomyous 



  namespace WignerDMatrixEnv
  {

    namespace
    {
      bool local_registration = false; 
    }

    bool registerAll(const int Jmax)
    {
      printer_function<debug_print_reg_all>(__PRETTY_FUNCTION__); 
      bool success = true; 

      if(!!! local_registration)
      {
        for(int J =0 ; J <= Jmax; ++J)
          success &= do_mom_reg(J);  

        local_registration = true; 
      }
      return success; 
    }

    WignerMatrix_t* 
      call_factory(const mom_t &p, const int J)
      {
        return query_factory(p,J); 
      }

  } // WignerDMatrixEnv


} // radmat

