/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : redstar_canonical_rotations.cc

 * Purpose :

 * Creation Date : 09-12-2013

 * Last Modified : Mon 28 Apr 2014 10:41:52 AM EDT

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/

#include "redstar_canonical_rotations.h"
#include <exception>
#include "radmat/utils/printer.h"
#include "radmat/utils/euler_angles.h"



namespace radmat
{

  namespace CanonicalRotationEnv
  {


    namespace
    {
      template<int X, int Y, int Z> 
        mom_t gen_mom(void)
        {
          mom_t p;
          p.resize(3); 
          p[0] = X;
          p[1] = Y;
          p[2] = Z; 
          return p;
        }


      template<int X, int Y, int Z> 
        RotationMatrix_t * 
        gen_rot_mat(void)
        {
          rEulAngles eul = get_euler_angles( gen_mom<X,Y,Z>() ); 
          return new RotationMatrix_t(genRotationMatrix(eul));  
        }

      // Table VI of 1107.1930
      RotationMatrix_t * 
        gen_rot_ref_mat_D4()
        {
          rEulAngles eul; 
          eul.alpha = 0;
          eul.beta = 0; 
          eul.gamma = 0; 
          return new RotationMatrix_t(genRotationMatrix(eul)); 
        }


      RotationMatrix_t * 
        gen_rot_ref_mat_D2()
        {
          rEulAngles eul; 
          double pi = acos(-1); 
          eul.alpha = pi/2.;
          eul.beta = pi/4.; 
          eul.gamma = -pi/2.; 
          return new RotationMatrix_t(genRotationMatrix(eul)); 
        }

      RotationMatrix_t * 
        gen_rot_ref_mat_D3()
        {
          rEulAngles eul; 
          eul.alpha = acos(-1)/4.;
          eul.beta = acos(1./sqrt(3)); 
          eul.gamma = 0; 
          return new RotationMatrix_t(genRotationMatrix(eul)); 
        }


      bool do_reg(const std::string &id, 
          RotationMatrix_t* (*f)(void) )
      {
        return  TheCanonicalRotationFactory::Instance().registerObject(id,f);
      }


      std::string gen_id(const mom_t &p)
      {
        std::stringstream ss; 
        ss << "p" << p[0] << p[1] << p[2]; 
        return ss.str(); 
      }

      template<int X, int Y, int Z> 
        bool do_reg(void)
        {
          bool success = true; 
          mom_t p;
          p.resize(3);
          p[0] = X; 
          p[1] = Y;
          p[2] = Z; 

          std::string id = gen_id(p);
          success &=  do_reg(id,gen_rot_mat<X,Y,Z>); 

          if( !!! success )
          {
            std::cout << __PRETTY_FUNCTION__ << ": reg failed for " << id << std::endl;
          }
          return success; 
        }

      bool local_registration = false; 
    } // anonomyous 

    bool registerAll(void)
    {
      printer_function<debug_print_reg_all>(__PRETTY_FUNCTION__); 
      bool success = true; 

      if ( !!! local_registration )
      {
        // ref rotations
        success &= do_reg(
            Hadron::generateLittleGroup(gen_mom<0,0,0>()),gen_rot_mat<0,0,0>);
        success &= do_reg(
            Hadron::generateLittleGroup(gen_mom<1,0,0>()),gen_rot_ref_mat_D4); 
        success &= do_reg(
            Hadron::generateLittleGroup(gen_mom<1,1,0>()),gen_rot_ref_mat_D2); 
        success &= do_reg(
            Hadron::generateLittleGroup(gen_mom<1,1,1>()),gen_rot_ref_mat_D3); 

        // Oh
        success &= do_reg<0,0,0>(); 

        // D4 
        success &= do_reg<1,0,0>(); 
        success &= do_reg<0,1,0>(); 
        success &= do_reg<0,0,1>(); 
        success &= do_reg<-1,0,0>(); 
        success &= do_reg<0,-1,0>(); 
        success &= do_reg<0,0,-1>(); 

        success &= do_reg<2,0,0>(); 
        success &= do_reg<0,2,0>(); 
        success &= do_reg<0,0,2>(); 
        success &= do_reg<-2,0,0>(); 
        success &= do_reg<0,-2,0>(); 
        success &= do_reg<0,0,-2>(); 

        // D2      
        success &= do_reg<1,1,0>(); 
        success &= do_reg<0,1,1>(); 
        success &= do_reg<1,0,1>(); 
        success &= do_reg<1,-1,0>(); 
        success &= do_reg<0,1,-1>(); 
        success &= do_reg<-1,0,1>(); 
        success &= do_reg<-1,1,0>(); 
        success &= do_reg<0,-1,1>(); 
        success &= do_reg<1,0,-1>(); 
        success &= do_reg<-1,-1,0>(); 
        success &= do_reg<0,-1,-1>(); 
        success &= do_reg<-1,0,-1>(); 


        // D3
        success &= do_reg<1,1,1>(); 
        success &= do_reg<-1,1,1>(); 
        success &= do_reg<1,-1,1>(); 
        success &= do_reg<1,1,-1>(); 
        success &= do_reg<-1,-1,1>(); 
        success &= do_reg<1,-1,-1>(); 
        success &= do_reg<-1,1,-1>(); 
        success &= do_reg<-1,-1,-1>(); 

        local_registration = true; 
      }

      if(!!!success)
        throw std::string("reg error in CanonicalRotationEnv");

      return success ; 
    }

    RotationMatrix_t* call_factory(const mom_t &p)
    {
      return call_factory(gen_id(p)); 
    }

    RotationMatrix_t* call_factory(const std::string &s)
    {
      return TheCanonicalRotationFactory::Instance().createObject(s); 
    }

  } // CanonicalRotationEnv

}



