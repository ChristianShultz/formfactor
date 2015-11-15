/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

* File Name : redstar_canonical_lattice_rotations.cc

* Purpose :

* Creation Date : 09-12-2013

* Last Modified : Sun 23 Feb 2014 10:35:08 AM EST

* Created By : shultz

_._._._._._._._._._._._._._._._._._._._._.*/



#include "redstar_canonical_rotations.h"
#include "redstar_canonical_lattice_rotations.h"
#include "radmat/utils/handle.h"
#include "hadron/irrep_util.h"
#include <exception>
#include "radmat/utils/printer.h"

namespace radmat
{

  namespace CanonicalLatticeRotationEnv
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
        generate_rotation_matrix(void)
      {
        mom_t p = gen_mom<X,Y,Z>(); 
        std::string LG = Hadron::generateLittleGroup(p);  

        rHandle<RotationMatrix_t> Rref(CanonicalRotationEnv::call_factory(LG)); 
        rHandle<RotationMatrix_t> Rp(CanonicalRotationEnv::call_factory(p)); 

        RotationMatrix_t *R = new Tensor<double,2>( (TensorShape<2>())[4][4], 0. );
        
        for(int mu = 0; mu < 4; ++mu)
          for(int nu = 0; nu < 4; ++nu)
            for(int sum = 0; sum < 4; ++sum)
            (*R)[mu][nu] += (*Rp)[mu][sum] * (*Rref)[nu][sum];

        R->lower_index(1); 

        return R; 
      }

      std::string gen_id(const mom_t &p)
      {
        std::stringstream ss;
        ss << Hadron::generateLittleGroup(p) 
          << ":p" << p[0] << p[1] << p[2];
        return ss.str(); 
      }

    } // anonomyous



    namespace
    {

      bool do_reg(const std::string &s, RotationMatrix_t* (*f)(void))
      {
        return TheCanonicalLatticeRotationFactory::Instance().registerObject(s,f); 
      }
      
      
      template<int X, int Y, int Z>
        bool do_reg(void)
        {
          std::string id =  gen_id(gen_mom<X,Y,Z>()); 
          bool success = do_reg(id,generate_rotation_matrix<X,Y,Z>); 
          if( !!! success )
          {
            std::cout << __PRETTY_FUNCTION__ << ": reg error for " << id << std::endl;
          }
          return success; 
        }

      bool local_registration = false; 
    }

    bool registerAll(void)
    {
      printer_function<debug_print_reg_all>(__PRETTY_FUNCTION__);
      bool success = true; 
      
      if(!!! local_registration )
      {
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
        throw std::string("reg error in CanonicalLatticeRotationEnv"); 

      return success;
    }


    RotationMatrix_t* call_factory(const mom_t &rot)
    {
      return TheCanonicalLatticeRotationFactory::Instance().createObject(gen_id(rot)); 
    }


  } // CanonicalLatticeRotationEnv

} // radmat
