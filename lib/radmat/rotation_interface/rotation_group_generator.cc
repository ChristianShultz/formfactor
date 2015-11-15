/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

* File Name : rotation_group_generator.cc

* Purpose :

* Creation Date : 14-04-2014

* Last Modified : Mon 28 Apr 2014 04:53:01 PM EDT

* Created By : shultz

_._._._._._._._._._._._._._._._._._._._._.*/

#include "rotation_group_generator.h"
#include "radmat/utils/utils.h"



namespace radmat
{

  std::ostream & operator<<(std::ostream &o, const mom_key &k)
  { 
    if( k.isOh )
      return o << "R" << k.x << k.y << k.z ;
    else 
      return o << "F" << k.x << k.y << k.z; 
  }

  std::ostream & operator<<(std::ostream &o, const mom_pair_key &k)
  { return o << "l" << k.l << "_r" << k.r; }


  namespace LatticeRotationEnv
  {

    namespace
    {
      bool local_registration = false; 
    }


    bool registerAll(void)
    {
      printer_function<debug_print_reg_all>(__PRETTY_FUNCTION__); 
      bool success = true; 

      if ( !!! local_registration)
      {
        success &= TheRotationGroupGenerator::Instance().registerAll(); 
        local_registration = true; 
      }

      if( !!! success )
      {
        throw std::string("reg error in LatticeRotationEnv");
      }

      return success; 
    }

    //////////////////////////////////////////////////////////////////////

    mom_pair_key 
      rotation_group_key(const mom_t &l, const mom_t &r)
      {
        return TheRotationGroupGenerator::Instance().canonical_frame(l,r);  
      }
    std::pair<mom_t,mom_t>  
      rotation_group_can_mom(const mom_t &l, const mom_t &r)
      {
        return TheRotationGroupGenerator::Instance().canonical_frame_moms(l,r);  
      }


    std::string 
      rotation_group_label(const mom_t &l, const mom_t &r)
      {
        mom_pair_key k = rotation_group_key(l,r); 
        std::stringstream ss; 
        ss << k;  
        return ss.str(); 
      }

  } // LatticeRotationEnv

} // radmat 

