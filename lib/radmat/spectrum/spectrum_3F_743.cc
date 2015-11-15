/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : spectrum_3F_743.cc

 * Purpose :

 * Creation Date : 06-05-2014

 * Last Modified : Tue 06 May 2014 06:35:53 PM EDT

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/

#include "spectrum_3F_743.h"
#include "spectrum_factory.h"

typedef ::radmat::TheSpectrumFactory Factory; 


namespace radmat
{

  namespace L_3F_743_StateFactoryEnv
  {
    namespace 
    {
      template<class T, class U> 
        T* upCast()
        {
          T *t = new U();
          POW2_ASSERT(t); 
          return t; 
        }

      template<class Derived>
        bool do_reg()
        {
          Derived d; 
          bool reg = Factory::Instance().registerObject(
              d.name() , upCast<SpectrumState_p,Derived> ); 
          if( !!! reg )
          {
            std::cout << d.name() << " failed to reg" << std::endl;
          }

          return reg;
        }

      bool do_reg_work()
      {
        bool success = true; 
        success &= do_reg<sPion>();
        success &= do_reg<sPionP>();
        success &= do_reg<sPionPP>();

        success &= do_reg<sRho>();
        success &= do_reg<sRhoP>();
        success &= do_reg<sRhoPP>();
        success &= do_reg<sRhoPPP>();

        success &= do_reg<sA1>(); 
        
        return success; 
      }

      bool local_reg = false; 
    } // anonomyous 

    bool registerAll()
    {
      bool success = true; 
      if( !!! local_reg)
      {
        success &= do_reg_work(); 
        local_reg = true; 
      }
      return success; 
    }

  } // L_3F... 
} // radmat 
