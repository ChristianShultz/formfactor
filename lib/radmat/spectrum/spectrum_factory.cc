/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : spectrum_factory.cc

 * Purpose :

 * Creation Date : 06-05-2014

 * Last Modified : Tue 06 May 2014 07:14:15 PM EDT

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/

#include "spectrum_factory.h"
#include "spectrum_3F_743.h"

typedef radmat::TheSpectrumFactory Factory; 

namespace radmat
{

  namespace 
  {

    bool do_reg_work(const int enumLat)
    {
      bool success = true; 

      switch( enumLat )
      {
        case LATTICE::L_3F_743:
          success &= L_3F_743_StateFactoryEnv::registerAll();

        default :
          success &= false; 
      }

      return success; 
    }


    bool local_reg = false; 
  } // anonomyous 

  namespace TheSpectrumFactoryEnv
  {

    bool registerAll( const int enumLat)
    {
      bool success = true; 

      if( !!! local_reg)
      {
        success &= do_reg_work(enumLat); 
        local_reg = true; 
      }

      return success; 
    }


    // make it blow up if anything goes wrong by wrapping another 
    // call around the factory.createObj method
    rHandle<SpectrumState_p > callFactory(const std::string &id)
    {
      SpectrumState_p *foo;
      foo = NULL;
      try
      {
        foo = Factory::Instance().createObject(id);
      }
      catch(std::exception &e)
      {
        std::cout << "id - " << id << std::endl;
        std::cout << __PRETTY_FUNCTION__ << e.what(); 
        throw e; 
      }
      catch(std::string &s)
      {
        std::cout << "id - " << id << std::endl;
        std::cout << __PRETTY_FUNCTION__ << s << std::endl;
        throw s;
      }
      catch(...)
      {
        std::cout << "id - " << id << std::endl;
        std::cout << __PRETTY_FUNCTION__ << ": some error" << std::endl;
        POW2_ASSERT(false); 
      }

      // not a null pointer
      POW2_ASSERT(foo);
      return rHandle<SpectrumState_p >(foo);
    }

    std::pair<double,double> pole_enhancement(const std::string &l, 
        const ADATXML::Array<int> &lmom, 
        const std::string &r, 
        const ADATXML::Array<int> &rmom, 
        const std::string &pole)
    {
      rHandle<SpectrumState_p> lefty,righty,poley; 
      lefty = callFactory(l);
      righty = callFactory(r);
      poley = callFactory(pole); 

#if 1
      std::cout << "left mass " << lefty->rest_mass()
        << "\nright mass " << righty->rest_mass() << std::endl;

      std::cout << "boost left mass " << lefty->boost_mass(lmom)
        << "\nboost right mass " << righty->boost_mass(rmom) << std::endl;

      std::cout << "pole mass " << poley->rest_mass() << std::endl;
      std::cout << "pole mass^2 " << (poley->rest_mass())*(poley->rest_mass()) << std::endl;
#endif 

      double q2 = Mink_qsq( lmom, lefty->rest_mass(), 
          rmom, righty->rest_mass(), 
          righty->mom_fac());

      double polish = q2 + (poley->rest_mass())*(poley->rest_mass()); 

      return std::make_pair(1./polish,q2); 
    }

  }

}; 



