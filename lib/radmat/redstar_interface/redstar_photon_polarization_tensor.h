#ifndef REDSTAR_PHOTON_POLARIZATION_TENSOR_H
#define REDSTAR_PHOTON_POLARIZATION_TENSOR_H 



#include "radmat/utils/tensor.h"
#include "radmat/utils/aux.h"
#include "radmat/utils/handle.h"
#include "radmat/utils/polarization_tensors.h"
#include "ensem/ensem.h"
#include "hadron/irrep_util.h"
#include "hadron/clebsch.h"
#include "itpp/itbase.h"
#include "xml_array.h"
#include <utility>
#include <map>



/*
NB: These rotate like Helicity creation operators or equivalently 
like helicity states using the rotation conventions from ADAT.
 */



namespace radmat
{


  /////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////
  // 3 Dimensional version

  // fwd of a recursive template pattern for generating polarisation tensors
  template<idx_t J>
    struct genPolTens3D;

  /*
     A genPolTens is a recursive way of generating polarisation tensors in a
     mildly efficient maner.  Use of handles allows sharing of data between the 
     parent instance (rank N) and its children (rank N-1 .. 1) so that we only have
     to generate the intermediate tensors one time at most.
   */

  // J = 0 specialization
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  template<>
    struct genPolTens3D<0>
    {  }; // empty

  // J = 1 specialization
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////


  // gcc is frustrating .. it won't let a many to one friendship pattern without
  // totally rewriting the templates so I'm leaving genPolTens<1> completly public
  template<>
    struct genPolTens3D<1>
    {
      // save some typing
      typedef XMLArray::Array<int> mom_t;
      typedef std::map<pKey_t,TensorBase*> map_t;
      typedef rHandle<map_t> map_handle;
      typedef rHandle<mom_t> mom_handle;

      genPolTens3D(void); // hidden

      //! constructor taking momentum -- tells us how to do the rotations via adat
      genPolTens3D(const mom_t &_p)
        : h_map(new map_t()) , h_mom(new mom_t(_p))
      { }

      ~genPolTens3D(void) { refresh_map(); }

      //! get a polarisation tensor for the inputs in the direction of p
      Tensor<std::complex<double>,1> operator()(const short hel)
      {
        return make(hel);
      }

      protected:
      //!  do not use!
      genPolTens3D(map_handle &_h_map, mom_handle & _h_mom )
        : h_map(_h_map) , h_mom(_h_mom)
      { }

      //! do not use!
      Tensor<std::complex<double> , 1> get(const short hel)
      {
        map_t::const_iterator it;
        it = h_map->find(pKey_t(1,hel));
        if(it != h_map->end())
          return m_downcast(it->second);

        h_map->insert(map_t::value_type(pKey_t(1,hel),make(hel).clone()));
        return get(hel);
      }

      //! do not use!
      Tensor<std::complex<double> , 1> m_downcast(TensorBase* bar)
      {
        Tensor<std::complex<double> , 1> * foo = dynamic_cast<Tensor<std::complex<double>, 1 >* >(bar);
        POW2_ASSERT(foo);
        return *foo;
      }

      // clean up the map each time we get a new input
      void refresh_map(void)
      {
        map_t::iterator it;
        for(it = h_map->begin(); it != h_map->end(); it++)
        {
          delete it->second;
          it->second = NULL;
        }
        h_map->clear();
      }

      //! make a J = 1 polarisation tensor
      Tensor<std::complex<double>,1> make( const short hel)
      {
        return creation_op_J1_3(*h_mom,hel); 
      }

      friend struct genPolTens3D<2>;

      // shared data store    
      map_handle h_map;
      mom_handle h_mom;
    };

  // recursive defs
  ////////////////////////////////////////////////////////////////////////////////////////////////////////

  template<idx_t J>
    struct genPolTens3D
    {
      typedef XMLArray::Array<int> mom_t;
      typedef std::map<pKey_t,TensorBase*> map_t;
      typedef rHandle<map_t> map_handle;
      typedef rHandle<mom_t> mom_handle;

      genPolTens3D(void); // hidden

      //! constructor taking momentum -- tells us how to do the rotations via adat
      genPolTens3D(const mom_t &_p)
        : h_map(new map_t()) , h_mom(new mom_t(_p))
      { }

      ~genPolTens3D(void) 
      {
        refresh_map(); 
      }

      //! return a polarisation tensor of rank J for the inputs
      Tensor<std::complex<double>, J> operator()(const short hel)
      {
        refresh_map();
        return get(hel);
      }

      private:
      // clean up the map each time we get a new input
      void refresh_map(void)
      {
        map_t::iterator it;
        for(it = h_map->begin(); it != h_map->end(); it++)
        {
          delete it->second;
          it->second = NULL;
        }
        h_map->clear();
      }

      // wrap make in a sensible form
      Tensor<std::complex<double> , J> get(const short hel)
      {
        map_t::const_iterator it;
        it = h_map->find(pKey_t(J,hel));
        if(it != h_map->end())
          return m_downcast(it->second);

        make(hel);
        return get(hel);
      }

      // downcast the base pointers to the derived type and return an actual obj
      Tensor<std::complex<double>, J> m_downcast(TensorBase * bar)
      {
        Tensor<std::complex<double> , J> * foo = dynamic_cast<Tensor<std::complex<double>, J >* >(bar);
        POW2_ASSERT(foo);
        return *foo;
      }

      // couple them together using the SU(2) CG found in adat
      void make(const short hel)
      {
        genPolTens3D< J - 1 > Jm(h_map, h_mom);
        genPolTens3D<1> J1(h_map, h_mom);
        double factor;

        Tensor<std::complex<double> , J> ptensor(std::vector<idx_t>(J, 3), 0.);

        for(short big_h = -(J-1); big_h < J ; ++big_h)
          for(short small_h = -1; small_h < 2; ++small_h)
          {
            factor = Hadron::clebsch(2 * (J - 1), 2 * big_h, 2, 2 * small_h, 2 * J, 2 * hel);

            if(factor != 0.)
              ptensor += (factor * (Jm.get(big_h) ^ J1.get(small_h)));
          }

        h_map->insert(map_t::value_type(pKey_t(J, hel), ptensor.clone()));
      }

      //! private constructor to share parent's information with children 
      genPolTens3D(map_handle &_h_map, mom_handle &_h_mom)
        : h_map(_h_map)  , h_mom(_h_mom)
      { }

      friend struct genPolTens3D<J+1>;

      map_handle h_map;
      mom_handle h_mom;
    };

} // radmat



#endif /* REDSTAR_PHOTON_POLARIZATION_TENSOR_H */
