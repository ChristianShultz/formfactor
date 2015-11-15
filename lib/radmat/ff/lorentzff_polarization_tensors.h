#ifndef LORENTZFF_POLARIZATION_TENSORS_H
#define LORENTZFF_POLARIZATION_TENSORS_H 


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




namespace radmat
{


  /////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////
  // 4 Dimensional version

  // fwd of a recursive template pattern for generating polarisation tensors
  template<idx_t J>
    struct ZAxisHelicityTensor;

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
    struct ZAxisHelicityTensor<0>
    { }; // empty

  // J = 1 specialization
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////


  // gcc is frustrating .. it won't let a many to one friendship pattern without
  // totally rewriting the templates so I'm leaving genPolTens<1> completly public
  template<>
    struct ZAxisHelicityTensor<1>
    {
      typedef std::map<pKey_t,TensorBase*> map_t;
      typedef rHandle<map_t> map_handle;

      ZAxisHelicityTensor(void)
        : h_map(new map_t()) 
      { }

      virtual ~ZAxisHelicityTensor(void) { refresh_map(); }

      //! get a polarisation tensor for the inputs in the direction of p
      Tensor<std::complex<double>,1> operator()(const double mod_p,
          const short hel,
          const double E)
      {
        return make(mod_p,hel,E);
      }

      protected:
      ZAxisHelicityTensor(map_handle &_h_map)
        : h_map(_h_map) 
      { }

      //! do not use!
      Tensor<std::complex<double> , 1> get(const double mod_p, 
          const short hel,
          const double E)
      {
        map_t::const_iterator it;
        it = h_map->find(pKey_t(1,hel));
        if(it != h_map->end())
          return m_downcast(it->second);

        h_map->insert(map_t::value_type(pKey_t(1,hel),make(mod_p,hel,E).clone()));
        return get(mod_p,hel,E);
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
      Tensor<std::complex<double>,1> make(const double mod_p, 
          const short hel,
          const double E)
      {
        return creation_op_J1_4z(mod_p,hel,E); 
      }

      friend struct ZAxisHelicityTensor<2>;
      friend struct ZAxisHelicityTensor<3>;
      friend struct ZAxisHelicityTensor<4>;
      friend struct ZAxisHelicityTensor<5>;
      friend struct ZAxisHelicityTensor<6>;
      friend struct ZAxisHelicityTensor<7>;
      friend struct ZAxisHelicityTensor<8>;
      friend struct ZAxisHelicityTensor<9>;
      friend struct ZAxisHelicityTensor<10>;

      // shared data store    
      map_handle h_map;
    };

  // recursive defs
  ////////////////////////////////////////////////////////////////////////////////////////////////////////

  template<idx_t J>
    struct ZAxisHelicityTensor
    {
      typedef std::map<pKey_t,TensorBase*> map_t;
      typedef rHandle<map_t> map_handle;

      ZAxisHelicityTensor(void)
        : h_map(new map_t()) 
      { }

      ~ZAxisHelicityTensor(void) 
      {
        refresh_map(); 
      }

      //! return a polarisation tensor of rank J for the inputs
      Tensor<std::complex<double>, J> operator()(const double mod_p,
          const short hel,
          const double E)
      {
        refresh_map();
        return get(mod_p,hel,E);
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
      Tensor<std::complex<double> , J> get(const double mod_p,
          const short hel,
          const double E)
      {
        map_t::const_iterator it;
        it = h_map->find(pKey_t(J,hel));
        if(it != h_map->end())
          return m_downcast(it->second);

        make(mod_p,hel,E);
        return get(mod_p,hel,E);
      }

      // downcast the base pointers to the derived type and return an actual obj
      Tensor<std::complex<double>, J> m_downcast(TensorBase * bar)
      {
        Tensor<std::complex<double> , J> * foo = dynamic_cast<Tensor<std::complex<double>, J >* >(bar);
        POW2_ASSERT(foo);
        return *foo;
      }

      // couple them together using the SU(2) CG found in adat
      void make(const double mod_p,
          const short hel,
          const double E)
      {
        ZAxisHelicityTensor< J - 1 > Jm(h_map);
        ZAxisHelicityTensor<1> J1(h_map);
        double factor;

        Tensor<std::complex<double> , J> ptensor(std::vector<idx_t>(J, 4), 0.);

        for(short big_h = -(J-1); big_h < J ; ++big_h)
          for(short small_h = -1; small_h < 2; ++small_h)
          {
            factor = Hadron::clebsch(2 * (J - 1), 2 * big_h, 2, 2 * small_h, 2 * J, 2 * hel);

            if(factor != 0.)
              ptensor += (factor * (Jm.get(mod_p,big_h,E) ^ J1.get(mod_p,small_h,E)));
          }

        h_map->insert(map_t::value_type(pKey_t(J, hel), ptensor.clone()));
      }

      //! private constructor to share parent's information with children 
      ZAxisHelicityTensor(map_handle &_h_map)
        : h_map(_h_map) 
      { }

      friend struct ZAxisHelicityTensor<J+1>;

      map_handle h_map;
    };

} // radmat






#endif /* LORENTZFF_POLARIZATION_TENSORS_H */
