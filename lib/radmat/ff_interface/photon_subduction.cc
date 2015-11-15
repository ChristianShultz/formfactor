/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

* File Name : photon_subduction.cc

* Purpose :

* Creation Date : 17-03-2014

* Last Modified : Wed 23 Apr 2014 02:31:17 PM EDT

* Created By : shultz

_._._._._._._._._._._._._._._._._._._._._.*/

#include "photon_subduction.h"
#include "radmat/utils/stringify.h"
#include "radmat/redstar_interface/redstar_cartesian_interface.h"
#include "radmat/redstar_interface/redstar_photon_props.h"
#include "formfactor_subduced_formfactors.h"
#include "radmat/data_representation/data_representation.h"
#include "itpp/itbase.h"



namespace radmat
{
  namespace
  {
    
  //////////////////////////////////////////////////////////////////
  lorentz_to_cubic_t 
    operator*(const std::complex<double> &c,
        const lorentz_to_cubic_t &l)
    {
      return c*l; 
    }

  //////////////////////////////////////////////////////////////////
  itpp::Vec<lorentz_to_cubic_t> 
    operator*(const std::complex<double> &d,
        const itpp::Vec<lorentz_to_cubic_t> &o)
    {
      itpp::Vec<lorentz_to_cubic_t> foo(o.size()); 
      for(int i = 0; i < o.size(); ++i)
        foo[i] = d * o[i];
      return foo; 
    }



    //////////////////////////////////////////////////////////////////
    itpp::Vec<lorentz_to_cubic_t> 
    operator*(const itpp::Mat<std::complex<double> > &m,
        const itpp::Vec<lorentz_to_cubic_t> &v)
    {
      POW2_ASSERT(v.size() == m.cols()); 
      itpp::Vec<lorentz_to_cubic_t> ret(m.rows()); 

      for(int row = 0; row < m.rows(); ++row)
        for(int col = 0; col < m.cols(); ++col)
        {
          if(std::norm(m(row,col)) > 0.0001)
          {
            // init 
            if(ret[row].m_expr.size() == 0)
              ret[row] =  m(row,col)*v(col);
            else
              ret[row] = ret[row] + m(row,col)*v(col);
          }
        }

      return ret; 
    }

    


    lorentz_to_cubic_t
      generate_list(const std::complex<double> & d , const int i )
      {
        return convert_to_list( lorentz_to_cubic_t::ListObj_t(d,i) ); 
      }



    lorentz_to_cubic_t
      spatial_insertions(const std::string &rep, 
          const int row, 
          const ADATXML::Array<int> &q)
      {
        if(rep == Stringify<J1m>() )
          return generate_list(std::complex<double>(1.,0.),row); 
       
        std::string key = Stringify<J1m>() + rep;  
        bool create = PHOTON_CREATE;
        
        itpp::Mat<std::complex<double> > eps = eps3d( q, create );  
        itpp::Vec<lorentz_to_cubic_t> cart(3),hel; 
        std::complex<double> one(1.,0.); 
        cart[0] = generate_list(one,1);
        cart[2] = generate_list(one,2);
        cart[3] = generate_list(one,3);

        // 0->+1 , 1->0 , 2->-1
        hel = eps * cart; 

        // determine if we have a valid key 
        SubduceTableMap::map_t::const_iterator subduce; 
        subduce = TheSmarterSubduceTableMap::Instance().mappy.find(key); 
        if( subduce == TheSmarterSubduceTableMap::Instance().mappy.end())
        {
          std::cout << __PRETTY_FUNCTION__ << ": error, key " << key 
            << " is not present" << std::endl;
          throw std::string("subduction error"); 
          exit(1); 
        }

        // pull down the subduction list 
        const SubduceTableMap::irrep_sub_table* sub = subduce->second; 
        SubduceTableMap::sub_list::const_iterator it; 

        // these guys are rows but 0 based so subtract 1 from the 
        // FORTRAN arrays that adat uses 
        lorentz_to_cubic_t ret; 
        for(it = sub->begin(row); it != sub->end(row); ++it)
          if( create )
            ret = ret + it->first * hel[1 - it->second];  
          else
            ret = ret + std::conj(it->first) * hel[1 - it->second];  

        // this is now a sum over cartesian indicies that takes us to 
        // a some row of some irrep of some cubic group  
        return ret; 
      }

    lorentz_to_cubic_t
      temporal_insertions(const std::string &rep, 
          const int row, 
          const ADATXML::Array<int> &q)
      {
        if(rep == Stringify<J0p>() )
          return generate_list(std::complex<double>(1.,0.),row); 

        std::string key = Stringify<J0p>() + rep;  
        bool create = PHOTON_CREATE;
        
        itpp::Mat<std::complex<double> > eps = eps3d( q, create );  

        lorentz_to_cubic_t ret,scal; 
        scal = generate_list(std::complex<double>(1.,0.),4); 

        // determine if we have a valid key 
        SubduceTableMap::map_t::const_iterator subduce; 
        subduce = TheSmarterSubduceTableMap::Instance().mappy.find(key); 
        if( subduce == TheSmarterSubduceTableMap::Instance().mappy.end())
        {
          std::cout << __PRETTY_FUNCTION__ << ": error, key " << key 
            << " is not present" << std::endl;
          throw std::string("subduction error"); 
          exit(1); 
        }

        // pull down the subduction list 
        const SubduceTableMap::irrep_sub_table* sub = subduce->second; 
        if(sub->size() != 1)
        {
          std::cout << __PRETTY_FUNCTION__ << ": the table is the wrong dimension"
            << " I think this ( " << key << " ) is a scalar key but the table has"
            << " dimension " << sub->size() << std::endl;
          throw std::string("subduction error"); 
          exit(1); 
        }

        std::complex<double> weight = sub->begin(1)->first; 
        
        if(create)
          ret = weight * scal; 
        else
          ret = std::conj(weight) * scal; 
    
        return ret; 
      }
  }


  lorentz_to_cubic_t
    photon_subduction(const std::string &rep, 
        const std::string &spher_rep, 
        const int row, 
        const ADATXML::Array<int> &q)
    {
      typedef lorentz_to_cubic_t (*ptr)(const std::string &,
          const int, 
          const ADATXML::Array<int> &);


      std::map<std::string,ptr> options_map; 
      options_map.insert(std::make_pair(Stringify<J0p>(),&temporal_insertions));
      options_map.insert(std::make_pair(Stringify<J1m>(),&spatial_insertions));

      std::map<std::string,ptr>::const_iterator foo; 
      foo = options_map.find(spher_rep); 
      if( foo == options_map.end() )
      {
        std::cout << __PRETTY_FUNCTION__ << ": error, " << spher_rep 
          << " is unsupported, try one of the following" << std::endl;
        for(foo = options_map.begin(); foo != options_map.end(); ++foo)
          std::cout << foo->first << std::endl;

        exit(1); 
      }
      
      ptr p = foo->second; 

      return (*p)(rep,row,q); 
    }


}
