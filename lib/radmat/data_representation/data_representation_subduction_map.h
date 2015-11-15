#ifndef DATA_REPRESENTATION_SUBDUCTION_MAP_H
#define DATA_REPRESENTATION_SUBDUCTION_MAP_H 

#include "data_representation_lorentz_groups.h"
#include "data_representation_cubic_groups.h"
#include "radmat/utils/pow2assert.h"

#include <map>
#include <vector>
#include <string>

#include "adat/singleton.h"


namespace radmat
{

  // the one map 
  struct SubduceTableMap
  {
    typedef std::pair<std::complex<double> , int > sub_chunk;
    typedef std::vector<sub_chunk> sub_list; 

    // hold the cont representation, the lattice representation, 
    // and the table that subduces from cont to lat
    //    -- build the subduction info into the table radmat uses 
    struct irrep_sub_table
    {
      typedef std::vector<sub_list> sub_table; 
      typedef rHandle<LorentzRep> cont_rep;
      typedef rHandle<CubicRep> lat_rep; 

      irrep_sub_table( const sub_table &s, const cont_rep &c , const lat_rep &l)
        : sub(s) , cont(c) , lat(l)
      { }


      // maps are 1 based to the outside world 
      // return a const reference to our member so we can
      // use begin and end functions w/o copy and in place
      const sub_list& query_1(const int row) const
      {
        POW2_ASSERT( row > 0 ); 
        return sub.at(row-1); 
      }

      sub_table::const_iterator begin() const { return sub.begin(); }
      sub_table::const_iterator end() const { return sub.end(); }

      sub_list::const_iterator begin(const int row) const 
      { return query_1(row).begin(); }
      sub_list::const_iterator end(const int row) const 
      { return query_1(row).end(); }

      int size() const { return sub.size(); }
      

      cont_rep get_c_rep() const { return rHandle<LorentzRep>(cont->clone()); }
      lat_rep get_l_rep() const { return rHandle<CubicRep>(lat->clone()); }

      private: 

      sub_table sub; 
      cont_rep cont; 
      lat_rep lat; 
    };

    const irrep_sub_table* get_table(const std::string &s)
    {
      map_t::const_iterator it; 
      it = mappy.find(s); 
      if(it == mappy.end())
      {
        std::cout << __PRETTY_FUNCTION__ << " key " << s << " not present" << std::endl;
        for(it = mappy.begin(); it != mappy.end(); ++it)
          std::cout << it->first << std::endl;

        throw std::string("subduce table key not present"); 
      }

      return it->second; 
    }

    typedef std::map<std::string,irrep_sub_table*> map_t; 

    ~SubduceTableMap() 
    {
      map_t::iterator it; 
      for(it = mappy.begin(); it != mappy.end(); ++it)
        delete it->second; 
    }

    map_t mappy; 
  };


  typedef Util::SingletonHolder< SubduceTableMap > TheSmarterSubduceTableMap; 

  std::string make_subduce_table_map_id(const rHandle<LorentzRep> &,
      const rHandle<CubicRep> &);

  namespace TheSmarterSubduceTableMapFactoryEnv
  {
    std::vector<std::string> all_keys(); 
    bool registerAll(); 
  }

} // radmat

#endif /* DATA_REPRESENTATION_SUBDUCTION_MAP_H */
