#ifndef MAP_OBJ_H
#define MAP_OBJ_H

#include <map>
#include <vector>


namespace radmat
{
  template<typename K, typename V> 
    struct radMap
    {
      typedef std::map<K,V> data_t;
      typedef typename data_t::const_iterator const_iterator;
      typedef typename data_t::iterator iterator; 
      typedef typename data_t::value_type value_type; 

      bool exists(const K &k) const { return dat.find(k) == dat.end();}
      void insert(const K &k , const V &v) { dat[k] = value_type(k,v); }
      std::vector<K> keys(void) const 
      {
        std::vector<K> v; 
        const_iterator it; 
        for(it = dat.begin(); it != dat.end(); ++it)
          v.push_back(it->first);
        return v;  
      }


      data_t dat;
    };

}



#endif /* MAP_OBJ_H */
