#ifndef ROTATION_GROUP_GENERATOR_H
#define ROTATION_GROUP_GENERATOR_H 




/*
 *
 *  Roughly this guy picks out canonical frames
 *  for evaluating matrix elements in
 *
 *  when we go through and evaluate helicity 
 *  matrix elements the manner in which we 
 *  build them attaches an extra phase
 *  corresponding to a rotation about the
 *  momentum direction of one of the particles
 *
 *  this class picks out which frames we want 
 *  and the canonical rotation classes then
 *  work together to apply a convention 
 *  ( an extra rotation ) to kill the phase
 *
 *  the choice of convention is arbitrary 
 *  the phase only appears when both paricles 
 *  are spin-1 or higher 
 *
 */


#include "rotation_utils.h"
#include "radmat/utils/handle.h"
#include "radmat/utils/stringify.h"
#include "io/adat_xmlio.h"
#include "hadron/irrep_util.h"
#include <map>
#include <vector>
#include "adat/singleton.h"
#include "adat/objfactory.h"
#include <complex>
#include "io/adat_xmlio.h"

// #define BREAK_ROTATIONS_REST


namespace radmat
{
  typedef ADATXML::Array<int> mom_t; 


  // a momentum key 
  struct mom_key
  {
    mom_key() {}
    mom_key(const int xx, const int yy, const int zz)
      : x(xx) , y(yy) , z(zz) , isOh(false) 
    {
      if( (x==0) && (y==0) && (z==0) )
      {
        isOh = true; 
        z = 1; 
      }
    }

    mom_key(const ADATXML::Array<int> &p)
      : x(p[0]) , y(p[1]), z(p[2]), isOh(false)
    { 
      if( (x==0) && (y==0) && (z==0) )
      {
        isOh = true; 
        z = 1; 
      }
    }

    mom_t mom() const
    {
      if( isOh )
        return gen_mom<0,0,0>(); 

      ADATXML::Array<int> p(3);
      p[0] = x; 
      p[1] = y; 
      p[2] = z; 
      return p; 
    }

    mom_t sort_mom() const
    {
      ADATXML::Array<int> p(3);
      p[0] = x; 
      p[1] = y; 
      p[2] = z; 
      return p; 
    }

    bool isOh; 

    int x,y,z; 
  }; 

  // a frame key 
  struct mom_pair_key
  {
    mom_pair_key() {}
    mom_pair_key(const int xl, const int yl, const int zl, 
        const int xr, const int yr, const int zr)
      : l(mom_key(xl,yl,zl)) , r(mom_key(xr,yr,zr))
    {}

    mom_pair_key( const mom_key &ll, const mom_key &rr)
      : l(ll) , r(rr)
    { }

    std::pair<mom_t,mom_t> moms() const
    {
      return std::make_pair( l.mom() , r.mom() ); 
    }

    std::pair<mom_t,mom_t> frame_orientation() const
    { return std::make_pair( l.sort_mom(), r.sort_mom() ); }

    mom_key l,r; 
  }; 



  std::ostream & operator<<(std::ostream &o, const mom_key &); 
  std::ostream & operator<<(std::ostream &o, const mom_pair_key &);


  // key comparison class, bung them all together  
  struct mom_key_comp
  {
    bool operator()(const mom_key &l, const mom_key &r) const
    {
      if(l.x != r.x)
        return l.x < r.x; 
      if(l.y != r.y)
        return l.y < r.y; 
      if(l.z != r.z)
        return l.z < r.z; 
      return l.isOh < r.isOh; 
    }

    bool operator()(const mom_pair_key &l, const mom_pair_key &r) const
    {
      if ( !!! exact_equivalence(l.l,r.l) )
        return this->operator()(l.l,r.l);

      return this->operator()(l.r,r.r); 
    }

    bool exact_equivalence( const mom_key &l, const mom_key &r) const 
    {
      return ((l.x==r.x) && (l.y==r.y) && (l.z==r.z) && (l.isOh==r.isOh));
    }

    bool exact_equivalence(const mom_pair_key &l , const mom_pair_key &r) const
    {
      return exact_equivalence(l.l,r.l) && exact_equivalence(l.r,r.r); 
    }

  };





  template<int MOM_MAX, int INSERTION_MAX_SQ>
    struct RotationGroupGenerator
    {

      bool is_related_by_rotation( const mom_pair_key &test, const mom_pair_key &can) const
      {
#ifdef BREAK_ROTATIONS_REST
        
        if( test.l.isOh != can.l.isOh ) 
          return false; 

        if( test.r.isOh != can.r.isOh )
          return false; 

        return related_by_rotation( test.l.sort_mom() , test.r.sort_mom() , can.l.sort_mom(), can.r.sort_mom() ); 
#else
        return related_by_rotation( test.l.mom() , test.r.mom() , can.l.mom(), can.r.mom() ); 
#endif
      }


      void insert(const mom_pair_key &k)
      {
        bool found = false; 
        std::vector<mom_pair_key>::const_iterator it; 
        for(it = canonical_frames.begin(); it != canonical_frames.end(); ++it)
          if( is_related_by_rotation( *it , k ) )
            break;

        if( it == canonical_frames.end() ) 
        {
          can_frame_map.insert(std::make_pair(k,k)); 
          canonical_frames.push_back(k); 
        }
        else
          can_frame_map.insert(std::make_pair(k,*it)); 

      }


      void initialize_frames(void)
      {
        // the loop ordering here sets the canonical names, 
        // use positive then lower so it is human readable
        for(int i = MOM_MAX; i >= -MOM_MAX; --i)
          for(int j = MOM_MAX; j >= -MOM_MAX; --j)
            for(int k = MOM_MAX; k >= -MOM_MAX; --k)
              for(int l = MOM_MAX; l >= -MOM_MAX; --l)
                for(int m = MOM_MAX; m >= -MOM_MAX; --m)
                  for(int n = MOM_MAX; n >= -MOM_MAX; --n)
                  {
                    if( (i-l)*(i-l) + (j-m)*(j-m) + (k-n)*(k-n) > INSERTION_MAX_SQ )
                      continue;
                    if( i*i + j*j + k*k > MOM_MAX * MOM_MAX )
                      continue; 
                    if( l*l + m*m + n*n > MOM_MAX * MOM_MAX )
                      continue; 

                    // flip ordering from loops so that the p_ref get preferentially selected 
                    insert( mom_pair_key(  k,j,i, n,m,l ) ); 
                  }
      }

      template<typename T> 
        void do_exit(const std::string &s, const T &t) const
        {
          std::cout << __PRETTY_FUNCTION__ << s << t << std::endl;
          exit(1); 
        }


      bool registerAll(void)
      {
        initialize_frames(); 
        return true; 
      }

      std::vector<mom_pair_key> unique_frames() const 
      { return canonical_frames; }

      std::pair<mom_t,mom_t> canonical_frame_moms(const mom_pair_key &p) const
      {
        mom_pair_key k = canonical_frame(p); 
        return k.moms(); 
      }

      mom_pair_key canonical_frame(const mom_pair_key &p) const
      {
        std::map<mom_pair_key,mom_pair_key,mom_key_comp>::const_iterator it; 
        it = can_frame_map.find(p);
        if( it == can_frame_map.end() )
          do_exit( "unknown key" , p ); 

        return it->second; 
      }

      mom_pair_key canonical_frame(const mom_t &l, const mom_t &r) const 
      {
        return canonical_frame( mom_pair_key( mom_key(l) , mom_key(r) ) ); 
      }

      std::pair<mom_t,mom_t> canonical_frame_moms(const mom_t &l, const mom_t &r) const
      {
        return canonical_frame_moms( mom_pair_key( mom_key(l) , mom_key(r)) ); 
      }


      typedef std::pair<mom_pair_key, mom_pair_key> data_type;  


      // key is frame label value is canonical frame
      std::map<mom_pair_key,mom_pair_key,mom_key_comp> can_frame_map; 
      std::vector<mom_pair_key> canonical_frames; 
    };



  namespace LatticeRotationEnv
  {
    // returns the canonical frame 
    mom_pair_key
      rotation_group_key(const mom_t &l, const mom_t &r);

    std::pair<mom_t,mom_t>
      rotation_group_can_mom(const mom_t &l, const mom_t &r);

    std::string 
      rotation_group_label(const mom_t &l, const mom_t &r);

    typedef Util::SingletonHolder<RotationGroupGenerator<2,4> > 
      TheRotationGroupGenerator;

    bool registerAll(); 
  }



} // radmat 





#endif /* ROTATION_GROUP_GENERATOR_H */
