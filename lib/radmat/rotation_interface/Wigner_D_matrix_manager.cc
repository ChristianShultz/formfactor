/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : Wigner_D_matrix_manager.cc

 * Purpose :

 * Creation Date : 14-04-2014

 * Last Modified : Fri 25 Jul 2014 01:06:47 PM EDT

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/



#include "Wigner_D_matrix_manager.h"
#include "hadron/clebsch.h"
#include "radmat/utils/pow2assert.h"
#include "semble/semble_semble.h"
#include "rotation_group_generator.h"
#include "radmat/utils/printer.h"
#include "Wigner_D_matrix_factory.h"


namespace radmat
{
  namespace
  {


    std::pair<mom_t,mom_t> bb_rot_rest(const mom_t &l, const mom_t &r)
    {
      std::pair<mom_t,mom_t> moms; 

#ifdef BREAK_ROTATIONS_REST
      mom_pair_key ky = mom_pair_key( mom_key(l) , mom_key(r) ); 
      moms =  ky.frame_orientation(); 
#else 
      moms = std::make_pair(l,r);  
#endif 

      return moms; 
    }




    std::complex<double> complex_zero(0.,0.); 

    std::complex<double> round_to_zero(const std::complex<double> &cd, const double thresh=1e-6)
    {return itpp::round_to_zero(cd,thresh); }

    std::ostream & operator<<(std::ostream &o, const mom_t &p)
    { return ( o << "p" << p[0] << p[1] << p[2]);}

    // a momentum key 
    struct wig_key
    {
      wig_key() {}
      wig_key(const int xx, const int yy, const int zz)
        : x(xx) , y(yy) , z(zz)
      { }

      wig_key(const ADATXML::Array<int> &p)
        : x(p[0]) , y(p[1]), z(p[2]) 
      { }

      int x,y,z; 
    }; 

    // error streaming 
    std::ostream & operator<<(std::ostream &o, const wig_key &k)
    {return (o << "p" << k.x << k.y << k.z);}

    // a frame key 
    struct wig_pair_key
    {
      wig_pair_key() {}
      wig_pair_key(const ADATXML::Array<int> &ll , 
          const ADATXML::Array<int> &rr, 
          const int JJ)
        : l(ll) , r(rr), J(JJ)
      {}

      wig_key l,r;
      int J;  
    }; 

    // error streaming 
    std::ostream &operator<<(std::ostream &o, const wig_pair_key &k)
    { return (o << "l" << k.l << "r" << k.r << "J" << k.J);}


    // key comparison class, bung them all together  
    struct wig_key_comp
    {
      bool operator()(const wig_key &l, const wig_key &r) const
      {
        if(l.x != r.x)
          return l.x < r.x; 
        else if(l.y != r.y)
          return l.y < r.y; 
        else
          return l.z < r.z; 
      }

      bool operator()(const wig_pair_key &l, const wig_pair_key &r) const
      {
        if ( !!! exact_equivalence(l.l,r.l) )
          return this->operator()(l.l,r.l);
        else if( !!! exact_equivalence(l.r,r.r) )
          return this->operator()(l.r,r.r); 
        else
          return l.J < r.J; 
      }

      bool exact_equivalence( const wig_key &l, const wig_key &r) const 
      {
        return ((l.x==r.x) && (l.y==r.y) && (l.z==r.z));
      }

      bool exact_equivalence(const wig_pair_key &l , const wig_pair_key &r) const
      {
        return exact_equivalence(l.l,r.l) && exact_equivalence(l.r,r.r) && (l.J == r.J); 
      }

    };


    typedef std::map<wig_pair_key,WignerMatrix_t,wig_key_comp> WignerMatrixMap_t; 
    WignerMatrixMap_t left_map, right_map; 

    WignerMatrixMap_t* call_left_map()
    {
      return &left_map;
    }

    WignerMatrixMap_t* call_right_map()
    {
      return &right_map;
    }

    struct injection_map_printer
    {
      static void print(const std::string &msg)
      {}
      // {std::cout << "injection_map_printer " << msg << std::endl;}
    };

    // pre compute avaliable wigner matricies  
    void inject_map(const int J)
    {
      typedef radmat::LatticeRotationEnv::TheRotationGroupGenerator RG; 
      std::map<mom_pair_key,mom_pair_key,mom_key_comp>::const_iterator it; 

      DMatrixManager Wigner; 

      for(it = RG::Instance().can_frame_map.begin(); it != RG::Instance().can_frame_map.end(); ++it)
      {
        mom_t l, r; 
        l = it->first.l.mom(); 
        r = it->first.r.mom(); 
        rCompEulAngles eul = Wigner.rotation_matrix(l,r); 
        WignerMatrix_t * left = Wigner.left_wigner_matrix(eul,l,r,J); 
        WignerMatrix_t *right = Wigner.right_wigner_matrix(eul,l,r,J); 

        wig_pair_key key(l,r,J); 

        //  std::stringstream ss; 
        //  ss << key; 
        //  printer_function<injection_map_printer>(ss.str()); 

        call_left_map()->insert(std::make_pair(key,*left)); 
        call_right_map()->insert(std::make_pair(key,*right)); 

        delete left; 
        delete right; 
      }
    }




    bool local_registration = false;  
    bool use_wigner_map = false; 

    bool init_maps(const int J)
    {
      if( !!! local_registration )
      {
        for(int i = 0; i <= J; ++i)
          inject_map(i); 

        local_registration = true; 
        use_wigner_map = true; 
      }
      return true; 
    }


    mom_t mul_mom(const RotationMatrix_t *R,
        const mom_t &x)
    {
      mom_t chk = gen_mom<0,0,0>();

      for(int i = 0; i < 3; ++i)
      {
        double res(0.); 
        for(int j = 0; j < 3; ++j)
          if ( fabs( (*R)[i+1][j+1] ) > 1e-6 )
            res += (*R)[i+1][j+1] * double(x[j]); 

        chk[i] = (res > 0.) ? int( fabs(res) + 0.5 ) : - int( fabs(res) + 0.5 ); 
      }

      return chk;  
    } 


  } // anonomyous 


  mom_pair_key
    DMatrixManager::get_frame(const mom_t &l, const mom_t &r) const
    {
      return radmat::LatticeRotationEnv::rotation_group_key(l,r); 
    }

  // there is a delta function that MUST be satisfied 
  void 
    DMatrixManager::check_throw_frame_err(const rCompEulAngles &eul, 
        const std::pair<mom_t,mom_t> &f, 
        const std::pair<mom_t,mom_t> &c) const
    {
      RotationMatrix_t *R( new RotationMatrix_t( genRotationMatrix(eul) ) ); 

      if( !!! check_total_frame_transformation(R,f.first,f.second,c.first,c.second,true) )
      {
        std::cout << __func__ << ": throwing string " << std::endl;
        delete R; 
        throw std::string("triad wigner rotation error"); 
      }

      delete R; 
    }

  WignerMatrix_t* 
    DMatrixManager::get_can_mat(const mom_t &p, const int J) const
    {
      return radmat::WignerDMatrixEnv::call_factory(p,J); 
    }  

  void 
    DMatrixManager::conjugate(WignerMatrix_t * D) const
    {
      WignerMatrix_t::iterator it;
      for(it = D->begin(); it != D->end(); ++it)
        *it = std::conj(*it); 
    }

  void 
    DMatrixManager::transpose(WignerMatrix_t *D) const
    {
      WignerMatrix_t foo(*D); 
      std::vector<idx_t> dimensions = D->getDim(); 
      POW2_ASSERT( dimensions.size() == 2 );
      POW2_ASSERT( dimensions[0] == dimensions[1] );
      int bound = dimensions[0]; 
      for(int i = 0; i < bound; ++i)
        for(int j =0; j < bound; ++j)
          (*D)[i][j] = foo[j][i]; 
    }

  void 
    DMatrixManager::dagger(WignerMatrix_t *D) const
    {
      conjugate( D ); 
      transpose( D ); 
    }


  void 
    DMatrixManager::clean(WignerMatrix_t *D, const double thresh) const
    {
      WignerMatrix_t::iterator it; 
      for(it = D->begin(); it != D->end(); ++it)
        *it = round_to_zero( *it , thresh ); 
    }

  rCompEulAngles 
    DMatrixManager::rotation_matrix(const mom_t &l, const mom_t &r) const
    {
      mom_pair_key frame = get_frame(l,r); 
      std::pair<mom_t,mom_t> f; 
#ifdef BREAK_ROTATIONS_REST 
      f = frame.frame_orientation(); 
#else
      f = frame.moms();
#endif
      std::pair<mom_t,mom_t> moms = bb_rot_rest(l,r); 

      return generate_rotation_matrix(moms.first,moms.second,f.first,f.second); 
    }

  WignerMatrix_t*
    DMatrixManager::wigner_matrix(const rEulAngles &eul,
        const int J) const
    {
      int bound = 2*J+1; 
      WignerMatrix_t *W  = new WignerMatrix_t( (TensorShape<2>())[bound][bound] , std::complex<double>(0.,0.));   

      for(int m1 = -J; m1 <= J; ++m1)
        for(int m2 = -J; m2 <= J; ++m2)
        {
          std::complex<double> cd = SEMBLE::toScalar(
              Hadron::Wigner_D(2*J,2*m1,2*m2,eul.alpha,eul.beta,eul.gamma));
          (*W)[J-m1][J-m2] = round_to_zero(cd,1e-6); 
        }

      // allow for us to treat this as tensor multiplication 
      W->lower_index(1); 

      return W; 
    }

  WignerMatrix_t*
    DMatrixManager::wigner_matrix(const rCompEulAngles &eul, 
        const int J) const
    {
      int bound = 2*J+1; 
      WignerMatrix_t *W  = new WignerMatrix_t( (TensorShape<2>())[bound][bound] , std::complex<double>(0.,0.));   
      WignerMatrix_t *tmp1; 
      WignerMatrix_t tmp2; 

      // identity 
      for(int i =0; i < bound; ++i)
        (*W)[i][i] = std::complex<double>(1.,0.); 

      // so we can abuse the tensor contraction code
      W->lower_index(1); 

      rCompEulAngles::const_iterator e; 
      for( e = eul.begin(); e != eul.end(); ++e)
      {
        // grab the pointer
        tmp1 = wigner_matrix( *e, J ); 

        // make an actual copy 
        tmp2 = *W ** tmp1 ; 

        // kill the new 
        delete tmp1; 

        // copy into W 
        *W = tmp2; 
      }

      return W; 
    }


  // notation follows notes
  WignerMatrix_t* 
    DMatrixManager::left_wigner_matrix(const rCompEulAngles &eul,
        const mom_t &l,
        const mom_t &r, 
        const int J,
        const bool use_map,
        const int print) const
    {
      int bound = 2*J+1; 
      WignerMatrix_t *W  = new WignerMatrix_t( (TensorShape<2>())[bound][bound] , std::complex<double>(0.,0.));   

      // use the precomputed wigner matricies
      if ( use_wigner_map && use_map)
      {
        wig_pair_key key(l,r,J); 
        WignerMatrixMap_t::iterator it = call_left_map()->find(key); 
        if( it == call_left_map()->end() )
        {
          std::cout << __PRETTY_FUNCTION__ << ": error missing key " 
            << key << std::endl;
          exit(1); 
        }

        *W = it->second; 

      }
      else
      {
        mom_pair_key ky = get_frame(l,r); 
        std::pair<mom_t,mom_t> can;
#ifdef BREAK_ROTATIONS_REST
        can = ky.frame_orientation(); 
#else
        can = ky.moms(); 
#endif 

        std::pair<mom_t,mom_t> moms = bb_rot_rest(l,r); 
        check_throw_frame_err(eul,moms,can);

        WignerMatrix_t *Wt,*Wn,*Wi; 

        Wt = wigner_matrix(eul,J); 
        Wn = radmat::WignerDMatrixEnv::call_factory(moms.first,J);
        Wi = radmat::WignerDMatrixEnv::call_factory(can.first,J);

        dagger(Wi); 
        dagger(Wt); 

        for(int i = 0; i < bound; ++i)
          for(int j = 0; j < bound; ++j)
            for(int k = 0; k < bound; ++k)
              for(int l = 0; l < bound; ++l)
                (*W)[i][l] += (*Wi)[i][j] * (*Wt)[j][k] * (*Wn)[k][l];

        clean(W); 

        if(print == 1)
        {
          clean(Wi); 
          clean(Wt); 
          clean(Wn); 
          //std::cout <<  "(*W)[i][l] += (*Wi)[i][j] * (*Wt)[j][k] * (*Wn)[k][l] "
          //  << "\nWi: " << *Wi << "\nWt:" << *Wt << "\nWn:" << *Wn << std::endl;
          std::cout << __func__ << " W:" << *W << std::endl;
        }

        delete Wt;
        delete Wn;
        delete Wi; 
      }

      return W; 
    }

  // notation follows notes
  WignerMatrix_t*
    DMatrixManager::right_wigner_matrix(const rCompEulAngles &eul, 
        const mom_t &l, 
        const mom_t &r, 
        const int J,
        const bool use_map,
        const int print) const
    {
      int bound = 2*J+1; 
      WignerMatrix_t *W  = new WignerMatrix_t( (TensorShape<2>())[bound][bound] , std::complex<double>(0.,0.));   

      // use the precomputed wigner matricies
      if( use_wigner_map && use_map )
      {
        wig_pair_key key(l,r,J); 
        WignerMatrixMap_t::iterator it = call_right_map()->find(key); 

        if( it == call_left_map()->end() )
        {
          std::cout << __PRETTY_FUNCTION__ << ": error missing key " 
            << key << std::endl;
          exit(1); 
        }

        *W =  it->second; 
      }
      else
      {
        mom_pair_key ky = get_frame(l,r); 
        std::pair<mom_t,mom_t> can;
#ifdef BREAK_ROTATIONS_REST
        can = ky.frame_orientation(); 
#else
        can = ky.moms(); 
#endif 

        std::pair<mom_t,mom_t> moms = bb_rot_rest(l,r); 
        check_throw_frame_err(eul,moms,can);

        WignerMatrix_t *Wt,*Wk,*Wl; 

        Wt = wigner_matrix(eul,J); 
        Wl = radmat::WignerDMatrixEnv::call_factory(moms.second,J);
        Wk = radmat::WignerDMatrixEnv::call_factory(can.second,J);

        dagger(Wl); 

        for(int i = 0; i < bound; ++i)
          for(int j = 0; j < bound; ++j)
            for(int k = 0; k < bound; ++k)
              for(int l = 0; l < bound; ++l)
                (*W)[i][l] += (*Wl)[i][j] * (*Wt)[j][k] * (*Wk)[k][l];

        clean(W); 

        if(print == 1)
        {
          clean(Wl); 
          clean(Wt); 
          clean(Wk); 
          // std::cout <<  "(*W)[i][l] += (*Wl)[i][j] * (*Wt)[j][k] * (*Wk)[k][l] "
          //   << "\nWl: " << *Wl << "\nWt:" << *Wt << "\nWk:" << *Wk << std::endl;
          std::cout << __func__ <<  " W:" << *W << std::endl;
        }

        delete Wt;
        delete Wl;
        delete Wk; 

      }

      return W; 
    }


  namespace WignerThreadMapEnv
  {
    bool registerAll(){ return init_maps(4); }
  }



}// radmat

