/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : redstar_three_point_utils.cc

 * Purpose :

 * Creation Date : 21-03-2014

 * Last Modified : Wed 13 Aug 2014 08:31:15 PM EDT

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/


#include "redstar_three_point_utils.h"
#include "redstar_improved_vector_current_xml_interface.h"
#include "redstar_cartesian_interface.h"
#include "adat/map_obj.h"
#include "semble/semble_meta.h"
#include "formfac/formfac_qsq.h"
#include "redstar_invert_subduction.h"
#include "hadron/irrep_util.h"
#include "hadron/ensem_filenames.h"
#include "redstar_photon_props.h"
#include "radmat/utils/tokenize.h"
#include "radmat/utils/printer.h"


#define ALLOW_CACHE_STUB


namespace radmat
{
  namespace
  {
    // returns true if momentum is conserved
    bool check_mom(const BlockData &l, const BlockData &g, const BlockData &r)
    {
      ADATXML::Array<int> ll,gg,rr;
      ll = l.data.begin()->m_obj.irrep.mom;
      gg = g.data.begin()->m_obj.irrep.mom;
      rr = r.data.begin()->m_obj.irrep.mom;

      int ls(1),gs(1),rs(1); 

      if( !!! l.data.begin()->m_obj.irrep.creation_op )
        ls = - 1 ; 

      if( !!! g.data.begin()->m_obj.irrep.creation_op )
        gs = - 1 ; 

      if( !!! r.data.begin()->m_obj.irrep.creation_op )
        rs = - 1 ; 

      for(int i =0; i < 3; ++i)
        if( ls*ll[i] + gs*gg[i] + rs*rr[i] != 0 )
          return false;  

      return true; 
    }

    // returns true if momentum is conserved
    bool check_mom(const BlockData &l, const VectorCurrentImprovedBlockData &g, const BlockData &r)
    {
      ADATXML::Array<int> ll,gg,rr;
      ll = l.data.begin()->m_obj.irrep.mom;
      gg = g.data.begin()->m_obj.irrep.mom;
      rr = r.data.begin()->m_obj.irrep.mom;

      int ls(1),gs(1),rs(1); 

      if( !!! l.data.begin()->m_obj.irrep.creation_op )
        ls = - 1 ; 

      if( !!! g.data.begin()->m_obj.irrep.creation_op )
        gs = - 1 ; 

      if( !!! r.data.begin()->m_obj.irrep.creation_op )
        rs = - 1 ; 

      for(int i =0; i < 3; ++i)
        if( ls*ll[i] + gs*gg[i] + rs*rr[i] != 0 )
          return false;  

      return true; 
    }


    EnsemRedstarNPtBlock 
      merge_ensem_blocks( const EnsemRedstarBlock &lefty,
          const EnsemRedstarBlock &gamma,
          const EnsemRedstarBlock &righty, 
          const std::string &ensemble)
      { 
        EnsemRedstarNPtBlock ret; 
        EnsemRedstarBlock::const_iterator l,g,r;

        for( l = lefty.begin(); l != lefty.end(); ++l)
          for( g = gamma.begin(); g != gamma.end(); ++g)
            for( r = righty.begin(); r != righty.end(); ++r)  
            {
              ENSEM::Complex coeff;
              Hadron::KeyHadronNPartNPtCorr_t npt; 
              coeff = l->m_coeff * g->m_coeff * r->m_coeff; 
              npt.npoint.resize(3); 

              // arrays are FORTRAN style 
              npt.npoint[1] = l->m_obj;
              npt.npoint[2] = g->m_obj;
              npt.npoint[3] = r->m_obj;
              npt.ensemble = ensemble; 

              ret = ret + EnsemRedstarNPtBlock::ListObj_t(coeff,npt); 
            }

        return ret; 
      } 


    ThreePointData
      make_data(const BlockData &l, 
          const BlockData &g, 
          const BlockData &r, 
          const std::string &ensemble)
      {
        ThreePointData ret; 
        ret.origin_rep = DataRep3pt(l.origin_rep,
            g.origin_rep,
            r.origin_rep);
        ret.data_rep = DataRep3pt(l.data_rep,
            g.data_rep,
            r.data_rep);
        ret.left_row = l.row;
        ret.gamma_row = g.row;
        ret.right_row = r.row; 

        ret.data = merge_ensem_blocks( l.data,g.data,r.data,ensemble);

        return ret; 
      }

    struct resum_printer
    {
      static void print(const std::string &msg)
      {}
      // {  std::cout << "resum_printer " << msg << std::endl; } 
    };

    EnsemRedstarNPtBlock 
      resum_ensem_npt_block(const EnsemRedstarNPtBlock &b)
      {
        EnsemRedstarNPtBlock ret; 

        typedef std::pair<EnsemRedstarNPtBlock::Coeff_t, 
                EnsemRedstarNPtBlock::Obj_t> coeff_object; 
        std::map<std::string,coeff_object> coeff_map; 
        EnsemRedstarNPtBlock::const_iterator it; 

        for(it = b.begin(); it != b.end(); ++it)
        {
          EnsemRedstarNPtBlock::Coeff_t coeff = it->m_coeff; 
          std::string key = Hadron::ensemFileName(it->m_obj); 

          if(coeff_map.find(key) != coeff_map.end())
            coeff += coeff_map[key].first; 
          
          // reinsert 
          coeff_map[key] = std::make_pair(coeff,it->m_obj); 
        }

        std::map<std::string,coeff_object>::const_iterator map_it; 
        for(map_it = coeff_map.begin(); map_it != coeff_map.end(); ++map_it)
          if( fabs(SEMBLE::toScalar(map_it->second.first) ) > 1e-6)
          {
            ret = ret + EnsemRedstarNPtBlock::ListObj_t(map_it->second.first, 
                map_it->second.second); 
            printer_function<resum_printer>(map_it->first); 
          }

        return ret; 
      }


    ThreePointData
      resum_three_point_data(const ThreePointData &d)
      {
        ThreePointData ret; 
        ret.origin_rep = d.origin_rep; 
        ret.data_rep = d.data_rep; 
        ret.left_row = d.left_row; 
        ret.right_row = d.right_row; 
        ret.gamma_row = d.gamma_row; 
        ret.data = resum_ensem_npt_block(d.data); 

        return ret; 
      }

    std::vector<ThreePointData>
      resum_three_point_data(const std::vector<ThreePointData> &d)
      {
        std::vector<ThreePointData> ret; 
        std::vector<ThreePointData>::const_iterator it; 
        ret.reserve(d.size()); 

        for(it = d.begin(); it != d.end(); ++it)
        {
          ThreePointData tmp = resum_three_point_data(*it); 
          if( tmp.data.size() == 0 )
            continue; 
          ret.push_back(tmp); 
        }

        return ret; 
      }


  } // anonomyous 


  //
  //
  // MERGE THE EASY UNINPROVED STUFF 
  //
  //
  std::vector<ThreePointData>
    merge_blocks(const std::vector<BlockData> &lefty, 
        const std::vector<BlockData> &gamma,
        const std::vector<BlockData> &righty,
        const std::string &ensemble)
    {
      std::vector<ThreePointData> ret;
      std::vector<BlockData>::const_iterator l,r,g;

      // some poor estimate of the size, lots of them are removed via 
      // momentum conservation
      ret.reserve( lefty.size() * righty.size() / 2 ); 

      for( l = lefty.begin(); l != lefty.end(); ++l)
        for( g = gamma.begin(); g != gamma.end(); ++g)
          for(r = righty.begin(); r != righty.end(); ++r)
          {
            if( check_mom(*l,*g,*r) )
              ret.push_back( make_data(*l,*g,*r,ensemble) ); 
          }

      return resum_three_point_data( ret ); 
    }


  // the harder problem 
  namespace 
  {

    typedef ADAT::MapObject<Hadron::KeyHadronNPartIrrep_t,double> singleThreadMassCache; 
    typedef std::pair<ThreePointData,RedstarImprovedVectorCurrentXML::improvement>
      ImprovedThreePointData; 


    ImprovedThreePointData 
      make_data(const BlockData &l, 
          const VectorCurrentImprovedBlockData &g, 
          const BlockData &r, 
          const std::string &ensemble)
      {
        ThreePointData ret; 
        ret.origin_rep = DataRep3pt(l.origin_rep,
            g.origin_rep,
            r.origin_rep);
        ret.data_rep = DataRep3pt(l.data_rep,
            g.data_rep,
            r.data_rep);
        ret.left_row = l.row;
        ret.gamma_row = g.row;
        ret.right_row = r.row; 

        ret.data = merge_ensem_blocks( l.data,g.data,r.data,ensemble);

        return std::make_pair(ret,g.imp); 
      }

    // the database type we will be using 
    typedef radmatAllConfDatabaseInterface< Hadron::KeyHadronNPartNPtCorr_t,
            ENSEM::EnsemVectorComplex,
            RadmatExtendedKeyHadronNPartIrrep_t,
            RadmatMassOverlapData_t> DBType_t;

    // key type
    typedef RadmatExtendedKeyHadronNPartIrrep_t DBKeyType_t; 

    // get all of the keys from a single data -- may be more than one pair  
    std::vector< std::pair<DBKeyType_t,DBKeyType_t> > 
      make_keys(const ThreePointData &d,const RedstarThreePointXMLInput &inp)
      {
        std::vector<std::pair<DBKeyType_t,DBKeyType_t> > ret; 

        EnsemRedstarNPtBlock::const_iterator it; 
        for(it = d.data.begin(); it != d.data.end(); ++it)
        {
          Hadron::KeyHadronNPartNPtCorr_t npt = it->m_obj; 

          ret.push_back( std::make_pair( DBKeyType_t(inp.pid_left,npt.npoint[1].irrep), 
                DBKeyType_t(inp.pid_right,npt.npoint[3].irrep)) ); 
        }

        return ret; 
      }


    singleThreadMassCache generate_mass_cache( const std::vector<ImprovedThreePointData> &d, 
        const RedstarThreePointXMLInput &inp)
    {
      singleThreadMassCache cache; 

      DBType_t db(inp.db_props); 
      std::vector<ImprovedThreePointData>::const_iterator it; 

      // loop the data and pull the relevant masses
      for(it = d.begin(); it != d.end(); ++it)
      {
        std::vector< std::pair<DBKeyType_t,DBKeyType_t> > keys = make_keys(it->first,inp); 
        std::vector< std::pair<DBKeyType_t,DBKeyType_t> >::const_iterator k;

        // loop the set of keys on this bit of data 
        for(k = keys.begin(); k != keys.end(); ++k)
        {
          DBKeyType_t left_cache_key = k->first; 
          DBKeyType_t right_cache_key = k->second; 

          DBKeyType_t left_db_key = left_cache_key; 
          DBKeyType_t right_db_key = right_cache_key; 

          // the radmat database keys are hacked to condense 
          // the number of copies, eg rows collapse to 0, 
          // and the stars of p collapse to a single ref
          left_db_key.doLG_symmetry(); 
          right_db_key.doLG_symmetry(); 

          // stupid large 
          double l = 1000; 
          double r = 2000;

          // check insert left
          if( !!! cache.exist(left_cache_key.m_basic_key) )
          {
            // update or use default
            if(db.exists(left_db_key) )
            {
              RadmatMassOverlapData_t tmp = db.fetch(left_db_key); 
              l = SEMBLE::toScalar(ENSEM::mean(tmp.E())); 
            }
            else
            {
              std::cout << "Error: missing db key " 
                  << Hadron::ensemFileName(left_db_key.m_basic_key) << std::endl;
              std::cout << "Key: \n" << left_db_key << std::endl;
#ifdef ALLOW_CACHE_STUB
              std::cout << "continuuing with a stub value of " << l << std::endl;
#else
              POW2_ASSERT(false); 
#endif
            }

            cache.insert(left_cache_key.m_basic_key,l); 
          }

          // check insert right
          if( !!! cache.exist(right_cache_key.m_basic_key) )
          {
            // update or use default 
            if(db.exists(right_db_key) )
            {
              RadmatMassOverlapData_t tmp = db.fetch(right_db_key); 
              r = SEMBLE::toScalar(ENSEM::mean(tmp.E())); 
            }
            else
            {
              std::cout << "Error: missing db key " 
                  << Hadron::ensemFileName(right_db_key.m_basic_key) << std::endl;
              std::cout << "Key: \n" << right_db_key << std::endl;
#ifdef ALLOW_CACHE_STUB
              std::cout << "continuuing with a stub value of " << r << std::endl;
#else
              POW2_ASSERT(false); 
#endif
            }

            cache.insert(right_cache_key.m_basic_key,r); 
          }

        } // k 
      } // it 

      return cache; 
    }

    /////////////////////////
    /////////////////////////

    Hadron::KeyHadronNPartNPtCorr_t::NPoint_t
      make_npt( const std::string &name, 
          const ADATXML::Array<int> &mom,
          const int row, 
          const int twoI_z,
          const bool creation_op,
          const bool smearedP,
          const int t_slice )
      {

        Hadron::KeyHadronNPartIrrep_t base; 

        //    IN WITH THAT NEW STUFF
        base.op.ops.resize(1); 
        base.op.ops[1].name = name;
        base.op.ops[1].mom_type = FF::canonicalOrder(mom); 

        base.row = row;
        base.twoI_z = twoI_z;
        base.mom = mom;
        base.creation_op = creation_op;
        base.smearedP = smearedP;

        Hadron::KeyHadronNPartNPtCorr_t::NPoint_t npt; 
        npt.t_slice = t_slice; 
        npt.irrep = base; 

        return npt;
      }

    /////////////////////////
    /////////////////////////

    EnsemRedstarBlock
      make_ensem_improvement_block( const std::string &name, 
          const ADATXML::Array<int> &mom,
          const int twoI_z,
          const bool creation_op,
          const bool smearedP,
          const int t_slice,
          const int hel, 
          const bool parity)
      {
        EnsemRedstarBlock ret; 

        // play the subduction dance
        ContinuumBosonExprPrimitive meson(1,parity,hel,Hadron::generateLittleGroup(mom)); 

        // the lattice version 
        ListLatticeIrrepExpr_t lattice_meson;

        // handle conventions
        if(creation_op)
          lattice_meson = invertSubduction(meson); 
        else
          lattice_meson = conj(invertSubduction(meson)); 

        // do this sum
        // O_{J,H} \sim \sum_{\lambda,\mu} S_{J,H}^{\lambda,\mu} O^{\lambda,\mu}
        ListLatticeIrrepExpr_t::const_iterator it; 
        for(it = lattice_meson.begin(); it != lattice_meson.end(); ++it)
        {
          std::stringstream op_name;
          op_name << name << "_" << it->m_obj.irrep; 

          Hadron::KeyHadronNPartNPtCorr_t::NPoint_t npt; 
          npt = make_npt( op_name.str(), 
              mom, 
              it->m_obj.row, 
              twoI_z,
              creation_op,
              smearedP,
              t_slice );

          ret = ret + EnsemRedstarBlock::ListObj_t(it->m_coeff,npt); 
        }

        return ret; 
      }

    /////////////////////////
    /////////////////////////
    itpp::Vec< EnsemRedstarBlock >
      make_ensem_improvement_blocks( const std::string &name, 
          const ADATXML::Array<int> &mom,
          const int twoI_z,
          const bool creation_op,
          const bool smearedP)
      {
        itpp::Vec<EnsemRedstarBlock> ret(3); 
        for(int i = 0; i < 3; ++i)
          ret[i] = make_ensem_improvement_block(name,mom,twoI_z,
              PHOTON_CREATE,smearedP,-3,1-i,false);
        return ret; 
      }

    itpp::Vec<EnsemRedstarBlock>
      make_cart_ensem_improvement_blocks(const std::string &name, 
          const Hadron::KeyHadronNPartIrrep_t &base,
          const int t_slice)
      {
        itpp::Vec<EnsemRedstarBlock> hel;
        hel = make_ensem_improvement_blocks( name, 
            base.mom, 
            base.twoI_z,
            base.creation_op,
            base.smearedP);

        itpp::Mat<std::complex<double> > M = invert2Cart(base.mom,PHOTON_CREATE); 

        return M * hel; 
      } 

    EnsemRedstarBlock 
      temporal_improvement( const std::string &name, 
          const Hadron::KeyHadronNPartNPtCorr_t::NPoint_t &base,
          const std::complex<double> pre_factor, 
          const double mom_fac)
      {
        EnsemRedstarBlock ret; 

        itpp::Vec<EnsemRedstarBlock> cart;
        cart = make_cart_ensem_improvement_blocks(name,base.irrep,base.t_slice); 
        ADATXML::Array<int> mom = base.irrep.mom; 

        for(int i = 0; i < 3; ++i)
        {
          ENSEM::Complex weight; 
          weight = SEMBLE::toScalar(pre_factor*mom_fac*double(mom[i])); 
          ret = ret + weight * cart[i]; 
        }

        return ret;  
      }

    struct temporal_printer
    {
      static void print(const std::string &msg)
      {}
      // {std::cout << "temporal_printer " << msg << std::endl;}
    };

    // this one is a bit harder, lets do it using cartesian components 
    ThreePointData apply_temporal_improvement(const ThreePointData &d, 
        const RedstarImprovedVectorCurrentXML::improvement &inp, 
        const double mom_fac)
    {
      ThreePointData ret; 
      ret.origin_rep = d.origin_rep; 
      ret.data_rep = d.data_rep; 
      ret.left_row = d.left_row; 
      ret.right_row = d.right_row; 
      ret.gamma_row = d.gamma_row; 
      ret.data = d.data; 

      double nu = inp.Nu_s; 
      double xi = inp.Xi_0; 
      double coeff_r = inp.coeff_r; 
      double coeff_i = inp.coeff_i; 
      std::string op_stem = inp.name; 
      std::complex<double> pre_factor = -0.25*nu*(1. - 1./xi)*std::complex<double>(coeff_r,coeff_i); 


      EnsemRedstarNPtBlock::const_iterator it; 
      // loop the data, find the terms then simply stick in the 
      // improvement name with an updated weight 
      for(it = d.data.begin(); it != d.data.end(); ++it)
      {
        // the key sits in the xml, use it to find the momentum 
        Hadron::KeyHadronNPartNPtCorr_t npt = it->m_obj; 
        ENSEM::Complex weight = it->m_coeff; 

        EnsemRedstarBlock improvement = temporal_improvement( op_stem, 
            npt.npoint[2], pre_factor,mom_fac); 

        EnsemRedstarBlock::const_iterator bit; 
        for(bit = improvement.begin(); bit != improvement.end(); ++bit)
        {
          ENSEM::Complex w = weight * bit->m_coeff; 
          
          if(std::norm( SEMBLE::toScalar(w) ) < 1e-6)
            continue; 

          Hadron::KeyHadronNPartNPtCorr_t corr = npt; 
          corr.npoint[2] = bit->m_obj; 
          ret.data.push_back( EnsemRedstarNPtBlock::ListObj_t(w,corr)); 
          printer_function<temporal_printer>( ensemFileName(corr) ); 
        }

      }


      return ret; 
    }


    struct spatial_printer
    {
      static void print(const std::string &msg)
      {}
      // { std::cout << "spatial_printer " << msg << std::endl;}
    };

    // do spatial improvement, the easy case 
    ThreePointData apply_spatial_improvement(const ThreePointData &d, 
        const RedstarImprovedVectorCurrentXML::improvement &inp, 
        const singleThreadMassCache &cache, 
        const double mom_fac)
    {
      ThreePointData ret; 
      ret.origin_rep = d.origin_rep; 
      ret.data_rep = d.data_rep; 
      ret.left_row = d.left_row; 
      ret.right_row = d.right_row; 
      ret.gamma_row = d.gamma_row; 
      ret.data = d.data; 

      double nu = inp.Nu_s; 
      double xi = inp.Xi_0; 
      double coeff_r = inp.coeff_r; 
      double coeff_i = inp.coeff_i; 
      std::string op_stem = inp.name; 


      std::complex<double> pre_factor = 0.25*( 1 - xi )*std::complex<double>(coeff_r,coeff_i); 

      EnsemRedstarNPtBlock::const_iterator it; 
      // loop the data, find the terms then simply stick in the 
      // improvement name with an updated weight 
      for(it = d.data.begin(); it != d.data.end(); ++it)
      {
        // the key sits in the xml, use it to find the zero component 
        // of the momentum transfer 
        Hadron::KeyHadronNPartNPtCorr_t npt = it->m_obj; 
        ENSEM::Complex weight = it->m_coeff; 
        double l = cache[npt.npoint[1].irrep];
        double r = cache[npt.npoint[3].irrep]; 
        double q0 = l - r; 

        if( !!! PHOTON_CREATE )
          q0 = -q0; 

        // update the weight to be prefactor * q_0
        weight = weight * SEMBLE::toScalar( pre_factor*q0); 

        if( std::norm( SEMBLE::toScalar(weight) ) < 1e-6)
          continue;

        // find the photon name 
        std::string phot_op_name = npt.npoint[2].irrep.op.ops[1].name; 

        // tokenize it, last one is the rep 
        std::vector<std::string> tokens = tokenize(phot_op_name,"_"); 
        std::string imp_op_name = op_stem + "_" + * ( tokens.end() -1); 

        // update npt to have the improvement as the name
        npt.npoint[2].irrep.op.ops[1].name = imp_op_name; 

        // throw it back into the data list 
        ret.data.push_back( EnsemRedstarNPtBlock::ListObj_t(weight,npt) ); 

        printer_function<spatial_printer>( Hadron::ensemFileName(npt) ); 
      }

      return ret; 
    }


    // break on time vs space  
    ThreePointData apply_improvement(const ImprovedThreePointData &d, 
        const singleThreadMassCache &cache, 
        const double mom_fac)
    {
      if( d.first.origin_rep.g == "J0p" )
        return apply_temporal_improvement(d.first,d.second,mom_fac); 
      else if (d.first.origin_rep.g == "J1m")
        return apply_spatial_improvement(d.first,d.second,cache,mom_fac); 
      else
      {
        std::cout << "unknown rep " << d.first.origin_rep.g << std::endl;
        POW2_ASSERT(false); 
      }

      return d.first; 
    }

  } // anonomyous 

  std::vector<ThreePointData> 
    merge_blocks(const std::vector<BlockData> &lefty, 
        const std::vector<VectorCurrentImprovedBlockData> &gamma, 
        const std::vector<BlockData> &righty, 
        const std::string &ensemble,
        const RedstarThreePointXMLInput &inp)
    {
      std::vector<ThreePointData> ret; 
      std::vector<ImprovedThreePointData> merge; 
      std::vector<BlockData>::const_iterator l,r;
      std::vector<VectorCurrentImprovedBlockData>::const_iterator g; 

      // some poor estimate of the size, lots of them are removed via 
      // momentum conservation
      ret.reserve( lefty.size() * righty.size() / 2 ); 

      for( l = lefty.begin(); l != lefty.end(); ++l)
        for( g = gamma.begin(); g != gamma.end(); ++g)
          for(r = righty.begin(); r != righty.end(); ++r)
          {
            if( check_mom(*l,*g,*r) )
              merge.push_back( make_data(*l,*g,*r,ensemble) ); 
          }

      ret.reserve(merge.size()); 
      singleThreadMassCache cache; 
      cache = generate_mass_cache(merge,inp); 

      std::vector<ImprovedThreePointData>::const_iterator it; 
      for(it = merge.begin(); it != merge.end(); ++it)
        ret.push_back( apply_improvement( *it, cache, inp.mom_fac ) ); 

      return resum_three_point_data( ret ); 
    }


}

