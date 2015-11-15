/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : redstar_particle_handler_utils.cc

 * Purpose :

 * Creation Date : 20-03-2014

 * Last Modified : Thu 01 May 2014 02:35:34 PM EDT

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/

#include "redstar_particle_handler_utils.h"
#include "redstar_invert_subduction.h"
#include "redstar_cartesian_interface.h"
#include "radmat/data_representation/data_representation.h"
#include "radmat/utils/printer.h"
#include "radmat/utils/tokenize.h"
#include "hadron/irrep_util.h"
#include "ensem/ensem.h"
#include "formfac/formfac_qsq.h"
#include "hadron/ensem_filenames.h"
#include "redstar_photon_props.h"
#include "semble/semble_meta.h"
#include "itpp/itbase.h"
#include <math.h>

namespace radmat
{
  // utility functions
  namespace
  {

    struct particle_name_printer
    {
      static void print(const std::string &msg)
      {}
      // { std::cout << "particle_name_printer " << msg << std::endl;}
    };

    struct generate_block_print_name
    {
      static void print(const std::string &msg)
      {}
      // { std::cout << "generate_block_print_name " << msg << std::endl;}
    };

    struct move_lorentz4_printer
    {
      static void print(const std::string &msg) 
      {}
      // { std::cout << "move_to_lorentz4_printer" << msg << std::endl; }
    };

    std::string toString(const int i)
    {std::stringstream ss; ss << i; return ss.str();}

    std::string stringy_mom(const ADATXML::Array<int> mom)
    {
      ADATXML::Array<int> can = FF::canonicalOrder(mom);
      std::stringstream ss;
      ss << "p" << can[0] << can[1] << can[2];
      return ss.str(); 
    }

    /////////////////////////
    /////////////////////////

    BlockData make_block_data(const std::string &origin, 
        const std::string &data, 
        const int row, 
        const EnsemRedstarBlock &d)
    {
      BlockData ret; 
      ret.origin_rep = origin;
      ret.data_rep = data; 
      ret.row = row; 
      ret.data = d;
      return ret; 
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
      make_ensem_block( const std::string &name, 
          const ADATXML::Array<int> &mom,
          const int twoI_z,
          const bool creation_op,
          const bool smearedP,
          const int t_slice,
          const int J, 
          const int hel, 
          const bool parity,
          const bool isProjected )
      {
        EnsemRedstarBlock ret; 

        // play the subduction dance
        ContinuumBosonExprPrimitive meson(J,parity,hel,Hadron::generateLittleGroup(mom)); 

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
          op_name << name << "_";
          if(isProjected)
            op_name << stringy_mom(mom) << "_";
          op_name << it->m_obj.irrep; 

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

    ADATXML::Array< ADATXML::Array<int> > 
      gen_acceptable_moms(const int pmin, const int pmax)
      {
        // we dont know the size ahead of time 
        int pbound = int( ceil( sqrt(pmax) ) ); 
        std::vector<ADATXML::Array<int> > constuct;
        constuct.reserve( pmax*pmax*pmax ); 
        for(int x = -pbound; x <= pbound; ++x)
          for(int y = -pbound; y <= pbound ; ++y)
            for(int z = -pbound; z <= pbound; ++z)
            {
              int modp = x*x + y*y +z*z; 
              if( modp >= pmin )
                if( modp <= pmax )
                {
                  ADATXML::Array<int> mom(3);
                  mom[0] = x;
                  mom[1] = y;
                  mom[2] = z;
                  constuct.push_back(mom); 
                }
            }

        int sz = constuct.size(); 
        ADATXML::Array< ADATXML::Array<int> > ret(sz); 
        for(int i = 0; i < sz; ++i)
          ret[i] = constuct[i];

        return ret; 
      }

    /////////////////////////
    /////////////////////////

    struct lor_2_cub_printer
    {
      static void print(const std::string &msg)
      {}
      // {std::cout << "lor_2_cub_printer " << msg << std::endl;}
    };

    std::vector<BlockData>
      conv_lorentz_to_cubic_block(const std::vector<BlockData> &d)
      {
        std::vector<BlockData> ret; 

        typedef std::pair<Hadron::KeyHadronNPartNPtCorr_t::NPoint_t,std::string> data_t; 

        std::map<std::string,data_t> seen;  
        std::vector<BlockData>::const_iterator big_it; 
        EnsemRedstarBlock::const_iterator little_it; 

        // need to be careful about the origin in case we pass two types 
        // when building xml files!!
        for(big_it = d.begin(); big_it != d.end(); ++big_it)
        {
          std::string origin = big_it->origin_rep; 
          for(little_it = big_it->data.begin(); little_it != big_it->data.end(); ++little_it)
          {
            seen[ Hadron::ensemFileName(little_it->m_obj.irrep) + origin ] 
              = std::make_pair(little_it->m_obj,origin); 
          }
        }

        ENSEM::Complex one( cmplx(ENSEM::Real(1.),ENSEM::Real(0.))); 
        std::map<std::string,data_t>::const_iterator it; 
        for(it = seen.begin(); it != seen.end(); ++it)
        {
          std::string name = it->second.first.irrep.op.ops[1].name; 
          int row = it->second.first.irrep.row; 
          std::vector<std::string> tokens = tokenize(name, "_" ); 
          std::string rep = * ( tokens.end() -1 ); 
          EnsemRedstarBlock data; 
          data = convert_to_list( EnsemRedstarBlock::ListObj_t(one,it->second.first) );
          ret.push_back( make_block_data(it->second.second,rep,row,data) ); 
        }

        return ret; 
      }

    /////////////////////////
    /////////////////////////

    // piggy back off the single particle stuff then apply the weight 
    std::vector<BlockData> 
      generate_photon_blocks( RedstarSingleParticleMesonXML * ptr, 
          const RedstarVectorCurrentXML::insertion &ins) 
      {
        std::vector<BlockData> ret; 

        int sz = ins.photons.size(); 
        for(int i = 0; i < sz; ++i)
        {
          RedstarVectorCurrentXML::pfrag phot = ins.photons[i]; 
          ENSEM::Complex weight( cmplx( ENSEM::Real(phot.coeff_r) , ENSEM::Real(phot.coeff_i) ) ); 
          ptr->name = phot.name; 
          std::vector<BlockData> d = generate_lorentz_block( ptr ); 
          std::vector<BlockData>::iterator update;
          for(update = d.begin(); update != d.end(); ++update)
            update->data = weight * update->data;

          ret.reserve( ret.size() + d.size() );
          ret.insert(ret.end(), d.begin(), d.end() ); 
        }

        return ret; 
      }

    /////////////////////////
    /////////////////////////

    std::pair<std::string,ADATXML::Array<int> > 
      tag_block(const BlockData &d)
      {
        Hadron::KeyHadronNPartIrrep_t irrep; 
        irrep = d.data.begin()->m_obj.irrep; 
        std::stringstream ss; 
        ss << irrep.mom[0] << irrep.mom[1] << irrep.mom[2];
        ss << "h" << d.row; 
        return std::make_pair(ss.str() , irrep.mom); 
      }

    // take the O(3) bit of O(3) + 1 data and move from
    // the 3 dimensional polarization to the vector index
    // of the vector curent
    //    NB: only use it for unimproved data 
    std::vector<BlockData>
      move_to_lorentz4(const std::vector<BlockData> &helicity3)
      {
        std::vector<BlockData> ret; 

        printer_function<move_lorentz4_printer>( 
            "helicity.size " + toString( int(helicity3.size())) ); 


        // all the data merged together 
        std::map<std::string,BlockData> merge_blocks; 

        // the unique momenta
        std::map<std::string, ADATXML::Array<int> > moms; 

        std::vector<BlockData>::const_iterator block_it; 
        std::map<std::string,BlockData>::iterator b; 

        // group them on helicities 
        for(block_it = helicity3.begin(); block_it != helicity3.end(); ++block_it)
        {
          std::pair<std::string,ADATXML::Array<int> > foo = tag_block(*block_it); 
          b = merge_blocks.find(foo.first); 

          // first time we see it add a momentum 
          if( b == merge_blocks.end() )
          {
            merge_blocks.insert(std::make_pair(foo.first,*block_it));
            std::stringstream ss;
            ss << foo.second[0] << foo.second[1] << foo.second[2]; 
            moms.insert(std::make_pair(ss.str(),foo.second));  
          }
          else
          {
            b->second.data = b->second.data + block_it->data; 
          }
        }

        std::map<std::string,ADATXML::Array<int> >::const_iterator it; 
        for(it = moms.begin(); it != moms.end(); ++it)
        {
          itpp::Vec<EnsemRedstarBlock> cart,hel(3); 
          itpp::Mat<std::complex<double> > M = invert2Cart(it->second,PHOTON_CREATE); 
          for(int h = -1; h < 2; ++h)
          {
            std::stringstream ss; ss << h; 
            b = merge_blocks.find( it->first + "h" + ss.str() ); 
            if( b == merge_blocks.end() )
            {
              std::cout << __PRETTY_FUNCTION__ 
                << " error, missing " << it->first << std::endl;
              exit(1); 
            }
            hel[1-h] = b->second.data; 
          }

          b = merge_blocks.find( it->first + "h0" ); 

          cart = M * hel; 

          ret.push_back( make_block_data( b->second.origin_rep , b->second.data_rep, 1, cart[0] )); 
          ret.push_back( make_block_data( b->second.origin_rep , b->second.data_rep, 2, cart[1] )); 
          ret.push_back( make_block_data( b->second.origin_rep , b->second.data_rep, 3, cart[2] )); 
        }

        return ret; 
      }

    /////////////////////////
    /////////////////////////

    struct resum_printer_2
    {
      static void print(const std::string &msg)
      {}
      // {std::cout << "resum_printer_2 " << msg << std::endl; }
    }; 

    BlockData
      resum_ensem_redstar_block(const BlockData &d)
      {
        BlockData ret; 
        ret = d; 

        printer_function<resum_printer_2>( " origin: " + d.origin_rep); 

        EnsemRedstarBlock::const_iterator it; 
        std::map<std::string, 
          std::pair<EnsemRedstarBlock::Coeff_t , EnsemRedstarBlock::Obj_t> > coeff_map; 

        for( it = d.data.begin(); it != d.data.end(); ++it)
        {
          ENSEM::Complex coeff = it->m_coeff; 
          printer_function<particle_name_printer>( 
              it->m_obj.irrep.op.ops[1].name );

          std::string key = Hadron::ensemFileName(it->m_obj.irrep); 

          if( coeff_map.find(key) != coeff_map.end() )
            coeff += coeff_map[key].first; 

          // reinsert 
          coeff_map[key] = std::make_pair(coeff,it->m_obj); 
        } 

        EnsemRedstarBlock sum; 
        std::map<std::string,
          std::pair<EnsemRedstarBlock::Coeff_t , EnsemRedstarBlock::Obj_t> 
            >::const_iterator map_it; 

        std::complex<double> zero(0.,0.); 
        for(map_it = coeff_map.begin(); map_it != coeff_map.end(); ++map_it)
          if( SEMBLE::toScalar(map_it->second.first) != zero) 
            sum = sum + EnsemRedstarBlock::ListObj_t(map_it->second.first,map_it->second.second); 

        // update the block
        ret.data = sum;

        return ret; 
      }

    /////////////////////////
    /////////////////////////

    std::vector<BlockData> 
      resum_ensem_redstar_blocks(const std::vector<BlockData> &d)
      {
        std::vector<BlockData> ret;
        ret.reserve(d.size());
        std::vector<BlockData>::const_iterator it; 
        for(it = d.begin(); it != d.end(); ++it)
          ret.push_back(resum_ensem_redstar_block(*it));
        return ret; 
      }


  } // anonomyous 


  std::vector<BlockData>
    generate_lorentz_block( const RedstarSingleParticleMesonXML * const ptr)
    {
      std::vector<BlockData> ret; 

      rHandle<LorentzRep> cont_rep;
      cont_rep = LorentzRepresentationFactoryEnv::callFactory(ptr->cont_rep); 
      int spin = cont_rep->rep_spin(); 
      bool parity = cont_rep->rep_parity() == 1 ? true : false;

      int psz = ptr->mom.size();
      int hsz = ptr->H.size();
      for(int p = 0; p < psz; ++p)
        for(int h = 0; h < hsz; ++h)
        {
          EnsemRedstarBlock block; 

          printer_function<generate_block_print_name>( ptr->name ); 

          block = make_ensem_block( ptr->name, 
              ptr->mom[p],
              ptr->twoI_z,
              ptr->creation_op,
              ptr->smearedP,
              ptr->t_slice,
              spin,
              ptr->H[h],
              parity,
              ptr->isProjected); 

          ret.push_back( make_block_data( cont_rep->rep_id(), 
                cont_rep->rep_id(),
                ptr->H[h],
                block)); 
        }

      return resum_ensem_redstar_blocks(ret); 
    }


  std::vector<BlockData>
    generate_cubic_block( const RedstarSingleParticleMesonXML * const ptr)
    {
      std::vector<BlockData> cont = generate_lorentz_block(ptr);
      return conv_lorentz_to_cubic_block(cont); 
    }

  struct vc_lor_printer
  {
    static void print(const std::string &msg)
    {}
    // {std::cout << "vc_lor_printer " << msg << std::endl;}
  };

  std::vector<BlockData>
    generate_lorentz_block( const RedstarVectorCurrentXML * const ptr)
    {
      std::vector<BlockData> ret; 
      ADATXML::Array< ADATXML::Array<int> > moms;
      moms = gen_acceptable_moms(ptr->pmin,ptr->pmax); 

      RedstarSingleParticleMesonXML * gamma = new RedstarSingleParticleMesonXML(); 

      gamma->mom = moms; 
      gamma->twoI_z = 0; 
      gamma->creation_op = PHOTON_CREATE; 
      gamma->isProjected = false; 
      gamma->t_slice = ptr->t_slice; 

      if( ptr->time.active )
      {
        ADATXML::Array<int> H(1);
        H[0] = 0; 
        gamma->H = H; 
        gamma->cont_rep = "J0p"; 
        gamma->smearedP = ptr->time.smearedP; 
        std::vector<BlockData> foo = generate_photon_blocks( gamma , ptr->time); 
        ret.reserve( ret.size() + foo.size() ); 
        ret.insert( ret.end() , foo.begin() , foo.end() ); 
      }
      if( ptr->space.active )
      {
        ADATXML::Array<int> H(3);
        H[0] = 1; 
        H[1] = 0;
        H[2] = -1; 
        gamma->H = H; 
        gamma->cont_rep = "J1m"; 
        gamma->smearedP = ptr->space.smearedP; 
        std::vector<BlockData> foo = generate_photon_blocks( gamma , ptr->space); 
        std::vector<BlockData> bar = move_to_lorentz4( foo ); 

        ret.reserve( ret.size() + bar.size() ); 
        ret.insert( ret.end() , bar.begin() , bar.end() ); 
      }

      std::vector<BlockData>::const_iterator it; 
      for(it = ret.begin(); it != ret.end(); ++it)
        printer_function<vc_lor_printer>(it->origin_rep); 

      return resum_ensem_redstar_blocks(ret); 
    }


  struct vc_cub_printer
  {
    static void print(const std::string &msg)
    {}
    // {std::cout << "vc_cub_printer " << msg << std::endl;}
  };

  std::vector<BlockData>
    generate_cubic_block( const RedstarVectorCurrentXML * const ptr)
    {
      // these two steps get us an unweighted list of everything 
      std::vector<BlockData> cont = generate_lorentz_block(ptr);
      std::vector<BlockData> resum = conv_lorentz_to_cubic_block(cont); 

      // now go back through and apply the photon weights 

      // first build a map of the weights 
      std::map<std::string,ENSEM::Complex> weights; 
      if(ptr->time.active)
      {
        for(int i = 0; i < ptr->time.photons.size(); ++i)
        {
          RedstarVectorCurrentXML::pfrag frag;
          frag = ptr->time.photons[i]; 
          ENSEM::Complex weight( cmplx(ENSEM::Real(frag.coeff_r) , ENSEM::Real(frag.coeff_i) )); 
          weights.insert(std::make_pair(frag.name,weight) );
        }
      }
      if(ptr->space.active)
      {
        for(int i = 0; i < ptr->space.photons.size(); ++i)
        {
          RedstarVectorCurrentXML::pfrag frag;
          frag = ptr->space.photons[i]; 
          ENSEM::Complex weight( cmplx(ENSEM::Real(frag.coeff_r) , ENSEM::Real(frag.coeff_i)) ); 
          weights.insert(std::make_pair(frag.name,weight) );
        }
      }

      // now loop over them and apply the weights
      std::vector<BlockData>::iterator it; 
      EnsemRedstarBlock::const_iterator block_it; 
      std::map<std::string,ENSEM::Complex>::const_iterator weight_it; 

      // reweight them for the photon weight 
      for( it = resum.begin(); it != resum.end(); ++it)
      {
        printer_function<vc_cub_printer>( it->origin_rep ); 

        EnsemRedstarBlock update; 
        for( block_it = it->data.begin(); block_it != it->data.end(); ++block_it)
        {
          std::string name = block_it->m_obj.irrep.op.ops[1].name; 
          ENSEM::Complex weight; 
          bool found = false; 

          for( weight_it = weights.begin(); weight_it != weights.end(); ++weight_it)
          {
            // is there a partial match ( doesnt care about the subduce part )
            if(name.find(weight_it->first) != std::string::npos)
            {
              weight = weight_it->second; 
              found = true; 
              break; 
            }
          }

          // we should never ever see this 
          if( !!! found )
          {
            printer_function<console_print>("missing " + name ); 
            throw std::string("photon weight error"); 
          } 

          // add it to the block update 
          update = update + weight * *block_it; 
        }

        // update the block to have the correct weights 
        it->data = update; 
      }

      return resum; 
    }


  std::vector<VectorCurrentImprovedBlockData>
    generate_lorentz_block(const RedstarImprovedVectorCurrentXML *const ptr)
    {
      std::vector<VectorCurrentImprovedBlockData> ret; 
      std::vector<BlockData>::const_iterator it; 
      std::vector<BlockData> piggy; 
      RedstarVectorCurrentXML back; 
      back.pmin = ptr->pmin; 
      back.pmax = ptr->pmax; 
      back.t_slice = ptr->t_slice; 
      back.time = ptr->time; 
      back.space = ptr->space;

      piggy = generate_lorentz_block( &back ); 
      ret.reserve(piggy.size()); 

      for(it = piggy.begin(); it != piggy.end(); ++it)
      {
        VectorCurrentImprovedBlockData foo; 
        foo.origin_rep = it->origin_rep; 
        foo.data_rep = it->data_rep; 
        foo.row = it->row; 
        foo.data = it->data; 
        foo.imp = ptr->imp; 
        ret.push_back(foo); 
      }

      return ret; 
    }


  struct improved_vc_printer
  {
    static void print(const std::string &msg)
    {}
    // {std::cout << "improved_vc_printer " << msg << std::endl;}
  };

  std::vector<VectorCurrentImprovedBlockData>
    generate_cubic_block(const RedstarImprovedVectorCurrentXML *const ptr)
    {
      std::vector<VectorCurrentImprovedBlockData> ret; 
      std::vector<BlockData>::const_iterator it; 
      std::vector<BlockData> piggy; 
      RedstarVectorCurrentXML back; 
      back.pmin = ptr->pmin; 
      back.pmax = ptr->pmax; 
      back.t_slice = ptr->t_slice; 
      back.time = ptr->time; 
      back.space = ptr->space;

      piggy = generate_cubic_block( &back ); 
      ret.reserve(piggy.size()); 

      for(it = piggy.begin(); it != piggy.end(); ++it)
      {
        VectorCurrentImprovedBlockData foo; 
        foo.origin_rep = it->origin_rep; 
        foo.data_rep = it->data_rep; 
        foo.row = it->row; 
        foo.data = it->data; 
        foo.imp = ptr->imp; 
        ret.push_back(foo); 

        printer_function<improved_vc_printer>( " origin " + foo.origin_rep ); 
      }

      return ret; 
    }
} 
