/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : lorentzff_canonical_rotations_checker.cc

 * Purpose :

 * Creation Date : 24-12-2013

 * Last Modified : Fri 03 Oct 2014 04:50:43 PM EDT

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/


#include "radmat/rotation_interface/rotation_interface.h"
#include "radmat/construct_data/construct_data.h"
#include "lorentzff_canonical_rotations_checker.h"
#include "lorentzff_canonical_frame_formfacs_rotation_manager.h"
#include "lorentzff_canonical_JJlist_ff.h"
#include "radmat/ff_interface/formfactor_kinematic_factors.h"
#include "semble/semble_semble.h"
#include <sstream>
#include <string> 
#include <complex>
#include <map>
#include <vector>
#include <exception>
#include "ensem/ensem.h"
#include "jackFitter/ensem_data.h"
#include "jackFitter/jackknife_fitter.h"
#include "jackFitter/ensem_data.h"
#include "jackFitter/three_point_fit_forms.h"
#include "jackFitter/plot.h"
#include <math.h>
#include <iostream>
#include <complex>
#include "adat/handle.h"


namespace radmat
{


  namespace 
  {
    ////////////////////////////////////////////////
    std::string 
      string_mom(const mom_t &p) 
      {
        std::stringstream ss; 
        ss << p[0] << p[1] << p[2] ;
        return ss.str(); 
      }


    int local_count = 0; 

    ////////////////////////////////////////////////
    bool 
      test_equivalence( const ENSEM::EnsemVectorComplex &A, 
          const ENSEM::EnsemVectorComplex &B)
      {
        ENSEM::EnsemVectorComplex C = A - B; 

        ENSEM::EnsemVectorReal real, imag;

        double const_real, const_imag; 

        real = ENSEM::real(C);
        imag = ENSEM::imag(C);

        std::vector<double> t;
        for(int i = 0; i < C.numElem(); ++i)
          t.push_back(double(i)); 

        int thigh = t.size() * 0.85; 
        int tlow = t.size() * 0.15; 

        EnsemData ereal(t,real),eimag(t,imag); 
        ereal.hideDataAboveX(thigh - 0.1); 
        eimag.hideDataBelowX(tlow -0.1); 

        ADAT::Handle<FitFunction> freal(new ThreePointConstant), fimag(new ThreePointConstant);  
        JackFit fit_real(ereal,freal), fit_imag(eimag,fimag); 

        fit_real.runAvgFit(); 
        fit_imag.runAvgFit(); 
        const_real = fit_real.getAvgFitParValue(0);
        const_imag = fit_imag.getAvgFitParValue(0);

        double const_real_var = fit_real.getAvgFitParError(0);
        double const_imag_var = fit_imag.getAvgFitParError(0);

        std::cout << __func__ << ": real = " << const_real << " +/- " << const_real_var << std::endl;
        std::cout << __func__ << ": imag = " << const_imag << " +/- " << const_imag_var << std::endl;

        // zero vy value
        if( fabs(const_real) < 1e-2 ) 
          if( fabs(const_imag) < 1e-2)
            return true; 

        // inverse covariance was singular
        if( isnan(const_real) )
          if( isnan(const_imag) )
            return true; 

        //   // zero w/ in err
        //   if( fabs(const_real) - 1.*sqrt(fabs(const_real_var)) < 0.)
        //     if( fabs(const_imag) - 1.*sqrt(fabs(const_imag_var)) < 0.)
        //       return true; 

        std::stringstream cnt; 
        cnt << ++local_count; 
        std::string count = cnt.str(); 

        std::cout << __func__ << ": reporting failure, C = A-B (should be zero) " << std::endl;
        std::cout << __func__ << ": real = " << const_real << " +/- " << const_real_var << std::endl;
        std::cout << __func__ << ": imag = " << const_imag << " +/- " << const_imag_var << std::endl;
        std::cout << "count " << count << std::endl;
        ENSEM::write( std::string("C_") + count + std::string(".jack") , C ); 
        ENSEM::write( std::string("A_") + count + std::string(".jack") , A ); 
        ENSEM::write( std::string("B_") + count + std::string(".jack") , B ); 


        return false; 
      }



    ////////////////////////////////////////////////

    typedef std::string KEY_t; 
    typedef ENSEM::EnsemVectorComplex DATA_t; 
    typedef std::map<KEY_t,DATA_t> MAP_t; 

    MAP_t init_map(const rHandle<LLSQLatticeMultiData> &d) 
    {
      MAP_t ret; 

      std::vector<ThreePointDataTag> tags = d->tags(); 
      for(int i = 0; i < tags.size(); ++i)
      {
        ADATIO::BinaryBufferWriter binBuff; 
        write(binBuff,tags[i]); 
        KEY_t tag_key = binBuff.str(); 
        DATA_t tag_data = d->get_row_ensem(i); 

        if(ret.find(tag_key) != ret.end())
        {
          std::cout << __func__ << ": key duplication error" << std::endl; 
          throw std::string("duplicate key error"); 
        }

        ret.insert(std::make_pair(tag_key,tag_data)); 
      }

      std::cout << __func__ << ": init map with " << ret.size() << " elements" << std::endl;

      return ret; 
    }

    bool my_is_nan(const std::complex<double> d)
    {
      return (isnan(std::real(d)) && isnan(std::imag(d))); 
    }

    bool check_relation(const ThreePointDataTag &tc, 
        const std::complex<double> &wc, 
        const ENSEM::EnsemVectorComplex &ec, 
        const ThreePointDataTag &t, 
        const std::complex<double> &w, 
        const ENSEM::EnsemVectorComplex &e)
    {
      // put the expected phase on the botom so that we divide / multiply by the conjugate 
      std::complex<double> factor = wc/w; 

      if( my_is_nan(factor) )
      {
        // both of them are zero -- dont bother comparing corrs  
        if(( fabs(std::real(wc)) < 1e-6) && ( fabs(std::imag(wc)) < 1e-6))
          return true; 
      } 

      ENSEM::Complex efactor = ENSEM::cmplx( ENSEM::Real(std::real(factor)) , ENSEM::Real(std::imag(factor)));

      // divide by the expected relative phase 
      return test_equivalence(ec,efactor*e); 
    } 



    void check_relations(const MAP_t &mappy)
    {
      FFKinematicFactors_t KGen;

      MAP_t::const_iterator it; 
      for(it = mappy.begin(); it != mappy.end(); ++it)
      {
        ADATIO::BinaryBufferReader BinBuffReader(it->first); 
        ThreePointDataTag t;
        read(BinBuffReader,t); 

        std::cout << __func__ << ": attempting " 
          << t.file_id << std::endl;

        // set up the canonical key  
        ThreePointDataTag canonical =t; 

        std::pair<mom_t,mom_t> canonical_mom_pair; 
        canonical_mom_pair = LatticeRotationEnv::rotation_group_can_mom(canonical.left_mom,canonical.right_mom); 
        mom_t qc = canonical_mom_pair.first - canonical_mom_pair.second; 
        canonical.left_mom = canonical_mom_pair.first; 
        canonical.q = qc; 
        canonical.right_mom = canonical_mom_pair.second; 
        canonical.file_id = generate_file_id(canonical); 

        // make a key and look for it 
        ADATIO::BinaryBufferWriter bar; 
        write(bar,canonical); 
        std::string canonical_key = bar.str(); 
        MAP_t::const_iterator can_it; 
        can_it = mappy.find(canonical_key); 

        if(can_it == mappy.end())
        {
          std::cout << __func__ <<": error canonical key missing for " 
            << canonical.rot_qsq_tag(true) << std::endl;

          std::cout << "missing " << canonical.file_id << "\n  avail:\n" << std::endl;; 

          // dump the map -- reuse the canonical iterator since 
          //    we will throw right after the map dump 
          for(can_it = mappy.begin(); can_it != mappy.end(); ++can_it)
          {
            ADATIO::BinaryBufferReader bark(can_it->first); 
            ThreePointDataTag bite;
            read(bark,bite); 
            std::cout << bite.file_id << std::endl;
          }

          // break the loop and exit code  
          throw std::string("canonical key missing"); 
        }

        ENSEM::EnsemVectorComplex data,can_data; 
        data = it->second; 
        can_data = can_it->second; 

        rHandle<Rep_p> lrep,rrep; 
        lrep = canonical.data_rep.lefty(); 
        rrep = canonical.data_rep.righty(); 

        // hardwire, recompile if you want to roll something else
        canonical.mat_elem_id = "J1mJ1m_test"; 

        // this convention must match the subduction convention 
        if(lrep->rep_type() == Stringify<CubicRep_t>()); 
        canonical.mat_elem_id += "__" + lrep->rep_id() + "," + rrep->rep_id(); 

        t.mat_elem_id = canonical.mat_elem_id; 


        FFKinematicFactors_t::KinematicFactorRow ffc,fft; 
        ffc = KGen.genFactors(&canonical); 
        fft = KGen.genFactors(&t); 

        if( !!! (check_relation( canonical, ffc.mean()[0] , can_data, t, fft.mean()[0] , data ) ) )
        {
          // bad stuff 
        }

      }
    }


    void handle_work(const rHandle<LLSQLatticeMultiData> &d)
    {
      MAP_t mappy = init_map(d); 
      check_relations(mappy); 
    }

  }



  ////////////////////////////////////////////////
  ////////////////////////////////////////////////
  ////////////////////////////////////////////////
  ////////////////////////////////////////////////
  ////////////////////////////////////////////////
  ////////////////////////////////////////////////
  ////////////////////////////////////////////////
  ////////////////////////////////////////////////
  void
    LatticeRotationRelationChecker::check(
        const rHandle<LLSQLatticeMultiData> &d) const
    {
      handle_work(d); 
    }



}

