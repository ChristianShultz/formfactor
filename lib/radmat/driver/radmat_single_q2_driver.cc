/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name :

 * Purpose :

 * Creation Date : 25-02-2013

 * Last Modified : Sat 18 Oct 2014 12:08:07 PM EDT

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/



#include "radmat/driver/radmat_single_q2_driver.h"
#include "radmat/llsq/llsq_multi_data_serialize.h"
#include "radmat/llsq/llsq_formfactor_data.h"
#include "radmat/ff/ff.h"
#include "radmat/utils/printer.h"
#include "semble/semble_semble.h"
#include "ensem/ensem.h"
#include <sstream>



namespace radmat
{



  namespace 
  {
    struct dimension_printer
    {
      static void print(const std::string &msg)
      {}
      //      { std::cout << msg << std::endl; }
    };

    template<typename T>
      std::string to_string( T t )
      {
        std::stringstream ss; 
        ss << t ;
        return ss.str(); 
      }


    template<typename T>
      typename SEMBLE::PromoteEnsemVec<T>::Type
      get_ensem_row(const int row, const SEMBLE::SembleMatrix<T> &in)
      {
        typename SEMBLE::PromoteEnsemVec<T>::Type out;
        out.resize(in.getB());
        out.resizeObs(in.getM());
        for(int elem = 0; elem < in.getM(); ++elem)
          ENSEM::pokeObs(out,in.getEnsemElement(row,elem),elem);

        return out;
      }

    void write_jackfile_fit_report ( const TinsFitter &fits, const std::string &bpth)
    {
      std::vector<std::string> keys = fits.ff_ids(); 
      std::vector<std::string>::const_iterator it; 
      for(it = keys.begin(); it != keys.end(); ++it)
      {
        std::stringstream ss;
        ss << bpth << *it << ".jack";
        ENSEM::write(ss.str(),fits.getFF(*it));  
      }

      std::stringstream ss; 
      ss << bpth << "Q2.jack";
      ENSEM::write(ss.str(),fits.getQ2()); 
    }

    void write_single_jackfile_fit_report ( const TinsFitter &fits, const std::string &bpth, const std::string  &ff)
    {
      std::stringstream ss;
      ss << bpth << ff<< ".jack";
      ENSEM::write(ss.str(),fits.getFF(ff));  

      std::stringstream s; 
      s << bpth << "Q2.jack";
      ENSEM::write(s.str(),fits.getQ2()); 
    }

  } // anonomyous 



  RadmatSingleQ2Driver::RadmatSingleQ2Driver(void)
  {
    init_false(); 
  }


  RadmatSingleQ2Driver::RadmatSingleQ2Driver(const RadmatSingleQ2Driver &o)
    : init_linear_system(o.init_linear_system) , init_fits(o.init_fits) ,
    linear_system(o.linear_system) , fit_across_time(o.fit_across_time)
  { rot_id = std::string(""); }


  RadmatSingleQ2Driver& RadmatSingleQ2Driver::operator=(const RadmatSingleQ2Driver &o)
  {
    if(this != &o)
    {
      init_linear_system = o.init_linear_system;
      init_fits = o.init_fits;
      linear_system = o.linear_system;
      fit_across_time = o.fit_across_time;
      rot_id = o.rot_id; 
    }
    return *this; 
  }

  bool RadmatSingleQ2Driver::load_llsq(const rHandle<LLSQLatticeMultiData> &d, 
      const double pole_mass_squared,
      const double tolerance, 
      const bool mix_irreps)
  {
    if(!!!linear_system.load_data(d,tolerance))
      return false;

    if(linear_system.peek_tags().empty())
    {
      std::cerr << __func__ << ": warning, no tags" << std::endl; 
      return false;
    }

    if( (pole_mass_squared + SEMBLE::toScalar(ENSEM::mean(linear_system.Q2())) < 0.))
    {
      std::cout << "Killing Q2 = " << SEMBLE::toScalar(ENSEM::mean(linear_system.Q2())) 
        << "  q2tag = " << linear_system.qsq_label() << "  beacause pole_mass^2 + Q^2 < 0 "
        << "\nvalue is " << pole_mass_squared + SEMBLE::toScalar(ENSEM::mean(linear_system.Q2())) << std::endl;
      return false; 
    }

    // append the label
    append_rotation_group_label(rotation_group_label(mix_irreps)); 


    init_linear_system = true; 

    SEMBLE::SEMBLEIO::makeDirectoryPath(base_path() + std::string("llsq"));
    linear_system.dump_llsq_lattice(base_path() + std::string("llsq/"));
    linear_system.save_llsq_state( base_path() + std::string("llsq/") ); 
    return true;
  }

  bool RadmatSingleQ2Driver::load_llsq(const rHandle<LLSQLatticeMultiData> &d, 
      const double tolerance,
      const bool save_state)
  {
    if(!!!linear_system.load_data(d,tolerance))
      return false;

    if(linear_system.peek_tags().empty())
    {
      std::cerr << __func__ << ": warning, no tags" << std::endl; 
      return false;
    }

    init_linear_system = true; 

    if(save_state)
    {
      SEMBLE::SEMBLEIO::makeDirectoryPath(base_path() + std::string("llsq"));
      linear_system.dump_llsq_lattice(base_path() + std::string("llsq/"));
    }
    return true;
  }

  void RadmatSingleQ2Driver::solve_llsq(const std::string &soln_ID)
  {
    check_exit_linear_system();
    std::cout << "Solving Q2 = " << linear_system.qsq_label() 
      << "  " << linear_system.peek_tags().begin()->mom_string() << std::endl;

    linear_system.solve_llsq(soln_ID);
    init_solved_llsq = true; 
  }

  void RadmatSingleQ2Driver::save_llsq_state(void) const
  {
    check_exit_solved_llsq(); 
    linear_system.save_ff_state( base_path() ); 
  }

  void RadmatSingleQ2Driver::save_ff_of_t(void) const
  {
    check_exit_solved_llsq(); 
    SEMBLE::SembleMatrix<std::complex<double> > FF_of_t = linear_system.peek_FF(); 
    std::string path;
    path = base_path() + std::string("ff_of_t/");
    SEMBLE::SEMBLEIO::makeDirectoryPath( path ); 

    for(int row = 0; row < FF_of_t.getN(); ++row)
    {
      std::stringstream id;
      id << path << "unphasedFF_" << row << ".jack"; 
      ENSEM::write( id.str() , get_ensem_row(row,FF_of_t)); 
    }
  }

  FormFacSolutions<std::complex<double> > 
    RadmatSingleQ2Driver::grab_ff_solution(void) const
  {
    check_exit_solved_llsq(); 
    return linear_system.grab_ff_solution(); 
  }


  void RadmatSingleQ2Driver::fit_data(const ThreePointComparatorProps_t &fit_props, 
      const int tsrc,
      const int tsnk)
  {
    check_exit_linear_system();
    check_exit_solved_llsq(); 

    std::cout << "Fitting Q2 = " << linear_system.qsq_label() 
      << "  " << linear_system.peek_tags().begin()->mom_string() << std::endl;

    SEMBLE::SembleMatrix<std::complex<double> > FF_of_t = linear_system.peek_FF(); 
    LLSQComplexFormFactorData_t tmp; 

    printer_function<dimension_printer>( "FF_of_t B="
        + to_string( FF_of_t.getB() ) );
    printer_function<dimension_printer>( "FF_of_t N="
        + to_string( FF_of_t.getN() ) );
    printer_function<dimension_printer>( "FF_of_tM=" 
        + to_string( FF_of_t.getM() ) );

    // load up data
    for(int row = 0; row < FF_of_t.getN(); ++row)
      tmp.insert(linear_system.ff_id(row),FF_of_t.getRow(row));

    tmp.Qsq = linear_system.Q2(); 
    std::string pth = base_path() + std::string("t_ins_fits/");
    SEMBLE::SEMBLEIO::makeDirectoryPath(pth);

    // run the fit
    fit_across_time.fit(pth,tmp,fit_props,tsrc,tsnk);
    init_fits = true; 
  }

  std::pair<ENSEM::EnsemReal,ENSEM::EnsemReal>
    RadmatSingleQ2Driver::fit_and_dump_single_ffs(const ThreePointComparatorProps_t &fit_props, 
        const FormFacSolutions<std::complex<double> > &ff_soln,
        const int tsrc,
        const int tsnk,
        const std::string &ff,
        const FitParValue &v) const
    {
      LLSQComplexFormFactorData_t tmp; 

      // load up data
      bool found = false; 
      for(int row = 0; row < ff_soln.FF_t.getN(); ++row)
      {
        if( ff_soln.Names[row] == ff ) 
        {
          tmp.insert(ff,ff_soln.FF_t.getRow(row));
          found = true; 
          break; 
        }
      }

      // std::cout << "fitting " << ff << std::endl;
      // std::cout << ff_soln.FF_t.mean() << std::endl;

      ThreePointDataTag t = *(ff_soln.Ingredients.begin()); 
      ENSEM::EnsemReal this_Q2 = t.Q2(); 

      if( !!! found )
      {
        std::cout << "couldnt find " << ff << " did you mean one of:" << std::endl;
        std::vector<std::string>::const_iterator nit; 
        for(nit = ff_soln.Names.begin(); nit != ff_soln.Names.end(); ++nit)
          std::cout << *nit << std::endl;
        std::cout << "skipping" << std::endl; 
        return std::make_pair(ENSEM::Real(0.)*this_Q2,ENSEM::Real(0.)*this_Q2);
      }

      std::string pth = std::string("");

      // run the fit
      TinsFitter lcl_fit_across_time;
      lcl_fit_across_time.fit(pth,tmp,fit_props,tsrc,tsnk,true,v);

      // overwrite fit result with component fit 
      lcl_fit_across_time.writeFitPlotsWithComponents(pth); 

      // report fits 
      lcl_fit_across_time.writeFitLogs(pth);
      ENSEM::write("Q2.jack",this_Q2);

      ENSEM::EnsemReal this_ff = lcl_fit_across_time.getFF(ff); 
      std::stringstream res; 
      res << "1 q2 " 
        << SEMBLE::toScalar(ENSEM::mean(this_Q2)) << " " 
        << sqrt(SEMBLE::toScalar(ENSEM::variance(this_Q2))); 
      res << " ff " 
        << SEMBLE::toScalar(ENSEM::mean(this_ff)) << " " 
        << sqrt(SEMBLE::toScalar(ENSEM::variance(this_ff))); 

      DataRep3pt dr = t.data_rep; 
      std::pair<mom_t,mom_t> c_mom; 
      c_mom = radmat::LatticeRotationEnv::rotation_group_can_mom(t.left_mom,t.right_mom);

      res << " | " << t.full_irrep_id(dr,dr.l) 
        << " " << c_mom.first[0]
        << " " << c_mom.first[1]
        << " " << c_mom.first[2];

      res << " " << t.full_irrep_id(dr,dr.r) 
        << " " << c_mom.second[0]
        << " " << c_mom.second[1]
        << " " << c_mom.second[2];

      res << " g " << t.full_irrep_id(dr,dr.g); 
      res << std::endl; 

      std::string fout; 
      fout = ff + ".refit.dat"; 
      std::ofstream out(fout.c_str());
      out << res.str(); 
      out.close(); 


      return std::make_pair(this_ff,this_Q2); 
    }


  void RadmatSingleQ2Driver::chisq_analysis(const int tlows, const int thighs)
  {
    check_exit_fits();

    int tlow = tlows; 
    int thigh = thighs; 

    // screw it, lets just always run the chisq since it is the correct thing to do
    //    std::vector<std::string> ids = fit_across_time.ff_ids(); 
    //    std::vector<std::string>::const_iterator it; 
    //    for(it = ids.begin(); it != ids.end(); ++it)
    //    {
    //      rHandle<FitThreePoint> some_fit = fit_across_time.getFit(*it); 
    //      if ( some_fit->tlow() > tlow) 
    //        tlow = some_fit->tlow(); 
    //
    //      if ( some_fit->thigh() < thigh)
    //        thigh = some_fit->thigh();
    //
    //    }

    if ( thigh < tlow ) 
    {
      std::cout << __func__ << ": range [" << tlow << "," << thigh << "]" << std::endl; 
      std::cout << "thigh < tlow, unable to continue " << std::endl; 
      return ; 
    }

    std::cout << __func__ << ": using range [" << tlow << "," << thigh << "]" << std::endl; 

    SEMBLE::SEMBLEIO::makeDirectoryPath(base_path() + std::string("chisq/"));
    linear_system.chisq_analysis( fit_across_time.fetchFF().second,
        base_path() + std::string("chisq/") , tlow , thigh,1e-6);
  }


  ENSEM::EnsemReal RadmatSingleQ2Driver::Q2(void) const
  {
    check_exit_linear_system();
    return linear_system.Q2(); 
  }


  std::pair<ENSEM::EnsemReal, SEMBLE::SembleVector<double> > 
    RadmatSingleQ2Driver::fetchFF(void) const
    { 
      check_exit_fits(); 
      return fit_across_time.fetchFF(); 
    }


  RadmatSingleQ2Solution RadmatSingleQ2Driver::fetchSolution(void) const
  {
    check_exit_fits(); 

    RadmatSingleQ2Solution ret; 

    // yoink the data from the fit across time class 
    std::pair<ENSEM::EnsemReal,std::map<std::string,ENSEM::EnsemReal> > d; 
    d = fit_across_time.fetchNamedFF(); 

    // bung it into this silly solution thingy 
    ret.Q2 = d.first; 
    ret.ff_map = d.second; 
    ret.tags = linear_system.peek_tags(); 

    return ret; 
  }


  std::string RadmatSingleQ2Driver::tags_at_this_Q2(void) const
  {
    check_exit_linear_system(); 
    double qq = SEMBLE::toScalar(ENSEM::mean(Q2())); 
    std::vector<ThreePointDataTag> tt = linear_system.peek_tags(); 
    std::vector<ThreePointDataTag>::const_iterator it;
    std::stringstream ss;

    for(it = tt.begin(); it != tt.end(); ++it)
      ss << qq << " (tag val = " <<  linear_system.qsq_sort() << ") " << it->file_id << std::endl;

    return ss.str(); 
  }

  std::string RadmatSingleQ2Driver::rotation_group_label(const bool mix_irreps) const
  {
    // check_exit_linear_system(); 
    std::vector<ThreePointDataTag> tt = linear_system.peek_tags(); 

    if( tt.empty() )
      return std::string(); 

    // have to do this complicated thing so that when we split on cubic 
    // irreps we dont have the individual single q2 drivers writing to the 
    // same directory, this effectively makes the output path a function 
    // of the data and the sorting method  
    std::stringstream id; 
    id << tt.begin()->rot_qsq_tag(mix_irreps); 
    id << "__" + radmat::LatticeRotationEnv::rotation_group_label(tt.begin()->left_mom,tt.begin()->right_mom);

    return id.str(); 
  }

  void RadmatSingleQ2Driver::dump_fits(void) 
  {
    check_exit_fits();
    SEMBLE::SEMBLEIO::makeDirectoryPath(base_path() + std::string("fit_logs/"));
    fit_across_time.writeFitLogs(base_path() + std::string("fit_logs/"));
    SEMBLE::SEMBLEIO::makeDirectoryPath(base_path() + std::string("component_fits/")); 
    fit_across_time.writeFitPlotsWithComponents(base_path() + std::string("component_fits/"));
    SEMBLE::SEMBLEIO::makeDirectoryPath(base_path() + std::string("jack_files/")); 
    write_jackfile_fit_report ( fit_across_time, base_path() + std::string("jack_files/"));
  }




  void RadmatSingleQ2Driver::dump_llsq(void)
  {
    check_exit_linear_system();
    SEMBLE::SEMBLEIO::makeDirectoryPath(base_path() + std::string("llsq"));
    linear_system.dump_llsq(base_path() + std::string("llsq/")); 
  }



  void  RadmatSingleQ2Driver::check_exit(const bool b, const char *c) const 
  {
    if(!!!b)
    {
      std::cerr << __func__ << ": error: called by " << c << ", exiting." << std::endl;
      exit(1); 
    }
  }


  void RadmatSingleQ2Driver::init_false(void)
  { 
    init_linear_system = false;
    init_fits = false;
    rot_id = std::string(""); 
  }


  std::string RadmatSingleQ2Driver::base_path(void) const
  {
    std::stringstream ss; 
    ss << SEMBLE::SEMBLEIO::getPath() << "Q2_" << linear_system.qsq_sort();

    // this is an ugly naming convention 
    if( rot_id != std::string("") )
    {
      SEMBLE::SEMBLEIO::makeDirectoryPath(ss.str());
      ss << "/" << rot_id; 
    }

    ss << "/";
    SEMBLE::SEMBLEIO::makeDirectoryPath(ss.str());
    return ss.str(); 
  }


}



