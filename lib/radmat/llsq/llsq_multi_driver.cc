/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : llsq_multi_driver.cc

 * Purpose :

 * Creation Date : 22-02-2013

 * Last Modified : Fri 03 Oct 2014 05:14:27 PM EDT

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/


#include "llsq_multi_driver.h"
#include "llsq_solvers.h"
#include "llsq_chisq.h"
#include "llsq_multi_data_serialize.h"
#include "llsq_solution.h"
#include "radmat/utils/pow2assert.h"
#include "radmat/ff_interface/formfactor_kinematic_factors.h"
#include "radmat/utils/printer.h"
#include "ensem/ensem.h"
#include <complex>
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <utility>
#include "radmat/utils/printer.h"

#include "adat/handle.h"
#include "jackFitter/ensem_data.h"
#include "jackFitter/jackknife_fitter.h"
#include "jackFitter/ensem_data.h"
#include "jackFitter/three_point_fit_forms.h"
#include "jackFitter/plot.h"
// #define DEBUG_AT_MAKE_MOM_INV_TAGS
// #define DEBUG_AT_ZERO_SORTING



namespace radmat
{

  namespace 
  {
    // pull the unique elems of the union -- guard with map 
    rHandle<LLSQLatticeMultiData>
      unique_data(const rHandle<LLSQLatticeMultiData> &lattice_data,
          const LLSQLatticeMultiData &zeroed_data)
      {
        rHandle<LLSQLatticeMultiData> ret(new LLSQLatticeMultiData() ); 

        std::map<std::string,int> data_map; 
        std::vector<LLSQLatticeMultiData::TT> tags = zeroed_data.tags(); 
        for(int i = 0; i < tags.size(); ++i)
        {
          std::string key = tags[i].file_id; 
          if(data_map.find(key) == data_map.end())
          {
            data_map[key] = 1; 
            ret->append_row_semble( zeroed_data.get_row_semble(i),tags[i]); 
          }
          else
          {
            std::cout << __func__ << ": warning duped data??" << std::endl;
          }
        }
        
        tags = lattice_data->tags(); 

        // add unique elems from lattice data 
        for(int i = 0; i < tags.size(); ++i)
        {
          std::string key = tags[i].file_id; 
          if(data_map.find(key) == data_map.end())
          {
            data_map[key] = 1; 
            ret->append_row_semble( lattice_data->get_row_semble(i),tags[i]); 
          }
          else
          {
            std::cout << __func__ << ": warning duped data??" << std::endl;
          }
        }

        return ret; 
      } 

    void init_dim(SEMBLE::SembleMatrix<std::complex<double> > &to_init, 
        const SEMBLE::SembleVector<std::complex<double> > &first_row)
    {
      SEMBLE::SembleMatrix<std::complex<double> > foo(first_row.getB(),1,first_row.getN());
      foo.zeros(); 
      for(int elem = 0; elem < first_row.getN(); ++elem)
        foo.loadEnsemElement(0,elem,first_row.getEnsemElement(elem)); 

      to_init = foo; 
    }



    std::string sort_string(const ThreePointDataTag &t)
    {
      std::stringstream ss; 
      ss << t.left_mom[0] << t.left_mom[1] << t.left_mom[2] 
        << t.right_mom[1] << t.right_mom[1] << t.right_mom[2] 
        << t.gamma_row; 
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


    template<typename T>
      void my_writer_mean(const std::string &fname, const SEMBLE::SembleMatrix<T> &in)
      {
        itpp::Mat<T> mean = itpp::round_to_zero( in.mean() , 1e-5 );
        std::ofstream out(fname.c_str());
        out << mean ;
        out.close(); 
      }

    template <typename T>
      std::string number_to_string ( T Number )
      {
        std::stringstream ss;
        ss << Number;
        return ss.str();
      }

    template<typename T>
      void my_writer_rows(const std::string &fname, const SEMBLE::SembleMatrix<T> &in)
      {
        const int nr = in.getN();

        for( int row = 0; row < nr; ++ row)
        {
          ENSEM::write(fname + number_to_string(row) + std::string(".jack"), 
              get_ensem_row(row,in)); 
        }
      }

    void my_writer_cont_expr(const std::string &fname, const std::vector<ThreePointDataTag> &tt)
    {
      std::ofstream out(fname.c_str());
      for(unsigned int idx = 0; idx < tt.size(); ++idx)
        out << idx << " " << tt[idx].file_id << std::endl;
      out.close(); 
    }

    template<typename T>
      void my_writer(const std::string &fname, const T &t)
      {
        std::ofstream out(fname.c_str());
        out << t;
        out.close(); 
      }


    // take the fit form factors and compute the chisq of the fit across time 
    std::string chisq_per_t(const SEMBLE::SembleVector<std::complex<double> > &FF,
        const rHandle<LLSQLatticeMultiData> &lattice_data,
        const SEMBLE::SembleMatrix<std::complex<double> > &K,
        const double tol)
    {
      SEMBLE::SembleMatrix<std::complex<double> > data = lattice_data->data(); 
      SEMBLE::SembleVector<std::complex<double> > thy = K * FF; 
      std::stringstream ss; 
      ss << "#t chisq/DoF nDoF \n";  // think # makes gnuplot ignore the line..

      const int Lt = data.getM();

      for(int t = 0; t < Lt; ++t)
      {
        std::pair<double,int> chisq = ChisqAndDoF(data.getCol(t),thy,tol);
        ss << t << " " << chisq.first / double(chisq.second) << " " << chisq.second << std::endl;
      }

      return ss.str(); 
    }


    std::string chisq_per_data_of_fit_range(const SEMBLE::SembleVector<std::complex<double> > &FF,
        const rHandle<LLSQLatticeMultiData> &lattice_data,
        const SEMBLE::SembleMatrix<std::complex<double> > &K,
        const double tol,
        const int tlow,
        const int thigh)
    {
      const int Lt = thigh - tlow; 
      SEMBLE::SembleMatrix<std::complex<double> > data = lattice_data->data();
      std::vector<ThreePointDataTag> tags = lattice_data->tags(); 
      SEMBLE::SembleVector<std::complex<double> > thy = K * FF; 
      const int sz = tags.size(); 
      std::stringstream ss;
      ss << "#id chisq/DoF nDoF \n";

      for(int corr = 0; corr < sz; ++corr)
      {
        SEMBLE::SembleVector<std::complex<double> > work_thy,work_data(data.getB(),Lt); 
        work_thy = work_data; 
        for(int t = 0; t < Lt; ++t)
        {
          work_thy.loadEnsemElement(t,thy.getEnsemElement(corr));
          work_data.loadEnsemElement(t,data.getEnsemElement(corr,t+tlow));
        }
        std::pair<double,int> chisq = ChisqAndDoF(work_data,work_thy,tol);
        ss << tags[corr].file_id << " " << chisq.first / double(chisq.second) 
          << " " << chisq.second << std::endl;
      }

      return ss.str(); 
    }


    std::string chisq_of_system_of_fit_range(const SEMBLE::SembleVector<std::complex<double> > &FF,
        const rHandle<LLSQLatticeMultiData> &lattice_data,
        const SEMBLE::SembleMatrix<std::complex<double> > &K,
        const double tol,
        const int tlow,
        const int thigh)
    {
      const int Lt = thigh - tlow; 
      SEMBLE::SembleMatrix<std::complex<double> > data = lattice_data->data();
      SEMBLE::SembleVector<std::complex<double> > thy = K * FF; 
      std::stringstream ss;
      const int ncorr = data.getN(); 
      const int rankV = ncorr * Lt; 
      SEMBLE::SembleVector<std::complex<double> > wdata(data.getB(),rankV),wthy(data.getB(),rankV); 

      // construct vectors for chisq
      for(int corr = 0; corr < ncorr; ++corr)      
        for(int t = 0; t < Lt; ++t)
        {
          wdata.loadEnsemElement(corr*Lt + t, data.getEnsemElement(corr,t+tlow));
          wthy.loadEnsemElement(corr*Lt + t, thy.getEnsemElement(corr)); 
        }

      std::pair<double,int> chisq = ChisqAndDoF(wdata,wthy,tol); 
      ss<< "tlow = " << tlow << " thigh = " << thigh << " tol = " << tol << std::endl;
      ss << "chisq = " << chisq.first / double(chisq.second) << " nDoF = " << chisq.second;
      return ss.str(); 
    }

    struct is_solveable_printer
    {
      static void print(const std::string &s)
      {}
      // { std::cout << "is_solveable_printer " + s << std::endl;} 
    };

    bool is_llsq_solveable( const itpp::Mat<std::complex<double> > &m, const double tolerance)
    {
      itpp::Mat<std::complex<double> > M = itpp::hermitian_transpose(m) * m ; 

      //  printer_function<is_solveable_printer>( "in rows " + number_to_string(m.rows()) ); 
      //  printer_function<is_solveable_printer>( "in cols " + number_to_string(m.cols()) ); 
      //  printer_function<is_solveable_printer>( "M.rows() " + number_to_string(M.rows()) ); 
      //  printer_function<is_solveable_printer>( "M.cols() " + number_to_string(M.cols()) ); 

      std::stringstream ss; 
      ss << "singular values[ ";

      itpp::Vec<double> s = itpp::svd(M); 
      int non_zero_singular_values(0); 
      for(int i = 0; i < s.size(); ++i)
      {
        ss << s(i) << ", ";
        if( s(i) > tolerance )         
          ++non_zero_singular_values; 
      }

      ss << "]\n";


      int possible_extraction(0); 
      int rows = m.rows(); 
      int cols = m.cols(); 

      for(int i = 0; i < rows; ++i)
      {
        int count(0); 
        for(int j = 0; j < cols; ++j)
          if( std::norm(m(i,j)) > tolerance )
            ++count; 

        possible_extraction = ( count > possible_extraction ) ? count : possible_extraction; 
      }

       std::cout << __func__ << ": " << ss.str() << "looking to extract " 
         << possible_extraction << " N non zero singular values = " << non_zero_singular_values << std::endl;

      //  printer_function<is_solveable_printer>( "nnz sings " + number_to_string(non_zero_singular_values)); 
      // printer_function<is_solveable_printer>( "possible sings " + number_to_string(possible_extraction)); 

      return non_zero_singular_values >= possible_extraction; 
    }


  } // anonomyous 





  LLSQMultiDriver_t::LLSQMultiDriver_t(void) 
  {
    init_false(); 
  }


  bool LLSQMultiDriver_t::load_data(const rHandle<LLSQLatticeMultiData> &d,
      const double tolerance)
  {
    init_false();
    lattice_data = d;
    POW2_ASSERT(&*d);
    init_lat = true; 

    sort_data(); 
    bool success = run_zero_filter(tolerance);
    return success; 
  }


  void LLSQMultiDriver_t::splash_tags(void) const
  {
    lattice_data->splash_tags(); 
  }


  void LLSQMultiDriver_t::sort_data(void)
  {
    rHandle<LLSQLatticeMultiData> sorted_data(new LLSQLatticeMultiData); 
    std::vector<ThreePointDataTag> tags = lattice_data->tags(); 
    std::map<std::string,std::vector<int> > elems;
    std::map<std::string,std::vector<int> >::iterator it; 

    unsigned int sz = tags.size(); 
    for(unsigned int elem = 0; elem < sz; ++elem)
    {
      std::string tt = sort_string(tags[elem]);
      it = elems.find(tt);
      if(it == elems.end())
      {
        elems[tt] = std::vector<int>(1,elem); 
      }
      else
      {
        it->second.push_back(elem); 
      }
    }


    for(it = elems.begin(); it != elems.end(); ++it)
    {
      std::vector<int>::const_iterator sorted; 
      for(sorted = it->second.begin(); sorted != it->second.end(); ++sorted)
        sorted_data->append_row_semble(lattice_data->get_row_semble(*sorted),tags[*sorted]);
    }

    lattice_data = sorted_data; 
  }



  bool LLSQMultiDriver_t::run_zero_filter(const double tolerance)
  {
    check_exit_lat(); 

    rHandle<LLSQLatticeMultiData> non_zero_data(new LLSQLatticeMultiData);
    std::vector<ThreePointDataTag> old_tags;
    SEMBLE::SembleVector<std::complex<double> > Zero; 
    old_tags = lattice_data->tags(); 

    const unsigned int sz = old_tags.size(); 

    if ( sz == 0 ) 
    {
      std::cerr << __func__ << "no tags?? " << std::endl; 
      return false; 
    }


    FFKinematicFactors_t Kt;

    Zero = Kt.genFactors(&(*old_tags.begin())); 
    Zero.zeros(); 


    bool first = true; 
    SEMBLE::SembleVector<std::complex<double> > workV;

#ifdef DEBUG_AT_ZERO_SORTING
    bool first_zero = true; 
    SEMBLE::SembleMatrix<T> Kzero;
#endif 

    // loop to determine what should be zero
    for(unsigned int elem = 0; elem < sz; ++elem)
    {
      // pull down this set of kinematic factors
#ifdef DEBUG_AT_ZERO_SORTING 
      SEMBLE::SembleMatrix<T> KMat; 
      KMat = Kt.genFactorsMat(&old_tags[elem]); 
      std::cout << __func__ << ": tag type " << (&old_tags[elem])->type() << std::endl; 
      std::cout << __func__ << ": tag -> " << old_tags[elem].mom_string() << std::endl;
      std::cout << __func__ << ": tag -> data_rep " << old_tags[elem].rot_qsq_tag(false) << std::endl;
      std::cout << __func__ << ": tag -> file_id " << old_tags[elem].file_id << std::endl;
      std::cout << __func__ << ": Kinematic Matrix \n" << SEMBLE::mean(KMat) << std::endl;

      workV = Kt.genFactors(&old_tags[elem]); 
      std::cout << __func__ << ": pre round " << SEMBLE::mean(workV) << std::endl;
      workV = SEMBLE::round_to_zero(workV,tolerance); 
      std::cout << __func__ << ": post round " << SEMBLE::mean(workV) << std::endl;
#else
      workV = SEMBLE::round_to_zero(
          Kt.genFactors(&old_tags[elem]), tolerance); 

#endif 

      // if it is zero push it into the zeroed data pile
      if(workV == Zero)
      {
        zeroed_data.append_row_semble(
            lattice_data->get_row_semble(elem),
            old_tags[elem]);    

#ifdef DEBUG_AT_ZERO_SORTING
        if(first_zero)
        {
          init_dim(Kzero,workV); 
          first_zero = false; 
        }
        else
          Kzero.append_row(workV); 
#endif 

      }  
      else
      {
        // initialize dimensions on the first append 
        if(first)
        {
          init_dim(K,workV); 
          first = false; 
        }
        else
          K.append_row(workV);

        // throw it into the good data pile 
        non_zero_data->append_row_semble(
            lattice_data->get_row_semble(elem),
            old_tags[elem]);
      }
    }

    // we build the kinematic decomp along the way, save that bit of work 
    init_K = true; 

    // update lattice data to that good stuff 
    lattice_data = non_zero_data;

    // sanity, turn this off later  
#ifdef DEBUG_AT_ZERO_SORTING
    if(first)
    {
      std::cout << __func__ << ": everything you passed me was zero, " 
        << "using tolerance = " << tolerance << std::endl;
      std::cout << __func__ << ": K-mean(zeroed data) \n" << SEMBLE::mean(Kzero) << std::endl;
    }
#endif

    // check that we can solve the llsq -- only move if we have some data 
    bool solveable = (!!!first) ?  is_llsq_solveable( SEMBLE::mean( peek_K() ) , 1e-5) : false; 

    // warn that we are killing this data point 
    if(!!!solveable)
    {
      std::cout << __func__ << ": not enough data points to solve the llsq" << std::endl;
      std::cout << "for " << old_tags.begin()->mom_string() << std::endl; 
    }

    return solveable;
  }



  // wrapper for the inversion method 
  void LLSQMultiDriver_t::solve_llsq(const std::string &soln_ID)
  {
    check_exit_lat();

    rHandle<LLSQBaseSolver_t<std::complex<double> > >
      foo = LLSQSolverFactoryEnv::callFactory(soln_ID);

    bool success(true); 

    if(foo->invertable())
      success &= solve_fast(soln_ID); 
    else
      success &= solve_slow(soln_ID);

    check_exit(success,__func__);
  }


  // legacy, this is old 
  void LLSQMultiDriver_t::chisq_analysis(const SEMBLE::SembleVector<double> &ff,
      const std::string &path,
      const int tlow,
      const int thigh,
      const double tol)
  {
    check_exit_lat();
    check_exit_K();
    check_exit_Kinv();
    check_exit_FF(); 



    // itpp/semble dont like mixed operations
    SEMBLE::SembleVector<std::complex<double> > ff_cmplx(ff.getB(),ff.getN());
    const int B = ff.getB();
    const int N = ff.getN();

    for(int bin = 0; bin < B; ++bin)
      for(int elem = 0; elem < N; ++elem)
        ff_cmplx.setElement(bin,elem,std::complex<double>(ff[bin][elem],0.));

    std::string chisq_per_tt = chisq_per_t(ff_cmplx,lattice_data,K,tol);
    std::string chisq_per_data = chisq_per_data_of_fit_range(ff_cmplx,
        lattice_data,K,tol,tlow,thigh);
    std::string chisq_of_sys = chisq_of_system_of_fit_range(ff_cmplx,
        lattice_data,K,tol,tlow,thigh); 

    std::stringstream ss,per_t,per_data,sys;
    ss << path; 
    per_t << ss.str() << "chisq_per_t.txt";
    per_data << ss.str() << "chisq_per_data.txt";
    sys << ss.str() << "chisq_of_sys_of_fit_range.txt"; 

    my_writer(per_t.str(),chisq_per_tt);
    my_writer(per_data.str(),chisq_per_data);
    my_writer(sys.str(),chisq_of_sys); 
  }


  void LLSQMultiDriver_t::dump_llsq_lattice(const std::string &path)
  {
    check_exit_lat();
    std::stringstream b,cont,zero,db;
    cont << path << "row_index_to_continuum_elem.txt";
    zero << path << "zeroed_matrix_elems";
    b << path << "lattice_mat_elem_";

    my_writer_rows(b.str(), lattice_data->data());
    my_writer_cont_expr(cont.str(), lattice_data->tags());

    SEMBLE::SEMBLEIO::makeDirectoryPath(zero.str()); 
    zero << "/zeroed_mat_elem_"; 
    my_writer_rows(zero.str(),zeroed_data.data()); 
    zero.str(std::string()); 
    zero << path << "zeroed_matrix_elems/row_index_to_continuum_elem.txt"; 
    my_writer_cont_expr(zero.str(),zeroed_data.tags()); 
  }



  void LLSQMultiDriver_t::dump_llsq(const std::string &path)
  {
    check_exit_K();
    check_exit_Kinv();
    check_exit_FF(); 

    std::stringstream Kpth; 
    Kpth <<  path << "K";
    SEMBLE::SEMBLEIO::makeDirectoryPath( Kpth.str() ); 
    Kpth << "/"; 



    std::stringstream ss, A,Ainv, x; 
    std::stringstream unity, Kstr,Kistr; 
    ss << path;
    A << ss.str() << "K.mean";
    Ainv << ss.str() << "Kinv.mean"; 
    unity << ss.str() <<  "Kinv_x_K.mean"; 
    x << ss.str() << "ff_";
    Kstr << Kpth.str() << "K_row_";
    Kistr << Kpth.str() << "Kinv_row_";


    my_writer_mean(A.str(),K);
    my_writer_mean(Ainv.str(),Kinv);
    my_writer_mean(unity.str(),Kinv*K); 
    my_writer_rows(x.str(), FF_t); 
    my_writer_rows(Kstr.str(),K); 
    my_writer_rows(Kistr.str(),Kinv); 

    std::string pth = path + std::string("solver.log"); 
    std::ofstream out( pth.c_str() ); 
    out << solver_log; 
    out.close(); 
  }



  void LLSQMultiDriver_t::save_llsq_state(const std::string &path) const
  {
    //    std::cout << __func__ << ": entering" << std::endl; 
    check_exit_lat();

    std::stringstream sslat,ssz; 
    sslat << path << "state_database.rad";
    ssz << path << "/zeroed_matrix_elems"; 
    SEMBLE::SEMBLEIO::makeDirectoryPath(ssz.str()); 
    ssz << "/zeroed_state_database.rad";

    std::stringstream all_data_name; 
    all_data_name << path << "all_data.rad"; 

    rHandle<LLSQLatticeMultiData> all_data = unique_data(lattice_data,zeroed_data); 


    ADATIO::BinaryFileWriter binlat(sslat.str()) , binz(ssz.str()), binall(all_data_name.str()); 
    write(binlat,*lattice_data); 
    write(binz,zeroed_data); 
    write(binall,*all_data); 

    binlat.close(); 
    binz.close(); 
    binall.close(); 

  }


  void LLSQMultiDriver_t::save_ff_state(const std::string &path) const
  {
    //    std::cout << __func__ << ": entering" << std::endl; 
    check_exit_lat();
    std::stringstream ss; 
    ss << path << "ff_database.rad"; 
    ADATIO::BinaryFileWriter bin(ss.str());

    std::vector<std::string> ffnames ;
    std::vector<ThreePointDataTag> fftags ; 
    SEMBLE::SembleMatrix<std::complex<double> > ffs; 

    ffs = peek_FF(); 
    fftags = peek_tags(); 

    for(int i = 0; i < ffs.getN(); ++i)
      ffnames.push_back(ff_id(i)); 

    FormFacSolutions<std::complex<double> >
      ff(ffs,fftags,ffnames); 

    write(bin,ff); 
    bin.close();
  }

  FormFacSolutions<std::complex<double> > LLSQMultiDriver_t::grab_ff_solution(void) const
  {
    check_exit_lat();

    std::vector<std::string> ffnames ;
    std::vector<ThreePointDataTag> fftags ; 
    SEMBLE::SembleMatrix<std::complex<double> > ffs; 

    ffs = peek_FF(); 
    fftags = peek_tags(); 

    for(int i = 0; i < ffs.getN(); ++i)
      ffnames.push_back(ff_id(i)); 

    FormFacSolutions<std::complex<double> >
      ff(ffs,fftags,ffnames); 

    return ff; 
  }

  SEMBLE::SembleMatrix<std::complex<double> > 
    LLSQMultiDriver_t::get_rephase_FF(void) const
    {
      typedef std::complex<double> T; 
      SEMBLE::SembleMatrix<T> un = peek_FF(); 
      ENSEM::EnsemVectorComplex in; 
      in = get_ensem_row(0,un); 

      ENSEM::EnsemVectorReal real, imag;

      double const_real, const_imag; 

      real = ENSEM::real(in);
      imag = ENSEM::imag(in);

      std::vector<double> t;
      for(int i = 0; i < in.numElem(); ++i)
        t.push_back(double(i)); 

      // take the middle 70 %
      int thigh = int( double(t.size()) * .85); 
      int tlow = int( double(t.size()) *.15);  
      EnsemData ereal(t,real),eimag(t,imag); 
      ereal.hideDataAboveX(thigh - 0.1); 
      eimag.hideDataBelowX(tlow -0.1); 

      ADAT::Handle<FitFunction> freal(new ThreePointConstant), fimag(new ThreePointConstant);  
      JackFit fit_real(ereal,freal), fit_imag(eimag,fimag); 

      fit_real.runAvgFit(); 
      fit_imag.runAvgFit(); 
      const_real = fit_real.getAvgFitParValue(0);
      const_imag = fit_imag.getAvgFitParValue(0);

      // apparently this guy is zero, just dump the whole thing 
      if( (const_real < 2e-2) && (const_imag < 2e-2) ) 
        return un;

      // something went wrong, deal with it higher up 
      if ( isnan(const_real) && isnan(const_imag))
        return un; 

      double fit_phase = std::arg(std::complex<double>(const_real,const_imag)); 
      double apply_phase = -fit_phase; 

      // apply this phase to everyone!
      std::complex<double> phase(cos(apply_phase),sin(apply_phase));  
      un *= phase; 

      return un; 
    }


  bool LLSQMultiDriver_t::solve_fast(const std::string &soln_ID)
  {
    check_exit_lat();    
    rHandle<LLSQBaseSolver_t<std::complex<double> > > 
      my_solver = LLSQSolverFactoryEnv::callFactory(soln_ID);
    POW2_ASSERT(&*my_solver); 

    bool success(true); 

    generate_kinematic_factors();

    Kinv = my_solver->inv(K); 
    solver_log = my_solver->solution_log(); 
    init_Kinv = true; 

    FF_t = Kinv * lattice_data->data(); 

    // pull out the name list here -- yuck
    rHandle<FormFactorBase_t> KK = 
      FormFactorDecompositionFactoryEnv::callFactory(
          lattice_data->tags().begin()->mat_elem_id);

    ff_ids = KK->ff_ids(); 

    init_FF = true; 

    /*
       std::cout << __func__ << std::endl;
       std::cout << "K = \n" << K.mean() << std::endl;
       std::cout << "Kinv = \n" << Kinv.mean() << std::endl;
       std::cout << "lat = \n" << SEMBLE::mean(lattice_data->data()) << std::endl;
       std::cout << "FF = \n" << FF_t.mean() << std::endl;
       */

    return success; 
  }


  bool LLSQMultiDriver_t::solve_slow(const std::string &soln_ID)
  {
    bool success(true); 
    check_exit(false,__func__); 
    return success;
  }


  void LLSQMultiDriver_t::generate_kinematic_factors(void)
  {
    check_exit_lat(); 

    // can do this step during the zero filter 
    if( !!! init_K )
    {
      std::vector<ThreePointDataTag> tags = lattice_data->tags(); 
      std::vector<ThreePointDataTag>::const_iterator it; 

      FFKinematicFactors_t KK;

      for(it = tags.begin(); it != tags.end(); ++it)
      {

        SEMBLE::SembleVector<std::complex<double> > work;
        work = KK.genFactors( &(*it) );

        if(it == tags.begin())
          init_dim(K,work); 
        else
          K.append_row(work);
      }

      init_K = true; 
    }
  }

  void LLSQMultiDriver_t::init_false(void)
  { 
    init_lat = false;
    init_K = false;
    init_Kinv = false;
    init_FF = false;
  }


  void LLSQMultiDriver_t::check_exit(const bool b, const char *c) const
  {
    if(!!!b)
    {
      std::cerr << __func__ << ": error: called by " << c << ", exiting." << std::endl;
      exit(1); 
    }
  }







}


#undef DEBUG_AT_MAKE_MOM_INV_TAGS
