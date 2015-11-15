/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : llsq_formfactor_data.cc

 * Purpose :

 * Creation Date : 21-02-2014

 * Last Modified : Mon 27 Oct 2014 01:39:46 PM EDT

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/

#include "llsq_formfactor_data.h"
#include "ensem/ensem.h"
#include "jackFitter/plot.h"
#include <string>
#include <math.h>
#include <vector>
#include <iostream>
#include <complex>
#include "radmat/utils/splash.h"
#include "radmat/utils/pow2assert.h"
#include "radmat/utils/printer.h"
#include "adat/handle.h"

#include "jackFitter/ensem_data.h"
#include "jackFitter/jackknife_fitter.h"
#include "jackFitter/three_point_fit_forms.h"



namespace radmat
{


  namespace 
  {
    // swapping in empty functions should make the compiler 
    // optimize these print statements away.. plus its fun
    struct dimension_printer
    {
      static void print(const std::string &msg)
      {}
      //      { std::cout << msg << std::endl; }
    };

    struct inp_dimension_printer
    {
      static void print(const std::string &msg)
      {}
      //      { std::cout << msg << std::endl; }
    };

    struct ret_dimension_printer
    {
      static void print(const std::string &msg)
      // {}
       { std::cout << msg << std::endl; }
    };

    struct mean_printer
    {
      static void print(const std::string &s)
      {}
      // { std::cout << s << std::endl;}
    };

    struct case_printer
    {
      static void print(const std::string &s)
       { std::cout << s << std::endl;}
    };

    struct ensem_printer
    {
      static void print(const std::string &s)
      {}
      // {std::cout << s << std::endl;}
    };

    template<typename T>
      std::string to_string( T t )
      {
        std::stringstream ss; 
        ss << t ;
        return ss.str(); 
      }



    enum PHASE
    {
      RP,
      RM,
      IP,
      IM,
      ZERO,
      ERROR
    };

    inline double rad(const double deg)
    {
      return 3.14159 * deg/180.;
    }

    // carpal tunnel prevention 
    typedef std::pair<PHASE,ENSEM::EnsemVectorReal> phase_pair;

    // run an avg fit on the real and imag bits
    //      to try to find a phase

    std::pair<phase_pair,std::string>
      convert_to_real(const ENSEM::EnsemVectorComplex & in, 
          const int tlow,
          const int thigh)
      {
        std::stringstream fit_log; 

        fit_log << "tlow " << tlow << "   thigh " << thigh << std::endl;

        ENSEM::EnsemVectorReal real, imag;

        real = ENSEM::real(in);
        imag = ENSEM::imag(in);

        //  be smarter and use a fit 
        //
        //        // mid plus offset (th-tl)/2 + tl
        //        int t_sample = int( double(thigh + tlow)/2. ); 
        //        double const_real = SEMBLE::toScalar( ENSEM::mean( ENSEM::peekObs(real,t_sample)));
        //        double const_imag = SEMBLE::toScalar( ENSEM::mean( ENSEM::peekObs(imag,t_sample)));
        //        double const_real_err = SEMBLE::toScalar( ENSEM::variance( ENSEM::peekObs(real,t_sample)));
        //        double const_imag_err = SEMBLE::toScalar( ENSEM::variance( ENSEM::peekObs(imag,t_sample)));


        std::vector<double> t; 
        for(int i = 0; i < in.numElem(); ++i)
          t.push_back(double(i)); 


        EnsemData ereal(t,real), eimag(t,imag); 

        ereal.hideDataAboveX(thigh + 0.1);
        ereal.hideDataBelowX(tlow -0.1);

        eimag.hideDataAboveX(thigh + 0.1);
        eimag.hideDataBelowX(tlow -0.1);

        // run avg fits
        ADAT::Handle<FitFunction> freal(new ThreePointConstant()), fimag(new ThreePointConstant());
        JackFit fit_real(ereal,freal) , fit_imag(eimag,fimag); 

        fit_real.runAvgFit();
        fit_imag.runAvgFit(); 

        // pull average 
        double const_real = fit_real.getAvgFitParValue(0);
        double const_real_err = fit_real.getAvgFitParError(0);

        double const_imag = fit_imag.getAvgFitParValue(0);
        double const_imag_err = fit_imag.getAvgFitParError(0);

        // sometimes needs some fine tuning
        double consistent_with_zero = 1.;

        fit_log << "consistency with zero set at " << consistent_with_zero << std::endl;
        fit_log << "const_real " << const_real 
          << "  +/-  " << const_real_err << std::endl;
        fit_log << "const_imag " << const_imag 
          << "  +/-  " << const_imag_err << std::endl;


        // what if it is a longitudial factor or crossing 
        //    we are only getting about 2% precision, 
        //    assume a FF of O(1)
        if( (fabs(const_real) < 2e-2) && (fabs(const_imag) < 2e-2) ) 
        {
          fit_log << "* decided too small to resolve automatically " << std::endl;
          fit_log << "making an arbitrary choice to return the larger" << std::endl;
          printer_function<case_printer>("ff is consistent with zero");

          std::stringstream corr_foor, corr_fooc; 

          for(int i = 0; i < in.numElem(); ++i)
          {
            ENSEM::EnsemReal foor,fooc; 
            foor = ENSEM::real(ENSEM::peekObs(in,i));
            fooc = ENSEM::imag(ENSEM::peekObs(in,i));

            corr_foor << i  << " " << ENSEM::toDouble(ENSEM::mean(foor)) << " " 
              << ENSEM::toDouble(ENSEM::sqrt(ENSEM::variance(foor))) << "\n"; 

            corr_fooc << i  << " " << ENSEM::toDouble(ENSEM::mean(fooc)) << " " 
              << ENSEM::toDouble(ENSEM::sqrt(ENSEM::variance(fooc))) << "\n"; 
          }

          fit_log << "corr.real = " << corr_foor.str() << std::endl; 
          fit_log << "corr.imag = " << corr_fooc.str() << std::endl;



          // return whichever is bigger?
          if( fabs(const_real) > fabs(const_imag) )
          {
            return std::make_pair(phase_pair(ZERO,real),fit_log.str());
          }
          else
          {
            return std::make_pair(phase_pair(ZERO,imag),fit_log.str());
          }
        }

        // is the constant consistent with zero? 
        if( fabs(const_real) - consistent_with_zero*fabs(const_real_err) < 0.)
        {

          // need to be careful about just returning imag here since it can 
          // muck up the phase check later, if in fact the correlator is consistent 
          // with zero via ERROR, not VALUE, then we can guard here 
          if( fabs(const_imag) - consistent_with_zero*fabs(const_imag_err) < 0.)
          {
            fit_log << "* decided overall consistent with zero" << std::endl;
            printer_function<case_printer>("ff is consistent with zero");
            printer_function<case_printer>("imag is consistent with zero");
            printer_function<case_printer>("imag = " + to_string(const_imag) 
                + "+/-" + to_string(sqrt(fabs(const_imag_err)))); 
            // doesn't matter, both are zero to precision  
            return std::make_pair(phase_pair(ZERO,real),fit_log.str());
          }


          fit_log << "* decided it was imag" << std::endl;
          printer_function<case_printer>("real is consistent with zero");
          printer_function<case_printer>("real = " + to_string(const_real) 
              + "+/-" + to_string(sqrt(fabs(const_real_err)))); 
          printer_function<case_printer>("imag = " + to_string(const_imag) 
              + "+/-" + to_string(sqrt(fabs(const_imag_err)))); 
          if (const_imag > 0. )
            return std::make_pair(phase_pair(IP,imag),fit_log.str()); 
          return std::make_pair(phase_pair(IM,imag),fit_log.str());
        }

        if( fabs(const_imag) - consistent_with_zero*fabs(const_imag_err) < 0.)
        {
          fit_log << "* decided it was real " << std::endl;
          printer_function<case_printer>("imag is consistent with zero");
          printer_function<case_printer>("imag = " + to_string(const_imag) 
              + "+/-" + to_string(sqrt(fabs(const_imag_err)))); 
          if( const_real > 0. )
            return std::make_pair(phase_pair(RP,real),fit_log.str()); 
          return std::make_pair(phase_pair(RM,real),fit_log.str());
        }


        // something unexpected happened
        if ( isnan(const_real) || isnan(const_imag))
        {
          fit_log << "foung nan, doing a dance " << std::endl;
          // the fitter seems to fail on ensembles with zero variance??
          // 
          // with svd resetting it is very easy to get an ensemble 
          // with mean zero and variance zero thus we have to 
          // work harder which makes christian cranky since its sunday 
          // 
          // the next time he read this it was a monday and he was still cranky

          // guard zero variance here.. 
          // this is completely nuts, someone needs to update the stupid fitter

          // need to use assignment operator, no builtin ensem constructors
          SEMBLE::SembleVector<double> foor; foor = real; 
          SEMBLE::SembleVector<double> fooi; fooi = imag; 

          // pull down the copies, grab the variance using semble
          // then compare it to zero, if the test passes check for 
          // either the real part or the imag part being explicitly 
          // zero, then return the other with the correct phase
          itpp::Vec<double> bar( foor.getN() ), bazr,bazi; 
          bar.zeros(); 
          bazr = foor.variance(); 
          bazi = fooi.variance(); 

          if( bazr == bar ) 
            if( bazi == bar )
            {
              printer_function<case_printer>("zero variance encountered"); 
              bazr = foor.mean(); 
              bazi = fooi.mean(); 
              if ( ( bazr == bar ) && ( bazi == bar ) )
                return std::make_pair(phase_pair(ZERO,real),fit_log.str()); 
              if ( bazr == bar )
                return bazi(0) > 0 ? std::make_pair(phase_pair(IP,imag),fit_log.str())
                  : std::make_pair(phase_pair(IM,imag),fit_log.str());
              if( bazi == bar )
                return bazr(0) > 0 ? std::make_pair(phase_pair(RP,real),fit_log.str()) 
                  : std::make_pair(phase_pair(RM,real),fit_log.str()); 
            }

          fit_log << " i quit " << std::endl;

          // otherwise die one function up since this is stupid 

          SPLASH("encountered nan: check bad_corr.nan.jack, bad_corr.nan.ax for the correlator"); 

          AxisPlot plot; 
          plot.addEnsemData(ENSEM::real(in),"\\sq",1);
          plot.addEnsemData(ENSEM::imag(in),"\\sq",2);
          plot.sendToFile("bad_corr.nan.ax");

          ENSEM::write("bad_corr.nan.jack",in); 

          return std::make_pair(phase_pair(ERROR,real),fit_log.str()); 
        }

        // move to simple things where we just check what axis the result is along 
        // and decide the phase based on the phase angle

        double fit_phase = std::arg(std::complex<double>(const_real,const_imag)); 

        if(fit_phase < -3.*3.14159/4.)
          fit_phase = - fit_phase; 


        double phase = fit_phase ; 

        fit_log << " general cases, found a phase of " << phase << std::endl;

        printer_function<case_printer>("general cases encountered"); 

        if( (phase < 0.174528) && (phase > -0.174708) ) // +/- 10 degree about 0 in rad
        {
          printer_function<case_printer>("RP encountered"); 
          fit_log << " decided RP " << std::endl;
          return std::make_pair(phase_pair(RP,real),fit_log.str());
        }
        else if( (phase > 1.39622) && (phase < 1.74528)) // 90deg
        {
          printer_function<case_printer>("IP encountered"); 
          fit_log << " decided IP " << std::endl;
          return std::make_pair(phase_pair(IP,imag),fit_log.str());
        }
        else if( (phase > 2.96697) || (phase < -2.96715))  //180deg // this is atan2 specific, it returns (-pi,pi)
        {
          printer_function<case_printer>("RM encountered"); 
          fit_log << " decided RM " << std::endl;
          return std::make_pair(phase_pair(RM,real),fit_log.str());
        }
        else if( (phase > -1.74546) && (phase < -1.3964)) // 270 deg
        {
          printer_function<case_printer>("IM encountered"); 
          fit_log << " decided IM " << std::endl;
          return std::make_pair(phase_pair(IM,imag),fit_log.str());
        }
        else
        {
          fit_log << " got a baddie " << std::endl;
          printer_function<case_printer>("unknown phase encountered"); 
          std::cout << "The calculated phase was " << phase*180./3.14159 << " (deg)" << std::endl;
          std::cout << "rl = " << const_real << " im = " << const_imag << std::endl;
          std::cout << "used tlow = " << tlow << " thigh = " << thigh << std::endl; 
          SPLASH("check bad_corr.jack, bad_corr.ax for the correlator, returning zero"); 

          AxisPlot plot; 
          plot.addEnsemData(ENSEM::real(in),"\\sq",1);
          plot.addEnsemData(ENSEM::imag(in),"\\sq",2);
          plot.sendToFile("bad_corr.ax");

          ENSEM::write("bad_corr.jack",in); 

        }

        fit_log << " you really mucked this up " << std::endl;

        return std::make_pair( phase_pair( ERROR, SEMBLE::toScalar( double(0.) ) * real ), fit_log.str()); 
      }


    // find the phase for a single vector
    std::pair<phase_pair,std::string>
      find_phase( const SEMBLE::SembleVector<std::complex<double> > &d, 
          const int tlow, 
          const int thigh)
      {
        ENSEM::EnsemVectorComplex e; 
        int bns = d.getB(); 
        int sz = d.getN(); 
        e.resize(bns); 
        e.resizeObs(sz); 
        for(int i = 0; i < sz; ++i)
          ENSEM::pokeObs( e, d.getEnsemElement(i) , i);

        // run the fits on the central 70% of the correlator, this 
        // is a completely arbitrary choice 
        return convert_to_real( e, int(sz*0.15) , int(sz*0.85) ); 
      }

    // find the phase for all of the vectors
    std::pair< std::map<std::string, phase_pair> , std::string >
      find_phases( const LLSQComplexFormFactorData_t &d , 
          const int tlow, 
          const int thigh)
      {
        std::stringstream log; 

        std::map<std::string,phase_pair> ret; 
        LLSQComplexFormFactorData_t::const_iterator it; 
        for(it = d.begin(); it != d.end(); ++it)
        {
          printer_function<inp_dimension_printer>(
              "input " +  it->first 
              + " N " + to_string( it->second.getN() )
              + " B " + to_string( it->second.getB() ) ); 

          printer_function<mean_printer>(
              "input_corr " +  it->first 
              + to_string( it->second.mean() ) ); 

          log << "input: " << it->first << std::endl;

          std::cout << "finding phase on " << it->first << std::endl;

          std::pair<phase_pair,std::string> foo = find_phase(it->second,tlow,thigh); 

          std::cout << "found phase for " << it->first << std::endl;

          log << foo.second << "\n\n" << std::endl;

          ret.insert( std::make_pair(it->first, foo.first ) ); 
        }

        return std::make_pair(ret,log.str()); 
      }

    // check that all phases are consistent
    bool check_phases( const std::map<std::string, phase_pair> &mappy)
    {
      std::map<std::string,phase_pair>::const_iterator it; 
      for(it = mappy.begin(); it != mappy.end(); ++it)
        if( it->second.first == ERROR )
        {
          std::cout << __PRETTY_FUNCTION__ << ": error encountered, exiting" << std::endl;
          std::cout << it->first << " was bad" << std::endl;
          return false; 
        }

      // above guards ERROR
      PHASE expectedA,expectedB;  
      for(it = mappy.begin(); it != mappy.end(); ++it)
        if( it->second.first != ZERO )
          expectedA = it->second.first; 

      // so it was something
      if( expectedA == RP )
        expectedB = RM; 
      if( expectedA == RM )
        expectedB = RP; 
      if( expectedA == IP )
        expectedB = IM; 
      if( expectedA == IM )
        expectedB = IP; 

      bool pass = true; 
      // check that they are all either real , imag , or zero
      for(it = mappy.begin(); it != mappy.end(); ++it)
      {
        printer_function<dimension_printer>( "check_phases" + it->first 
            + " numElem = " + to_string(it->second.second.numElem()) 
            + " size = " + to_string(it->second.second.size()) ); 
        if( (it->second.first != expectedA )
            &&(it->second.first != expectedB) 
            &&(it->second.first != ZERO) )
        {
          std::cout << __PRETTY_FUNCTION__ 
            << "\nerror: encountered, unexpected phases, returning failure" << std::endl;
          std::cout << "expected " << expectedA << " " << expectedB 
            << " or " << ZERO << " got " << it->second.first << std::endl;
          std::cout << "keys \nRP->" << RP 
            << "\nRM->" << RM
            << "\nIP->" << IP
            << "\nIM->" << IM
            << "\nZZ->" << ZERO 
            << std::endl;
          std::cout << it->first << " was bad" << std::endl;
          pass = false; 
        }
      }

      return pass; 
    }

    // run checks then push the result into the return data
    LLSQRealFormFactorData_t
      do_work_local(const LLSQComplexFormFactorData_t &d, 
          const int tlow,
          const int thigh, 
          const std::string &fbase)
      {
        printer_function<inp_dimension_printer>("entering rephase"); 

        std::pair<std::map<std::string, phase_pair>,std::string> mappy = find_phases(d,tlow,thigh); 
        bool passed = check_phases(mappy.first); 

        if( !!! passed ) 
        {
          std::cout << __PRETTY_FUNCTION__ << ": received failure \n Q2 = " 
            << SEMBLE::toScalar(ENSEM::mean(d.Q2())) << " exiting.. " << std::endl;
          std::cout << "fbase = " << fbase << std::endl;
          std::cout << "logs \n\n" << mappy.second << std::endl; 
          std::cout << "dumping results to run dir " << std::endl;


          LLSQComplexFormFactorData_t::const_iterator corr_it; 

          // write out what we had
          for(corr_it = d.begin(); corr_it != d.end(); ++corr_it)
          {
            ENSEM::EnsemVectorComplex e; 
            int bns = corr_it->second.getB(); 
            int sz = corr_it->second.getN(); 
            e.resize(bns); 
            e.resizeObs(sz); 
            for(int i = 0; i < sz; ++i)
              ENSEM::pokeObs( e, corr_it->second.getEnsemElement(i) , i);
            ENSEM::write( corr_it->first , e );
          }

          exit(1); 
        }

        printer_function<ret_dimension_printer>(
            " map size " + to_string(mappy.first.size()) ); 

        LLSQRealFormFactorData_t ret; 

        // update momentum transfer        
        ret.Qsq = d.Q2(); 

        // load elements for fit
        std::map<std::string,phase_pair>::const_iterator it; 
        for(it = mappy.first.begin(); it != mappy.first.end(); ++it)
        {
          printer_function<ret_dimension_printer>( 
              "output " + it->first 
              + " numElem " + to_string(it->second.second.numElem()) 
              + " size " + to_string(it->second.second.size()));  
          ret.mappy[it->first] = it->second.second; 
        }

        // possibly print output result
        LLSQRealFormFactorData_t::const_iterator pit; 
        for(pit = ret.begin(); pit != ret.end(); ++pit)
          printer_function<mean_printer>(
              "output_corr " +  pit->first 
              + to_string( pit->second.mean() ) ); 

        printer_function<ret_dimension_printer>("exiting rephase");

        std::stringstream output_f;
        output_f << fbase + "fit_log.txt";
        std::ofstream out(output_f.str().c_str());
        out << mappy.second;
        out.close(); 

        return ret; 
      }


  } // anonymous



  // callback 
  LLSQRealFormFactorData_t 
    rephase_formfactor_data( const LLSQComplexFormFactorData_t &d, 
        const int tlow, 
        const int thigh,
        const std::string &fbase)
    {
      return do_work_local(d,tlow,thigh,fbase); 
    }

} // radmat
