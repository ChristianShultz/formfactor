/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name :

 * Purpose :

 * Creation Date : 17-07-2012

 * Last Modified : Thu Nov 29 16:40:54 2012

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/
#include "headers/test_common.h"
#include "headers/tester.h"
#include "radmat/fake_data/fake_3pt_function_aux.h"
#include "radmat/fake_data/fake_data_ini.h"
#include "semble/semble_meta.h"
#include <complex>
#include <string>
#include <sstream>
#include "ensem/ensem.h"
#include "semble/semble_semble.h"
#include "radmat/utils/aux.h"
#include "radmat/utils/harmonic.h"
#include "radmat/utils/pow2assert.h"
#include "radmat/utils/splash.h"
#include "radmat/ff/formfactor_factory.h"
#include "radmat/ff/ff_gen_llsq_row.h"
#include "adat/handle.h"
#include "radmat/fake_data/fake_spectrum.h"
#include "radmat/fake_data/fake_overlaps.h"
#include <vector>
#include <string>
#include <sstream>
#include <iostream>


#include "radmat/fake_data/fake_3pt_function_aux.h"
#include "radmat/fake_data/fake_3pt_function.h"
#include "radmat/fake_data/fake_data_ini.h"
#include "semble/semble_meta.h"
#include <complex>

#include "radmat/utils/splash.h"
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <exception>
#include "io/adat_xmlio.h"
#include "ensem/ensem.h"




namespace
{
  std::string n = std::string("\n");
  std::string blank = std::string("");
}



#define DEBUG_3PT_CJS 


#ifdef DEBUG_3PT_CJS

#include <iostream>
#include <vector>

using namespace radmat;

  template<typename T>
void print(const std::string &foo, const typename ADAT::Handle<FakeDataInputs_p<T> > &input)
{


  std::cout << foo << std::endl;

  std::cout << "zsource: " << std::endl;
  typename std::vector<SEMBLE::SembleVector<T> >::const_iterator it;

  int ct(0);  
  for(it = input->zsource.begin(); it != input->zsource.end(); it++, ct++)
    std::cout << "t = " << ct << n << it->mean() << std::endl;

  std::cout << "zsink: " << std::endl;

  ct = 0;
  for(it = input->zsink.begin(); it != input->zsink.end(); it++, ct++)
    std::cout << "t = " << ct << n << it->mean() << std::endl;

  std::vector<SEMBLE::SembleVector<double> >::const_iterator dit; 

  std::cout << "specsource: " << std::endl;

  ct = 0;
  for(dit = input->specsource.begin(); dit != input->specsource.end(); dit++,ct++)
    std::cout << "t = " << ct << n << dit->mean() << std::endl;

  std::cout << "specsink: " << std::endl;
  ct = 0;
  for(dit = input->specsink.begin(); dit != input->specsink.end(); dit++,ct++)
    std::cout << "t = " << ct << n << dit->mean() << n ;

  std::cout << "momfactor: " << input->mom_factor << std::endl;


}


  template<typename T>
void print(const std::string &foo, const typename ADAT::Handle<FakeDataInputs<T> > &input)
{
  std::cout << foo;
  std::cout << "working: " << std::endl;
  print(blank,input->working);
  std::cout << "orig: " << std::endl;
  print(blank,input->original);
}

#else

using namespace radmat;

  template<typename T>
void print(const std::string &,const typename ADAT::Handle<FakeDataInputs_p<T> > &)
{ }



  template<typename T>
void print(const std::string &, const typename ADAT::Handle<FakeDataInputs<T> > &)
{ }

#endif






using namespace radmat;
using namespace ENSEM;
using namespace ADAT;
using namespace ADATIO;

  int 
main(int argc, char *argv[])
{
  if(argc != 2)
  {
    SPLASH("usage: test_fake_ini : <xmlinifile> ");
    exit(1);
  }

  std::string xmlinifile;
  std::istringstream val(argv[1]);
  val >> xmlinifile;
  std::cout << "Loading xmlinifile: " << xmlinifile << std::endl;

  FakeDataIni_t fakeinikeys;

  try
  {
    XMLReader xml(xmlinifile);
    read(xml,"/FakeDataIni",fakeinikeys);
  }
  catch(std::exception &e)
  {
    std::cout << "std exception: " << e.what();
  }
  catch(std::string &e)
  {
    std::cout << __func__ << ": ERROR: can't read xmlinifile (" << xmlinifile << "): " << e << std::endl;
    exit(1);
  }
  catch(...)
  {
    SPLASH("An error occured while loading the fake data inifile");
    exit(1);
  }



  ADAT::Handle<FakeDataInputs_p<std::complex<double> > > 
    orig = generateOriginalInputs<std::complex<double> >(fakeinikeys);

  //   print(std::string("genOriginalInputs"),orig);


  ADAT::Handle<FakeDataInputs<std::complex<double> > > input = copyFakeInput(orig);

  //  print (std::string("copy") , input);

  applyZSuppression(input->working);

  //print(std::string("Zsup"), input);

  applyDispersion(input->working,fakeinikeys.dataProps.momenta[0]);

  //print(std::string("dispersion"), input);

  //  SEMBLE::PromoteEnsem<std::complex<double> >::Type foobar;
  //  foobar =  makeFakeDataPoint(input,fakeinikeys.dataProps.momenta[0],0,0,15,0);


  //////////////////////////////////////////////////////
  //////////////////////////////////////////////////////
  //////////////////////////////////////////////////////




#ifdef DEBUG_3PT_CJS

  //////////////////////////////////////////////////////
  //////////////////////////////////////////////////////
  //////////////////////////////////////////////////////




  ADAT::Handle<FakeDataInputs_p<std::complex<double> > > handle = input->working;

  ENSEM::EnsemVectorComplex tpt,matEnsem,factorEnsem;
  tpt.resize(handle->specsink[0].getB());
  tpt.resizeObs(handle->specsink.size());
  tpt = SEMBLE::toScalar(std::complex<double>(0.,0.));
  matEnsem = tpt;
  factorEnsem = tpt;

  for(int t_ins = handle->ini.timeProps.tsource; t_ins <= fakeinikeys.timeProps.tsink; t_ins++)
  {
    int lorentz(0);

    pProp_t mom = handle->ini.dataProps.momenta[0];
    ThreePtPropagationFactor<std::complex<double> > facGen;
    FakeMatrixElement matGen;

    const SEMBLE::SembleVector<std::complex<double> > *zsource, *zsink;
    const SEMBLE::SembleVector<double> *specsink, *specsource;
    const SEMBLE::SembleVector<double> *specsink_ins, *specsource_ins;




    zsource = &(handle->zsource[handle->ini.timeProps.tsource]);
    zsink = &(handle->zsink[handle->ini.timeProps.tsink]);
    specsource = &(handle->specsource[handle->ini.timeProps.tsource]);
    specsink = &(handle->specsink[handle->ini.timeProps.tsink]);
    specsource_ins = &(handle->specsource[t_ins]);
    specsink_ins = &(handle->specsink[t_ins]);


    const int nsource = specsource->getN();
    const int nsink = specsink->getN();
    ENSEM::EnsemComplex bar;

    bar.resize(specsink->getB());
    bar = SEMBLE::toScalar(std::complex<double>(0.,0.));

    ENSEM::EnsemComplex baz,boo;

    baz = bar;
    boo = bar;

    std::stringstream ss;
    ss << handle->ini.matElemProps.diag << "_" << 0 << "_" << 0;


    for(int source = 0; source < nsource; source++)
      for(int sink = 0; sink < nsink; sink++)
      {

        SemblePInv mom_inv =   makeMomInvariants( (*specsink_ins)(sink),
            (*specsource_ins)(source),
            mom.momSink,
            mom.momSource,
            handle->mom_factor);

               ENSEM::EnsemReal Q2 = mom_inv.Q2(); 

        ENSEM::EnsemComplex mat;
        mat = matGen(ss.str(),
            lorentz,
            mom_inv,
            handle->ffgenerator(source,sink),
            SEMBLE::toScalar(ENSEM::mean(Q2)));

        ENSEM::EnsemComplex factor;
        factor = facGen( (*specsink)(sink), (*zsink)(sink), handle->ini.timeProps.tsink,t_ins,
            (*specsource)(source), (*zsource)(source), handle->ini.timeProps.tsource);

        bar += (mat * factor);
        boo += mat;
        baz += factor;


        /*//print a load of crap
          std::cout << "source: " << source << "  sink: " << sink << std::endl;
          std::cout << "Q2: " << SEMBLE::toScalar(ENSEM::mean(Q2)) << std::endl;
          std::cout << "mat: " << SEMBLE::toScalar(ENSEM::mean(mat)) << std::endl;
          std::cout << "fac: " << SEMBLE::toScalar(ENSEM::mean(factor)) << std::endl; 
          std::cout << "mat*fac: " << SEMBLE::toScalar(ENSEM::mean(mat*factor)) << std::endl;
         */

      }
    //    __builtin_trap();

    pokeObs(tpt,bar,t_ins);
    pokeObs(matEnsem,boo,t_ins);
    pokeObs(factorEnsem,baz,t_ins);
  }

  std::cout << "three-pt" << std::endl;
  for(unsigned int t = 0; t < handle->specsink.size(); t++)
    std::cout << t << " " << ENSEM::mean(peekObs(tpt,t)) << std::endl;

  std::cout << "mat elem" << std::endl;
  for(unsigned int t = 0; t < handle->specsink.size(); t++)
    std::cout << t << " " << ENSEM::mean(peekObs(matEnsem,t)) << std::endl;

  std::cout << "factor elem" << std::endl;
  for(unsigned int t = 0; t < handle->specsink.size(); t++)
    std::cout << t << " " << ENSEM::mean(peekObs(factorEnsem,t)) << std::endl;


#endif


  return 0;
}





#undef  DEBUG_3PT_CJS


