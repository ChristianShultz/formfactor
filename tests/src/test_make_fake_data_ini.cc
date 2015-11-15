/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : test_make_fake_data_ini.cc

 * Purpose :

 * Creation Date : 13-07-2012

 * Last Modified : Tue Jul 17 16:19:31 2012

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/

#include "../headers/test_common.h"
#include "radmat/fake_data/fake_data_ini.h"
#include <string>


namespace radmat
{

  FakeDataIni_t makeFakeIni(void)
  {
    FakeDataIni_t ini;

    // stateProps 
    ini.stateProps.sameOp = true;
    ini.stateProps.readZ = false;

    ini.stateProps.mProps.sourceMasses.resize(3);
    ini.stateProps.mProps.sourceMasses[0] = 0.1;
    ini.stateProps.mProps.sourceMasses[1] = 0.2;
    ini.stateProps.mProps.sourceMasses[2] = 0.3;
    ini.stateProps.mProps.sourceVarO = 0.1;
    ini.stateProps.mProps.sourceUpdateCovariance = true;
    ini.stateProps.mProps.sourceUpdateVariance = true;
    ini.stateProps.mProps.sinkMasses = ini.stateProps.mProps.sourceMasses;
    ini.stateProps.mProps.sinkVarO = ini.stateProps.mProps.sourceVarO;
    ini.stateProps.mProps.sinkUpdateCovariance = ini.stateProps.mProps.sourceUpdateCovariance;
    ini.stateProps.mProps.sinkUpdateVariance = ini.stateProps.mProps.sourceUpdateVariance;
    
    ini.stateProps.zProps.overlapGenerator = std::string("unitary");
    ini.stateProps.zProps.suppress = true;
    ini.stateProps.zProps.targetZ_at_order1 = true;
    ini.stateProps.zProps.suppressionOrder = 0.1;
    ini.stateProps.zProps.updateCovariance = true;
    ini.stateProps.zProps.updateVariance = true;
    ini.stateProps.zProps.varianceOrder = 0.1;


    // matElemProps
    ini.matElemProps.diag = std::string("PiPi");
    ini.matElemProps.off = ini.matElemProps.diag;
      ini.matElemProps.left_target = 0;
    ini.matElemProps.right_target = ini.matElemProps.left_target;

    // timeProps
    ini.timeProps.tsource = 0;
    ini.timeProps.tsink = 25;

    // dataProps
    ini.dataProps.ncfg = 200;
    ini.dataProps.hel_source.resize(1);
    ini.dataProps.hel_source[0] = 0;
    ini.dataProps.hel_sink = ini.dataProps.hel_source;
    ini.dataProps.momenta.resize(4);
    ini.dataProps.momenta[0].momSource.resize(3);
    ini.dataProps.momenta[0].momSource[0] = 0;
    ini.dataProps.momenta[0].momSource[1] = 0;
    ini.dataProps.momenta[0].momSource[2] = 0;
    ini.dataProps.momenta[0].momSink = ini.dataProps.momenta[0].momSource;
    ini.dataProps.momenta[1] = ini.dataProps.momenta[0];
    ini.dataProps.momenta[2] = ini.dataProps.momenta[0];
    ini.dataProps.momenta[3] = ini.dataProps.momenta[0];
    ini.dataProps.momenta[0].momSource[2] = 1;
    ini.dataProps.momenta[1].momSink = ini.dataProps.momenta[0].momSource;
    ini.dataProps.momenta[2].momSource[0] = 1;
    ini.dataProps.momenta[3].momSink = ini.dataProps.momenta[2].momSource;


    // dispersionProps
    ini.dispersionProps.a_t_inverse = 5.6;
    ini.dispersionProps.xi = 3.444;
    ini.dispersionProps.L_s = 2;


    return ini;
  }


}
