#ifndef RADMAT_DRIVER_PROPS_H
#define RADMAT_DRIVER_PROPS_H

#include <string>
#include "io/adat_xmlio.h"
#include "jackFitter/three_point_fit_forms.h"
#include "radmat/construct_data/construct_correlators_xml.h"

namespace radmat
{
  struct RDriverProps_t
  {
    ThreePointComparatorProps_t threePointComparatorProps;
    ThreePointCorrIni_t threePointIni;  
    
    std::string chisq;
    int version; 
    int bigPhase;                         // an overall plus or minus 
    int maxThread;                        // # < 1 does nothing , #>= 1 sets nthread to # 
    double poleMass;                      // does the analytic continuation exist?  this is the lightest vector state
    double tolerance; 
  };                                       //  or the branch point at the energy of two pions -- nb this is the mass not
                                           //  the square of the mass


  std::string toString (const RDriverProps_t &);
  std::ostream& operator<<(std::ostream& , const RDriverProps_t &);
  void read(ADATXML::XMLReader &xml, const std::string &path, RDriverProps_t &); 

}// namespace radmat


#endif /* RADMAT_DRIVER_PROPS_H */
