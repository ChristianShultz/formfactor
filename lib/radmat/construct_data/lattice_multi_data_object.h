#ifndef LATTICE_MULTI_DATA_OBJECT_H
#define LATTICE_MULTI_DATA_OBJECT_H

#include "radmat/llsq/llsq_multi_data.h"
#include "radmat/data_representation/data_tag_three_point.h"

namespace radmat
{

  // the definition of the linear system we are creating
  typedef LLSQMultiData<ThreePointDataTag,std::complex<double> > LLSQLatticeMultiData; 

} // radmat

#endif /* LATTICE_MULTI_DATA_OBJECT_H */
