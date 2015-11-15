#ifndef REDSTAR_PARTICLE_HANDLER_UTILS_H
#define REDSTAR_PARTICLE_HANDLER_UTILS_H 


#include "hadron/hadron_npart_npt_corr.h"
#include "radmat/utils/obj_expr_t.h"
#include "ensem/ensem.h"
#include <string>

#include "redstar_single_particle_meson_xml_interface.h"
#include "redstar_vector_current_xml_interface.h"
#include "redstar_improved_vector_current_xml_interface.h"

namespace radmat
{

  // the list object expressions that we will be playing with, 
  // these represent linear combinations of operators
  typedef ListObjExpr_t<ENSEM::Complex,
          Hadron::KeyHadronNPartNPtCorr_t::NPoint_t> EnsemRedstarBlock;

  // this is a primitive data tag
  // origin_rep is where it lives in lorentz world
  // data_rep is where it lives ( cubic or lorentz )
  // row is the row of the data rep 
  // data is a list of coefficient times NPoint_t 
  struct BlockData
  {
    std::string origin_rep;
    std::string data_rep;
    int row; 
    EnsemRedstarBlock data; 
  };



  
  std::vector<BlockData>
    generate_lorentz_block( const RedstarSingleParticleMesonXML * const );
  
  std::vector<BlockData>
    generate_cubic_block( const RedstarSingleParticleMesonXML * const );

  std::vector<BlockData>
    generate_lorentz_block( const RedstarVectorCurrentXML * const );
  
  std::vector<BlockData>
    generate_cubic_block( const RedstarVectorCurrentXML * const );
  
  struct VectorCurrentImprovedBlockData
  {
    std::string origin_rep;
    std::string data_rep;
    int row; 
    EnsemRedstarBlock data; 
    RedstarImprovedVectorCurrentXML::improvement imp; 
  };

  std::vector<VectorCurrentImprovedBlockData>
    generate_lorentz_block( const RedstarImprovedVectorCurrentXML *const ); 

  std::vector<VectorCurrentImprovedBlockData>
    generate_cubic_block( const RedstarImprovedVectorCurrentXML *const); 

} 


#endif /* REDSTAR_PARTICLE_HANDLER_UTILS_H */
