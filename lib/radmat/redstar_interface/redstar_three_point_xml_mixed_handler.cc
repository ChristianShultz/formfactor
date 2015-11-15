/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

* File Name : redstar_three_point_xml_mixed_handler.cc

* Purpose :

* Creation Date : 23-04-2014

* Last Modified : Wed 30 Apr 2014 02:17:17 PM EDT

* Created By : shultz

_._._._._._._._._._._._._._._._._._._._._.*/


#include "redstar_three_point_xml_mixed_handler.h"
#include "redstar_three_point_utils.h"
#include "radmat/utils/utils.h"

namespace radmat
{



  std::vector<ThreePointData>
    RedstarThreePointXMLMixedHandler::handle_work(const RedstarThreePointXMLInput &inp) 
    {
      ADATXML::Array<AbstractNamedObject<AbsRedstarXMLInterface_t> > xml_int; 
      std::string ensemble; 

      xml_int = this->get_npt(); 
      ensemble = this->get_ensemble(); 



      // grab left 
      const RedstarSingleParticleMesonXML * lefty;
      POW2_ASSERT( xml_int[0].param->type() == Stringify<RedstarSingleParticleMesonXML>() );      
      lefty = dynamic_cast< const RedstarSingleParticleMesonXML * >( xml_int[0].param.get_ptr() ); 
      std::vector<BlockData> left = generate_cubic_block(lefty);

      // grab right 
      const RedstarSingleParticleMesonXML * righty;
      POW2_ASSERT( xml_int[2].param->type() == Stringify<RedstarSingleParticleMesonXML>() );      
      righty = dynamic_cast< const RedstarSingleParticleMesonXML * >( xml_int[2].param.get_ptr() ); 
      std::vector<BlockData> right = generate_cubic_block(righty); 

      ADATXML::Array<int> timesliz(3); 
      timesliz[0] = lefty->t_slice; 
      timesliz[2] = righty->t_slice; 

      // figure out what to do with the photon 
      if( xml_int[1].param->type() == Stringify<RedstarVectorCurrentXML>() )
      {
        const RedstarVectorCurrentXML *gammay; 
        gammay = dynamic_cast< const RedstarVectorCurrentXML * >( xml_int[1].param.get_ptr() ); 
        std::vector<BlockData> gamma = generate_lorentz_block(gammay);
        timesliz[1] = gammay->t_slice; 
        this->set_timeslice(timesliz); 
        return merge_blocks( left, gamma, right, ensemble ); 
      }
      else if( xml_int[1].param->type() == Stringify<RedstarImprovedVectorCurrentXML>() )
      {
        const RedstarImprovedVectorCurrentXML *gammay; 
        gammay = dynamic_cast<const RedstarImprovedVectorCurrentXML *>(xml_int[1].param.get_ptr());
        std::vector<VectorCurrentImprovedBlockData> gamma = generate_lorentz_block(gammay); 
        timesliz[1] = gammay->t_slice; 
        this->set_timeslice(timesliz); 
        return merge_blocks(left,gamma,right,ensemble,inp); 
      }
      else
      {
        POW2_ASSERT(false); 
      }

      // make compiler happy 
      return std::vector<ThreePointData>(); 
    }


} // radmat 


