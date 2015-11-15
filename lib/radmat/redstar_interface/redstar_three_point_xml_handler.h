#ifndef REDSTAR_THREE_POINT_XML_HANDLER_H
#define REDSTAR_THREE_POINT_XML_HANDLER_H 


#include "redstar_three_point_xml_interface.h"
#include "redstar_three_point_utils.h"
#include "ensem/ensem.h"
#include "hadron/hadron_npart_npt_corr.h"
#include "radmat/utils/stringify.h"
#include <vector>


namespace radmat
{


  struct RedstarThreePointXMLHandler;
  REGISTER_STRINGIFY_TYPE(RedstarThreePointXMLHandler);

  // generate data and primitive tag types
  struct RedstarThreePointXMLHandler
    : public RedstarThreePointXML
  {

    virtual std::string type() const
    {
      return Stringify<RedstarThreePointXMLHandler>();
    }

    // need to pass names in the case of improvement
    virtual std::vector<ThreePointData> handle_work(const RedstarThreePointXMLInput &) = 0; 
  }; 


} // radmat 



#endif /* REDSTAR_THREE_POINT_XML_HANDLER_H */
