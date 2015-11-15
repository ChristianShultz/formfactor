#ifndef REDSTAR_THREE_POINT_XML_MIXED_HANDLER_H
#define REDSTAR_THREE_POINT_XML_MIXED_HANDLER_H 


#include "redstar_three_point_xml_handler.h"
#include "radmat/utils/stringify.h"

namespace radmat
{

  struct RedstarThreePointXMLMixedHandler;
  REGISTER_STRINGIFY_TYPE( RedstarThreePointXMLMixedHandler ); 

  struct RedstarThreePointXMLMixedHandler
    : public RedstarThreePointXMLHandler
  {
    virtual ~RedstarThreePointXMLMixedHandler() {}

    virtual std::string type() const
    { return Stringify<RedstarThreePointXMLMixedHandler>(); }

    virtual std::vector<ThreePointData>
      handle_work(const RedstarThreePointXMLInput &);
  }; 
  

} // radmat 




#endif /* REDSTAR_THREE_POINT_XML_MIXED_HANDLER_H */
