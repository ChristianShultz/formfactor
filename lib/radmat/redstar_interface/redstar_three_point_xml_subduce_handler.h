#ifndef REDSTAR_THREE_POINT_XML_SUBDUCE_HANDLER_H
#define REDSTAR_THREE_POINT_XML_SUBDUCE_HANDLER_H 

#include "redstar_three_point_xml_handler.h"
#include "radmat/utils/stringify.h"


namespace radmat
{

  struct RedstarThreePointXMLSubduceHandler;
  REGISTER_STRINGIFY_TYPE( RedstarThreePointXMLSubduceHandler );  

  struct RedstarThreePointXMLSubduceHandler
    : public RedstarThreePointXMLHandler
  {
    virtual ~RedstarThreePointXMLSubduceHandler() {}

    virtual std::string type() const 
    { return Stringify<RedstarThreePointXMLSubduceHandler>(); }

    virtual std::vector<ThreePointData>
      handle_work(const RedstarThreePointXMLInput &); 
  }; 


} // radmat 


#endif /* REDSTAR_THREE_POINT_XML_SUBDUCE_HANDLER_H */
