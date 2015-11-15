#ifndef REDSTAR_THREE_POINT_XML_LORENTZ_HANDLER_H
#define REDSTAR_THREE_POINT_XML_LORENTZ_HANDLER_H 

#include "redstar_three_point_xml_handler.h"
#include "radmat/utils/stringify.h"

namespace radmat
{

  struct RedstarThreePointXMLLorentzHandler;
  REGISTER_STRINGIFY_TYPE( RedstarThreePointXMLLorentzHandler ); 

  struct RedstarThreePointXMLLorentzHandler
    : public RedstarThreePointXMLHandler
  {
    virtual ~RedstarThreePointXMLLorentzHandler() {}

    virtual std::string type() const
    { return Stringify<RedstarThreePointXMLLorentzHandler>(); }

    virtual std::vector<ThreePointData>
      handle_work(const RedstarThreePointXMLInput &);
  }; 
  

} // radmat 


#endif /* REDSTAR_THREE_POINT_XML_LORENTZ_HANDLER_H */
