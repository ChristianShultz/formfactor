#ifndef REDSTAR_SINGLE_PARTICLE_MESON_XML_INTERFACE_H
#define REDSTAR_SINGLE_PARTICLE_MESON_XML_INTERFACE_H 

#include "redstar_abstract_xml_interface.h"
#include "radmat/utils/stringify.h"

  // 
  //
  // SINGLE PARTICLE MESON XML INTERFACE
  //
  //

namespace radmat
{

  struct RedstarSingleParticleMesonXML; 
  REGISTER_STRINGIFY_TYPE(RedstarSingleParticleMesonXML); 


  struct RedstarSingleParticleMesonXML
    : public AbsRedstarXMLInterface_t
  {

    virtual std::string type(void) const 
    {
      return std::string(Stringify<RedstarSingleParticleMesonXML>()); 
    }

    virtual void read(ADATXML::XMLReader &xml, 
        const std::string &path); 

    virtual std::string write(void) const;  

    std::string cont_rep;                               // J0p , J1m etc 
    ADATXML::Array<int> H;                              // helicity 
    bool fill_star;                                     // consider all rotations?
    ADATXML::Array< ADATXML::Array<int> > mom;          // momentum 
    int twoI_z;                                         // twice Isospin 
    std::string name;                                   // particle stem 
    bool creation_op;                                   // is it a creation operaor
    bool smearedP;                                      // is it a smeared operator
    bool isProjected;                                   // is it a projected operator
    int t_slice;                                        // t_slice
  }; 


}




#endif /* REDSTAR_SINGLE_PARTICLE_MESON_XML_INTERFACE_H */
