/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

* File Name : redstar_abstract_xml_factory.cc

* Purpose :

* Creation Date : 12-11-2013

* Last Modified : Wed 30 Apr 2014 02:09:37 PM EDT

* Created By : shultz

_._._._._._._._._._._._._._._._._._._._._.*/


#include "redstar_abstract_xml_factory.h"
#include "radmat/utils/pow2assert.h"
#include "radmat/utils/stringify.h"
#include "radmat/utils/printer.h"

#include <exception>

#include "redstar_three_point_xml_lorentz_handler.h"
#include "redstar_three_point_xml_mixed_handler.h"
#include "redstar_three_point_xml_subduce_handler.h"
#include "redstar_single_particle_meson_xml_interface.h"
#include "redstar_vector_current_xml_interface.h"
#include "redstar_improved_vector_current_xml_interface.h"

namespace FacEnv = radmat::TheRedstarAbstractXMLFactoryEnv;
typedef radmat::TheRedstarAbstractXMLFactory Factory; 

namespace radmat
{

  namespace TheRedstarAbstractXMLFactoryEnv
  {
    bool local_registration = false; 

    namespace
    {
      template<typename T> 
        AbsRedstarXMLInterface_t* 
        callback(void)
        {
          AbsRedstarXMLInterface_t *foo = new T; 
          return foo; 
        }
    }



    template<typename T> 
      bool do_reg(void)
      {
        bool foo = Factory::Instance().registerObject(Stringify<T>(),callback<T>); 
        if (!!!foo)
          std::cerr << __PRETTY_FUNCTION__ 
            << ":error reging " << Stringify<T>() << std::endl;

        return foo; 
      }


    bool registerAll(void)
    {
      printer_function<debug_print_reg_all>(__PRETTY_FUNCTION__); 
      bool success = true; 
      if ( !!! local_registration ) 
      {
        success &= do_reg<RedstarSingleParticleMesonXML>();
        success &= do_reg<RedstarVectorCurrentXML>(); 
        success &= do_reg<RedstarThreePointXMLLorentzHandler>(); 
        success &= do_reg<RedstarThreePointXMLSubduceHandler>(); 
        success &= do_reg<RedstarThreePointXMLMixedHandler>(); 
        success &= do_reg<RedstarImprovedVectorCurrentXML>(); 
        local_registration = true; 
      }

      if(!!! success)
        throw std::string("reg error in TheRedstarAbstractBlockFactoryEnv"); 

      return success; 
    }

    rHandle<AbsRedstarXMLInterface_t>
      callFactory(const std::string &id)
      {
        AbsRedstarXMLInterface_t *foo; 
        try
        {
          foo = Factory::Instance().createObject(id);  
        }
        catch(...)
        {
          POW2_ASSERT(false); 
        }

        POW2_ASSERT(foo); 

        return rHandle<AbsRedstarXMLInterface_t>(foo); 
      }

  } // TheRedstarAbstractBlockFactoryEnv

} // radmat


