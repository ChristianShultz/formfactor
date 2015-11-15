/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name :

 * Purpose :

 * Creation Date : 10-12-2013

 * Last Modified : Tue 06 May 2014 06:38:36 PM EDT

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/


#include "radmat/data_representation/data_representation_factory.h"
#include "radmat/data_representation/data_representation_subduction_map.h"
#include "radmat/ff/lorentzff_formfactor_factory.h"
#include "radmat/rotation_interface/rotation_group_generator.h"
#include "radmat/rotation_interface/Wigner_D_matrix_factory.h"
#include "radmat/rotation_interface/Wigner_D_matrix_manager.h"
#include "radmat/ff_interface/formfactor_factory.h"
#include "radmat/llsq/llsq_solvers.h"
#include "radmat/redstar_interface/redstar_abstract_xml_factory.h"
#include "radmat/redstar_interface/redstar_canonical_lattice_rotations.h"
#include "radmat/redstar_interface/redstar_canonical_rotations.h"
#include "radmat/redstar_interface/redstar_invert_subduction.h"
#include "radmat/spectrum/spectrum_factory.h"
#include "radmat/utils/printer.h"
#include <exception>

namespace radmat
{
  namespace AllFactoryEnv
  {

    namespace 
    {
      bool local_registration = false; 
    }


    bool registerAll(void)
    {
      printer_function<console_print>("Registering All Factories"); 
      bool success = true; 
      if(!!! local_registration ) 
      {

        try
        {
          // build subduction tables
          success &= radmat::InvertSubductionEnv::registerAll(); 
          // build lorentz helicity decompositions
          success &= radmat::LorentzffFormFactorDecompositionFactoryEnv::registerAll(); 
          // build  reps 
          success &= radmat::DataRepresentationFactoryEnv::registerAll(); 
          // register all subduction tables
          success &= radmat::TheSmarterSubduceTableMapFactoryEnv::registerAll();
          // build radmat's internal decompositions
          success &= radmat::FormFactorDecompositionFactoryEnv::registerAll(); 
          // build llsq solution classes
          success &= radmat::LLSQSolverFactoryEnv::registerAll(); 
          // redstar xml
          success &= radmat::TheRedstarAbstractXMLFactoryEnv::registerAll(); 
          // redstar canonical rotations 
          success &= radmat::CanonicalRotationEnv::registerAll(); 
          // redstar all lattice rotations
          success &= radmat::CanonicalLatticeRotationEnv::registerAll(); 
          // some more rotations
          success &= radmat::LatticeRotationEnv::registerAll(); 
          // d matricies 
          success &= radmat::WignerDMatrixEnv::registerAll(4); // up to J = 2 
          // pre allocate thread storage 
          success &= radmat::WignerThreadMapEnv::registerAll(); 
          // spectral info 
          success &= radmat::TheSpectrumFactoryEnv::registerAll( LATTICE::L_3F_743);
        }
        catch(std::exception &e)
        {
          std::cout << "std::exception: " << e.what() << std::endl; 
          std::cout << "dying in " << __PRETTY_FUNCTION__ << ":" 
            <<__FILE__ << ":" << __LINE__ << std::endl;
          exit(1); 
        }
        catch(std::string &s)
        {
          std::cout << "string exception: " << s << std::endl; 
          std::cout << "dying in " << __PRETTY_FUNCTION__ << ":" 
            <<__FILE__ << ":" << __LINE__ << std::endl;
          exit(1); 
        }
        catch(...)
        {
          std::cout << "some non standard non string exception" << std::endl; 
          std::cout << "dying in " << __PRETTY_FUNCTION__ << ":" 
            <<__FILE__ << ":" << __LINE__ << std::endl;
          exit(1); 
        }

        local_registration = true;  
      }

      return success; 
    }


  } 


}

