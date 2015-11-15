/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name :

 * Purpose :

 * Creation Date : 24-10-2012

 * Last Modified : Thu Nov 29 16:51:04 2012

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/


#include "../headers/tester.h"
#include "../headers/test_common.h"
// -- doesnt exist anymore #include "radmat/load_data/continuum_projection.h"
#include "hadron/hadron_npart_npt_corr.h"
#include "io/adat_xmlio.h"
#include "semble/semble_meta.h"
#include <iostream>
#include <fstream>



namespace radmat
{

  namespace 
  {
    // Write a key
    void write(ADATXML::XMLWriter& xml, 
        const std::string& path, 
        const Hadron::KeyHadronNPartNPtCorr_t::NPoint_t& param)
    {
      push(xml, path);

      write(xml, "t_slice", param.t_slice);
      write(xml, "Irrep", param.irrep);

      pop(xml);
    }
  } // namespace anonomyous


  tester test_xml_to_redstar(void)
  {
    tester m_test(__func__);

#if 0 // --- old junk 
    simpleWorld::ContinuumStatePrimitive csxml;

    csxml.J = 1;
    csxml.H = 1;
    csxml.parity = true;
    csxml.mom.resize(1); 
    csxml.mom.resize(3);
    csxml.mom[0] = 1;
    csxml.mom[1] = 0;
    csxml.mom[2] = -1;
    csxml.twoI_z = 2;
    csxml.name = "rho_rhoxD0_J0__J1";
    csxml.creation_op = true;
    csxml.smearedP = true;

    simpleWorld::ContinuumMatElem::State source;

    source.state = csxml;
    source.t_slice = 0;


    std::stringstream log;

    ContinuumProjector foobar(source);
    ContinuumProjector::const_iterator it;
    for(it = foobar.begin(); it != foobar.end(); ++it)
    {
      ADATXML::XMLBufferWriter xml; 
      write(xml,"Key",it->m_obj);
      log << SEMBLE::toScalar(it->m_coeff);
      xml.print(log);
    }


    log << "----------------------------------------------------\n";
    log << "----------------------------------------------------\n";
    log << "----------------------------------------------------\n";
    log << "----------------------------------------------------\n";


    simpleWorld::ContinuumMatElem matrix_element;
    matrix_element.source = source; 
    matrix_element.sink = source; 
    matrix_element.insertion.resize(1);
    matrix_element.insertion[0] = source;


    ContinuumThreePointProjector barbar(matrix_element);
    ContinuumThreePointProjector::const_iterator itt; 

    for(itt = barbar.begin(); itt != barbar.end(); ++itt)
    {
      ADATXML::XMLBufferWriter xml;
      write(xml,"Key",itt->m_obj);
      log << SEMBLE::toScalar(itt->m_coeff);
      xml.print(log);
    }


    std::ofstream out("ContinuumProjector_output.test.out");
    out << log.str();
    out.close();

#endif
    return m_test;
  }   


} // namespace ramat


