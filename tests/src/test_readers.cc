/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name :

 * Purpose :

 * Creation Date : 15-08-2012

 * Last Modified : Mon Dec 10 09:17:59 2012

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/

#include "../headers/tester.h"
#include "../headers/test_common.h"
#include <stdio.h>
#include <fstream>
#include <string>
#include <vector>
#include <iostream>
#include "radmat/load_data/simple_world.h"
#include "radmat/load_data/build_correlators.h"
#include "io/adat_xmlio.h"

namespace radmat
{

  //void make_test_file(const std::string &filename);

  tester test_readers(void)
  {
    tester m_test(__func__);

    simpleWorld::ContinuumStateXML csxml;

    csxml.J = 0;
    csxml.parity = true;
    csxml.mom.resize(1); 
    csxml.mom[0].resize(3);
    csxml.mom[0][0] = 1;
    csxml.mom[0][1] = 0;
    csxml.mom[0][2] = -1;
    csxml.twoI_z = 2;
    csxml.op_stem = "pion_pionxD0_J0__J0_";
    csxml.creation_op = true;
    csxml.smearedP = true;


    simpleWorld::ContinuumInsertionXML::Insertion csxml2; 
    simpleWorld::ContinuumInsertionXML insertion; 

    csxml2.J = 0;
    csxml2.parity = true;
    csxml2.op_stem = "b_b0xD0_J0__J0_";
    csxml2.twoI_z = 0;
    csxml2.creation_op = true;
    csxml2.smearedP = true;

    insertion.t_slice = -3; 
    insertion.time = csxml2;
    csxml2.op_stem = "rho_rhoxD0_J0__J1_";
    csxml2.J = 1; 
    insertion.space = csxml2;


    simpleWorld::ContinuumMatElemXML CLME, CLME_in;
    simpleWorld::ContinuumMatElemXML::State source,sink,ins;
    
    source.state = csxml;
    source.t_slice = 0;

    CLME.source = source;
    CLME.sink = source;
    CLME.insertion = insertion; 

    ADATXML::XMLBufferWriter xml;
    simpleWorld::write(xml,"Key",CLME);
    
   // xml.print(std::cout);
    std::string xml_name("simpleWorldXMLout.xml");
    std::ofstream out(xml_name.c_str());
    xml.print(out);
    out.close();

    ADATXML::XMLReader xml_in(xml_name);
    ADATXML::XMLBufferWriter xml2;
    simpleWorld::read(xml_in,"/Key",CLME_in);
    simpleWorld::write(xml2,"Key",CLME_in);
    std::string xml_name2 = xml_name + std::string("_2");
    std::ofstream out2(xml_name2.c_str());
    xml2.print(out2);
    out2.close();


    std::vector<simpleWorld::ContinuumMatElem> lots_of_stuff = getContinuumMatElemFromXML(CLME);
    std::vector<simpleWorld::ContinuumMatElem>::const_iterator print_lots_of_stuff;

    for(print_lots_of_stuff = lots_of_stuff.begin(); print_lots_of_stuff != lots_of_stuff.end(); ++print_lots_of_stuff)
      std::cout << *print_lots_of_stuff << std::endl;

    // user should run diff on two output files if they care
    TESTER_TEST(m_test,true,"foobar");


    ThreePointCorrIni_t empty_ini; 
    ADATXML::XMLBufferWriter empty_writer; 
    write(empty_writer,"inikeys",empty_ini);
    std::string empty_name("radmat.inikeys.xml");
    std::ofstream empty_out(empty_name.c_str());
    empty_writer.print(empty_out);
    empty_out.close(); 
    return m_test;
  }


} // namespace radmat


