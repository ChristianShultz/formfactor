/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : redstar_three_point_xml_interface.cc

 * Purpose :

 * Creation Date : 20-03-2014

 * Last Modified : Wed 26 Mar 2014 01:25:09 PM EDT

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/


#include "redstar_three_point_xml_interface.h"
#include "redstar_abstract_xml_factory.h"
#include "radmat/utils/printer.h"


namespace radmat
{

  struct abs_obj_read_printer
  {
    static void print(const std::string &msg)
    { std::cout << "abs_obj_read " << msg << std::endl;}
  };

  namespace
  {
    template<typename T>
      void doXMLRead(ADATXML::XMLReader &ptop, 
          const std::string &path, 
          T &place, const char * f)
      {
        if(ptop.count(path) > 0)
          read(ptop,path,place);
        else
        {
          std::cerr << __PRETTY_FUNCTION__ << ": Error, called by " 
            << f << " trying to read path, " << path
            << ", path was empty, exiting" << std::endl;
          exit(1);
        }
      }

  } // anonomyous 


  void read(ADATXML::XMLReader &xml, 
      const std::string &path, 
      AbstractNamedObject<AbsRedstarXMLInterface_t> &obj)
  {
    ADATXML::XMLReader ptop(xml,path); 
    if(ptop.count("object_name") > 0)
      read(ptop,"object_name",obj.object_name);
    else
      throw std::string( "AbsXML error rad" ); 

    try
    {
      obj.param = TheRedstarAbstractXMLFactoryEnv::callFactory(obj.object_name); 
      printer_function<abs_obj_read_printer>(obj.object_name); 
      obj.param->read(ptop,std::string("param"));
    }
    catch(std::exception &e) 
    {
      std::cout << __PRETTY_FUNCTION__ 
        << ": error, e.what() = " << e.what() << std::endl; 
      throw e; 
    }
    catch(std::string &s)
    {
      std::cout << __PRETTY_FUNCTION__ 
        << ": error, " << s << std::endl; 
    }
    catch(...)
    {
      std::cout << __PRETTY_FUNCTION__ 
        << ": some non standard error" << std::endl; 
      throw std::string("in") + std::string(__PRETTY_FUNCTION__); 
    }
  }




  void 
    RedstarThreePointXML::read(ADATXML::XMLReader &xml, 
        const std::string &path)
    { 
      ADATXML::XMLReader ptop(xml,path); 
      doXMLRead(ptop,"npt",npt,__PRETTY_FUNCTION__); 
      doXMLRead(ptop,"ensemble",ensemble,__PRETTY_FUNCTION__); 

      if( npt.size() != 3 )
      {
        printer_function<console_print>("incorrect npt size"); 
        exit(1); 
      }
    }

  std::string 
    RedstarThreePointXML::write() const
    {
      return npt[0].param->write() 
        + npt[1].param->write()
        + npt[2].param->write() 
        + "\n ensembel = " + ensemble; 
    }


  namespace TheRedstarThreePointAbstractXMLFactoryEnv
  {
    rHandle<RedstarThreePointXML> callFactory(const std::string &id)
    {
      rHandle<AbsRedstarXMLInterface_t> base_h; 
      base_h = TheRedstarAbstractXMLFactoryEnv::callFactory(id); 

      // this will bomb if we can't do the cast
      RedstarThreePointXML* derived;
      derived = dynamic_cast<RedstarThreePointXML*>(base_h.get_ptr());

      // allocate a new guy for the ouput handle since the handle we 
      // called from TheRedstarAbstractXMLFactoryEnv will delete whatever
      // derived was pointing at on exit of this scope   
      return rHandle<RedstarThreePointXML>( new RedstarThreePointXML(*derived) );  
    }
  }


} // radmat 

