/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : redstar_vector_current_xml_interface.cc

 * Purpose :

 * Creation Date : 20-03-2014

 * Last Modified : Wed 26 Mar 2014 11:35:04 AM EDT

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/

#include "redstar_vector_current_xml_interface.h"
#include "redstar_photon_props.h"


namespace radmat
{

  namespace
  {

    template<typename T>
      void doXMLRead(ADATXML::XMLReader &ptop, 
          const std::string &path, 
          T &place, 
          const char * f)
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

    std::string toString(const RedstarVectorCurrentXML::insertion &i)
    {
      std::stringstream ss; 
      ss << "active= " << i.active << " create= " << PHOTON_CREATE
        << " smear= " << i.smearedP << " photons: ";
      for(int j = 0; j < i.photons.size(); ++j)
        ss << "(" << i.photons[j].coeff_r 
          << " + " << i.photons[j].coeff_i 
          << "i) x" << i.photons[j].name << "   ";
      return ss.str(); 
    }


  } // anonomyous  




  void 
    RedstarVectorCurrentXML::read(ADATXML::XMLReader &xml, const std::string &path)
    {
      ADATXML::XMLReader ptop(xml,path); 
      doXMLRead(ptop,"pmin",pmin,__PRETTY_FUNCTION__);       
      doXMLRead(ptop,"pmax",pmax,__PRETTY_FUNCTION__);       
      doXMLRead(ptop,"t_slice",t_slice,__PRETTY_FUNCTION__);       
      doXMLRead(ptop,"time",time,__PRETTY_FUNCTION__);       
      doXMLRead(ptop,"space",space,__PRETTY_FUNCTION__);       


      if ( time.active && space.active ) 
      {
        std::cout << __PRETTY_FUNCTION__ << __FILE__ << __LINE__
          << ": Warning, you decided to mix time and" 
          << " space and I don't think you should " << std::endl; 
      }
    }

  std::string 
    RedstarVectorCurrentXML::write(void) const
    {
      std::stringstream ss;
      ss << "pmin= " << pmin << " pmax=" << pmax << " t_slice= " << t_slice; 
      ss << "\ntime:\n" << toString(time) << std::endl;
      ss << "\nspace:\n" << toString(space) << std::endl;
      return ss.str(); 
    }


    void write(ADATXML::XMLWriter &xml, 
        const std::string &path,
        const RedstarVectorCurrentXML::pfrag &p)
    {
      ADATXML::push(xml,path);
      ADATXML::write(xml,"coeff_r",p.coeff_r);
      ADATXML::write(xml,"coeff_i",p.coeff_i); 
      ADATXML::write(xml,"name",p.name); 
      ADATXML::pop(xml);
    }

    void write(ADATXML::XMLWriter &xml,
        const std::string &path, 
        const RedstarVectorCurrentXML::insertion &i)
    {
      ADATXML::push(xml,path);
      ADATXML::write(xml,"active",i.active);
      ADATXML::write(xml,"smearedP",i.smearedP); 
      write(xml,"photons",i.photons); 
      ADATXML::pop(xml);
    }

    void read(ADATXML::XMLReader &xml, 
        const std::string &path,
        RedstarVectorCurrentXML::pfrag &p)  
    {
      ADATXML::XMLReader ptop(xml,path); 
      doXMLRead(ptop,"coeff_r",p.coeff_r,__PRETTY_FUNCTION__); 
      doXMLRead(ptop,"coeff_i",p.coeff_i,__PRETTY_FUNCTION__); 
      doXMLRead(ptop,"name",p.name,__PRETTY_FUNCTION__); 
    }

    void read(ADATXML::XMLReader &xml, 
        const std::string &path,
        RedstarVectorCurrentXML::insertion &i)  
    {
      ADATXML::XMLReader ptop(xml,path); 
      doXMLRead(ptop,"active",i.active,__PRETTY_FUNCTION__); 
      doXMLRead(ptop,"smearedP",i.smearedP,__PRETTY_FUNCTION__); 
      doXMLRead(ptop,"photons",i.photons,__PRETTY_FUNCTION__); 
    }


} // radmat
