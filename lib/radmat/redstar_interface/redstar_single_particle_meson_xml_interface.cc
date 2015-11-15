/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : redstar_single_particle_meson_xml_interface.cc

 * Purpose :

 * Creation Date : 20-03-2014

 * Last Modified : Thu 24 Apr 2014 06:31:30 PM EDT

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/


#include "redstar_single_particle_meson_xml_interface.h"
#include "radmat/utils/pow2assert.h"
#include "radmat/data_representation/data_representation_lorentz_groups.h"
#include "formfac/formfac_qsq.h"
#include "hadron/irrep_util.h"
#include <exception>
#include <iostream>
#include <sstream>


namespace radmat
{

  namespace
  {

    template<typename T>
      void doXMLRead(ADATXML::XMLReader &ptop, const std::string &path, T &place, const char * f)
      {
        if(ptop.count(path) > 0)
          read(ptop,path,place);
        else
        {
          std::cerr << __PRETTY_FUNCTION__ << ": Error, called by " << f << " trying to read path, " << path
            << ", path was empty, exiting" << std::endl;
          exit(1);
        }
      }


    // stringy mom
    std::string canonical_stringmom(const ADATXML::Array<int> &inp)
    {
      ADATXML::Array<int> can = FF::canonicalOrder(inp);
      std::stringstream ss;
      ss << can[0] << can[1] << can[2];
      return ss.str(); 
    }


    // get the star
    ADATXML::Array<ADATXML::Array<int> > star_p(const ADATXML::Array<ADATXML::Array<int> > &pin)
    {

      std::list<ADATXML::Array<int> > work;
      std::map<std::string,std::list<ADATXML::Array<int> > > seen;
      std::map<std::string,std::list<ADATXML::Array<int> > >::const_iterator it;
      std::list<ADATXML::Array<int> >::const_iterator listit; 

      int nele(0); 

      const unsigned int sz = pin.size(); 

      for(unsigned int elem = 0; elem < sz; ++elem)
      {
        std::string p = canonical_stringmom(pin[elem]);
        it = seen.find(p);
        if (it == seen.end())
        {
          work = Hadron::generateLittleGroupMom(Hadron::generateLittleGroup(pin[elem]),
              FF::canonicalOrder(pin[elem]));
          nele += work.size();
          seen.insert(std::pair<std::string,std::list<ADATXML::Array<int> > >(p,work)); 
        }        
      }


      ADATXML::Array<ADATXML::Array<int> > ret(nele);
      nele = 0;

      for(it = seen.begin(); it != seen.end(); ++it)
        for(listit = it->second.begin(); listit != it->second.end(); ++listit)
        {
          ret[nele] = *listit; 
          ++nele; 
        }

      return ret; 
    }


    std::string doPrint(const ADATXML::Array<int> &t)
    {
      std::stringstream ss; 
      for(int i = 0; i < t.size(); ++i)
        ss << " " << t[i];
      return ss.str();  
    } 

    std::string doPrint(const ADATXML::Array< ADATXML::Array<int> > &t)
    {
      std::stringstream ss; 
      for(int i = 0; i < t.size(); ++i)
        ss << doPrint( t[i] ) << std::endl;
      return ss.str();  
    } 



  } // anonomyous 



  std::string 
    RedstarSingleParticleMesonXML::write(void) const
    {
      std::stringstream ss; 
      ss << cont_rep << " H= " << doPrint(H) 
        << " fill_star= " << fill_star 
        << " twoI_z= " << twoI_z << " name = " << name 
        << " creation_op= " << creation_op << " smearedP= " << smearedP
        << " isProjected= " << isProjected << " t_slice= " << t_slice
        << " momlist: \n " << doPrint(mom); 
      return ss.str(); 
    }

  void 
    RedstarSingleParticleMesonXML::read(ADATXML::XMLReader &xml, const std::string &path)
    {
      ADATXML::XMLReader ptop(xml,path);
      doXMLRead(ptop,"cont_rep",cont_rep,__PRETTY_FUNCTION__);
      doXMLRead(ptop,"H",H,__PRETTY_FUNCTION__);   
      doXMLRead(ptop,"fill_star",fill_star,__PRETTY_FUNCTION__);   
      doXMLRead(ptop,"mom",mom,__PRETTY_FUNCTION__);   
      doXMLRead(ptop,"twoI_z",twoI_z,__PRETTY_FUNCTION__);   
      doXMLRead(ptop,"name",name,__PRETTY_FUNCTION__);   
      doXMLRead(ptop,"creation_op",creation_op,__PRETTY_FUNCTION__);   
      doXMLRead(ptop,"smearedP",smearedP,__PRETTY_FUNCTION__);   
      doXMLRead(ptop,"isProjected",isProjected,__PRETTY_FUNCTION__);   
      doXMLRead(ptop,"t_slice",t_slice,__PRETTY_FUNCTION__);   


      rHandle<LorentzRep> cont_rep_h;
      cont_rep_h = LorentzRepresentationFactoryEnv::callFactory(cont_rep); 
      int spin = cont_rep_h->rep_spin(); 

      if ( H.size() == 0) 
      {
        H.resize(2*spin + 1);
        for(int h = -spin; h <= spin; ++h)
          H[spin-h] = h; 
      }

      if(fill_star)
        mom = star_p(mom); 
    }

  //  //! xml writer
  //  void 
  //    RedstarSingleParticleMesonXML::write(ADATXML::XMLWriter &xml, const std::string &path ) const
  //    {
  //      ADATXML::push(xml,path);
  //      ADATXML::write(xml,"cont_rep",cont_rep);
  //      ADATXML::write(xml,"H",H);
  //      ADATXML::write(xml,"parity",parity);
  //      ADATXML::write(xml,"mom",mom);
  //      ADATXML::write(xml,"twoI_z",twoI_z);
  //      ADATXML::write(xml,"name",name);
  //      ADATXML::write(xml,"creation_op",creation_op);
  //      ADATXML::write(xml,"smearedP",smearedP);
  //      ADATXML::write(xml,"isProjected",isProjected); 
  //      ADATXML::write(xml,"t_slice",t_slice);
  //      ADATXML::pop(xml);
  //    }

} // radmat




