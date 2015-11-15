/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : radmat_overlap_key_val_db.cc

 * Purpose :

 * Creation Date : 02-10-2012

 * Last Modified : Sun 04 May 2014 04:20:22 PM EDT

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/


#include "radmat_overlap_key_val_db.h"
#include "ensem/ensem.h"
#include "hadron/ensem_filenames.h"

#include <sstream>

//#define DEBUG_RADMAT_KEY_VAL_DB

#ifdef DEBUG_RADMAT_KEY_VAL_DB

#include <iostream>
#include "semble/semble_meta.h"

using namespace SEMBLE;

#endif


namespace
{
  //! Binary output
  // --- this is copied from ensem_scalar_specific.h, line 636 -- in order to avoid changing ADAT 
  // I have just provided writer that works for our purposes. This is a hack and may break if the serialization 
  // routines in ADAT change at some point. 
  template<class T>  
    inline
    void stupid_writer(ADATIO::BinaryWriter& bin, const ENSEM::Ensem< ENSEM::OScalar<T> >& d)
    {
      d.checkSize(__func__);

      int num  = d.size();
      // original line is commented  -- method elem() doesn't exist in my version of ADAT..
      // int len  = d.elem().numElem();
      int len = d.numElem();
      int type = ENSEM::EnsbcIO< ENSEM::EScalar< ENSEM::OScalar<T> > >::type;

      ADATIO::write(bin, num);
      ADATIO::write(bin, len);
      ADATIO::write(bin, type);

      for(int n=0; n < num; ++n) 
        ENSEM::write(bin, d.elem(n));
    }

}



namespace radmat
{

  std::ostream& operator<<(std::ostream& os, const RadmatExtendedKeyHadronNPartIrrep_t &op)
  {
    os << "particle id - " << op.m_particle_id << " fn->  " 
      << Hadron::ensemFileName(op.m_basic_key) << "\n";
    os << op.m_basic_key;
    return os; 
  }

  void read(ADATIO::BinaryReader &bin, RadmatExtendedKeyHadronNPartIrrep_t &param)
  {
    readDesc(bin,param.m_particle_id);

#ifdef  DEBUG_RADMAT_KEY_VAL_DB
    std::cout << __func__ << " done reading particle id - " << param.m_particle_id << std::endl;
#endif
    ::Hadron::read(bin,param.m_basic_key);
#ifdef  DEBUG_RADMAT_KEY_VAL_DB
    std::cout << __func__ << " done reading HadronNPartIrrep_t - " << param.m_basic_key << std::endl;
#endif

  }


  void write(ADATIO::BinaryWriter &bin, const RadmatExtendedKeyHadronNPartIrrep_t &param)
  {
    writeDesc(bin,param.m_particle_id);
    ::Hadron::write(bin,param.m_basic_key);
  }


  //! read a key from xml
  void read(ADATXML::XMLReader &xml, const std::string &path,  RadmatExtendedKeyHadronNPartIrrep_t &param)
  {
    ADATXML::XMLReader ptop(xml,path);

    ADATXML::read(ptop, "particle_id" , param.m_particle_id);
    ::Hadron::read(ptop, "basic_key" , param.m_basic_key);
  }


  //! write a key to xml
  void write(ADATXML::XMLWriter &xml, const std::string &path, const RadmatExtendedKeyHadronNPartIrrep_t &param)
  {
    ADATXML::push(xml,path);
    ADATXML::write(xml,"particle_id",param.m_particle_id);
    ::Hadron::write(xml,"basic_key",param.m_basic_key);
    ADATXML::pop(xml);
  }

    //! a uniqe key
  std::string fileName(const RadmatExtendedKeyHadronNPartIrrep_t &param)
  {
    std::stringstream ss; 
    ss << param.m_particle_id << "." << Hadron::ensemFileName(param.m_basic_key);
    return ss.str();
  }




  void read(ADATIO::BinaryReader & bin, RadmatMassOverlapData_t &param)
  {
    ENSEM::read(bin,param.E());
    ENSEM::read(bin,param.Z());
  }


  void write(ADATIO::BinaryWriter &bin , const RadmatMassOverlapData_t &param)
  {
    stupid_writer(bin,param.E());
    stupid_writer(bin,param.Z());
  }

  bool operator==(const  RadmatMassOverlapData_t & a, const  RadmatMassOverlapData_t &b)
  {
    bool success = true;

    if( (a.E().size() != b.E().size()) | (a.Z().size() != b.Z().size()))
    {

#ifdef  DEBUG_RADMAT_KEY_VAL_DB
      std::cout << __func__ << " a.E().size() " << a.E().size() << " b.E().size() " << b.E().size() << std::endl;
      std::cout << __func__ << " a.Z().size() " << a.Z().size() << " b.Z().size() " << b.Z().size() << std::endl;
      std::cout << __func__ << ": returning false" << std::endl;
#endif


      return false;
    }

    ENSEM::BinaryReturn<RadmatMassOverlapData_t::DATA, RadmatMassOverlapData_t::DATA,ENSEM::OpEQ>::Type_t E,Z;

    E = (a.E() == b.E());
    Z = (a.Z() == b.Z());

    for(int i = 0; i < a.E().size(); ++i)
    {
#ifdef  DEBUG_RADMAT_KEY_VAL_DB

      if(!!!(ENSEM::toBool(E.elem(i) ) && ENSEM::toBool(Z.elem(i))))
      {
        if(!!!ENSEM::toBool(E.elem(i)))
        {
          std::cout << __func__ << " E elem "  << i << " failed equality \n" << std::endl;
          std::cout << "a " << toScalar(a.E().elem(i)) << "   b " << toScalar(b.E().elem(i)) << std::endl;   
        }

        if(!!!ENSEM::toBool(Z.elem(i)))
        {
          std::cout << __func__ << " Z elem "  << i << " failed equality \n" << std::endl;
          std::cout << "a " << toScalar(a.Z().elem(i)) << "   b " << toScalar(b.Z().elem(i)) << std::endl;  
        }
        return false;
      }
#endif

      success &= ENSEM::toBool(E.elem(i));
      success &= ENSEM::toBool(Z.elem(i));
    }

    return success; 
  }


}


#undef  DEBUG_RADMAT_KEY_VAL_DB

