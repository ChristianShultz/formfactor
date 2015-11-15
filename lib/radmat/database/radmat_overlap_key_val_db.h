#ifndef RADMAT_OVERLAP_KEY_VAL_DB_H
#define RADMAT_OVERLAP_KEY_VAL_DB_H


#include "hadron/hadron_npart_irrep.h"
#include "io/adat_io.h"
#include "io/adat_xmlio.h"
#include "ensem/ensem.h"
#include <string>
#include <iostream>


namespace radmat
{

  //! an extended key specifying the name of the particle to avoid ambiguity involved in state numbering
  struct RadmatExtendedKeyHadronNPartIrrep_t
  {

    RadmatExtendedKeyHadronNPartIrrep_t(void)
      : m_particle_id(std::string("default"))
    {}


    RadmatExtendedKeyHadronNPartIrrep_t(const std::string particle_id, const Hadron::KeyHadronNPartIrrep_t basic_key)
      : m_particle_id(particle_id) , m_basic_key(basic_key) 
    {}

    void dagger(void)
    {
      m_basic_key.creation_op = !!!m_basic_key.creation_op; 
    }

    void doLG_symmetry(void)
    {
      m_basic_key.mom[0] = 0;
      m_basic_key.mom[1] = 0;
      m_basic_key.mom[2] = 0;
      m_basic_key.row = 0; 
    }

    std::string m_particle_id;    
    Hadron::KeyHadronNPartIrrep_t m_basic_key;
  };

  //! for error messages
  std::ostream& operator<<(std::ostream& , const RadmatExtendedKeyHadronNPartIrrep_t &);

  //! read a serialized key
  void read(ADATIO::BinaryReader &bin, RadmatExtendedKeyHadronNPartIrrep_t &param);

  //! write a key to binary
  void write(ADATIO::BinaryWriter &bin, const RadmatExtendedKeyHadronNPartIrrep_t &param);

  //! read a key from xml
  void read(ADATXML::XMLReader &xml, const std::string &path,  RadmatExtendedKeyHadronNPartIrrep_t &param);

  //! write a key to xml
  void write(ADATXML::XMLWriter &xml, const std::string &path, const RadmatExtendedKeyHadronNPartIrrep_t &param);

  //! a uniqe key
  std::string fileName(const RadmatExtendedKeyHadronNPartIrrep_t &param);





  //! a serializable data container holding the energy and overlap
  struct RadmatMassOverlapData_t
  {
    typedef ENSEM::EnsemReal DATA;

    RadmatMassOverlapData_t(void)
    {}

    RadmatMassOverlapData_t(DATA &E, DATA &Z)
      : m_E(E) , m_Z(Z) 
    { }

    const DATA& E(void) const {return m_E;}
    DATA& E(void) {return m_E;}

    const DATA& Z(void) const {return m_Z;}
    DATA& Z(void) {return m_Z;}

    DATA m_E;
    DATA m_Z;
  };

  //! read serialized data
  void read(ADATIO::BinaryReader &bin, RadmatMassOverlapData_t &param);

  //! write serialized data
  void write(ADATIO::BinaryWriter &bin, const RadmatMassOverlapData_t &param);

  //! for testing routines
  bool operator==(const  RadmatMassOverlapData_t & , const  RadmatMassOverlapData_t &);


} // namespace radmat

#endif /* RADMAT_OVERLAP_KEY_VAL_DB_H */
