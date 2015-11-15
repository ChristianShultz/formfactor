#ifndef DATA_TAG_THREE_POINT_H
#define DATA_TAG_THREE_POINT_H 


#include "data_representation_three_point.h"
#include "data_tag_primitive.h"
#include "ensem/ensem.h"
#include "io/adat_xmlio.h"
#include "radmat/utils/stringify.h"



namespace radmat
{
  
  struct ThreePointDataTag; 
  REGISTER_STRINGIFY_TYPE(ThreePointDataTag); 
  

  struct ThreePointDataTag
    : public DataTagPrimitive
  {
    ThreePointDataTag();
  
    virtual ~ThreePointDataTag() {}
    virtual std::string type() const { return Stringify<ThreePointDataTag>(); }

    ENSEM::EnsemReal Q2() const; 
    void set_qsq_label(const double &q) {qsq_label = q;}
    double get_qsq_label() const { return qsq_label; }
    std::string mom_string(const ADATXML::Array<int> &p) const
    { std::stringstream ss; ss << p[0] << p[1] << p[2]; return ss.str();}
    std::string mom_string() const 
    { 
      return "left = " + mom_string(left_mom) 
      + " right = " + mom_string(right_mom);
    }

    std::string rot_qsq_tag(const bool b) const
    {
      return b ? rot_qsq_tag(origin_rep) : rot_qsq_tag(data_rep); 
    }

    std::string full_irrep_id( const DataRep3pt &d, const std::string &s) const
    {
      std::stringstream ss; 

      if( d.is_cubic( s ) )
      {
        rHandle<CubicRep> foo = d.call_rep_cub(s); 
        ss <<  foo->rep_id(); 
      }
      else
      {
        rHandle<LorentzRep> foo = d.call_rep_lor(s); 
        ss <<  foo->rep_id(); 
      }

      return ss.str(); 
    }

    std::string irrep_id( const DataRep3pt &d, const std::string &s) const
    {
      std::stringstream ss; 

      if( d.is_cubic( s ) )
      {
        rHandle<CubicRep> foo = d.call_rep_cub(s); 
        ss <<  foo->rep_irrep(); 
      }
      else
      {
        rHandle<LorentzRep> foo = d.call_rep_lor(s); 
        ss <<  foo->rep_id(); 
      }

      return ss.str(); 
    }

    std::string rot_qsq_tag( const DataRep3pt &d) const
    {
      std::stringstream ss; 
      ss << "l" +  irrep_id(d,d.l) + "xr" + irrep_id(d,d.r) + "xg" + irrep_id(d,d.g);
      return ss.str();  
    }


    // specify some representation information  
    DataRep3pt origin_rep; 
    DataRep3pt data_rep; 

    // unique file id 
    std::string file_id; 

    // matrix decomp key
    std::string mat_elem_id;

    // the invariants
    int left_row, gamma_row, right_row; 
    ADATXML::Array<int> left_mom,q,right_mom; 
    ENSEM::EnsemReal left_E, right_E; 

    // sorting label 
    double qsq_label; 
    double mom_fac; 
  }; 

  // binary read write for serialization 
  void write(ADATIO::BinaryWriter &, const ThreePointDataTag &);
  void read(ADATIO::BinaryReader &, ThreePointDataTag &); 

} 



#endif /* DATA_TAG_THREE_POINT_H */
