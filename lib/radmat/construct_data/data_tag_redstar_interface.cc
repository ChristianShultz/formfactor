/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

* File Name : data_tag_redstar_interface.cc

* Purpose :

* Creation Date : 24-03-2014

* Last Modified : Fri 03 Oct 2014 01:54:46 PM EDT

* Created By : shultz

_._._._._._._._._._._._._._._._._._._._._.*/

#include "data_tag_redstar_interface.h"
#include "radmat/utils/utils.h"



namespace radmat
{

  namespace
  {

    std::string string_mom(const ADATXML::Array<int> &p)
    {
      std::stringstream ss; 
      ss << p[0] << p[1] << p[2]; 
      return ss.str(); 
    }
  
  std::string gen_file_id(const ThreePointDataTag &t)
  {
    std::stringstream ss;  
    ss << t.data_rep.l << "r" << t.left_row 
      << "_p" << string_mom(t.left_mom)
      << "." << t.data_rep.g << "r" << t.gamma_row 
      << "_p" << string_mom(t.q)
      << "." << t.data_rep.r << "r" << t.right_row 
      << "_p" << string_mom(t.right_mom); 
    return ss.str(); 
  }


  // energies get filled when we start playing with databases later 
  ThreePointDataTag
    generate_tag(const ThreePointData &d, 
        const double mom_factor, 
        const double m_lefty, 
        const double m_righty, 
        const std::string &elem_id)
    {
      ThreePointDataTag ret; 
      ret.origin_rep = d.origin_rep; 
      ret.data_rep = d.data_rep;
      ret.mat_elem_id = elem_id; 
      ret.left_row = d.left_row; 
      ret.gamma_row = d.gamma_row; 
      ret.right_row = d.right_row; 
      ret.left_mom = d.data.begin()->m_obj.npoint[1].irrep.mom;  
      ret.q = d.data.begin()->m_obj.npoint[2].irrep.mom;  
      ret.right_mom = d.data.begin()->m_obj.npoint[3].irrep.mom;  

      ret.mom_fac = mom_factor; 
      ret.qsq_label = Mink_qsq( ret.left_mom, m_lefty, 
          ret.right_mom, m_righty, mom_factor);

      ret.file_id = gen_file_id(ret); 


      return ret; 
    }

  ///
  TaggedEnsemRedstarNPtBlock
    generate_tagged_data( const ThreePointData &d, 
        const double mom_factor, 
        const double m_lefty, 
        const double m_righty, 
        const std::string &elem_id)
    {
      return TaggedEnsemRedstarNPtBlock( d.data , 
        generate_tag(d,mom_factor,m_lefty,m_righty,elem_id)); 
    }
   
  } // anonomyous 

  std::vector<TaggedEnsemRedstarNPtBlock> 
    tag_lattice_xml(const std::vector<ThreePointData> &d, 
        const double mom_factor, 
        const double m_lefty, 
        const double m_righty, 
        const std::string &elem_id)
    {
      std::vector<TaggedEnsemRedstarNPtBlock> ret; 
      ret.reserve(d.size()); 

      std::vector<ThreePointData>::const_iterator it; 
      for(it = d.begin(); it != d.end(); ++it)
        ret.push_back( generate_tagged_data( *it, 
              mom_factor, 
              m_lefty,
              m_righty, 
              elem_id) );

      return ret; 
    }

  std::string generate_file_id(const ThreePointDataTag &t)
  {
    return gen_file_id(t); 
  }

} // radmat 



