/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

* File Name : data_tag_three_point.cc

* Purpose :

* Creation Date : 24-03-2014

* Last Modified : Wed 20 Aug 2014 08:37:22 PM EDT

* Created By : shultz

_._._._._._._._._._._._._._._._._._._._._.*/


#include "data_tag_three_point.h"
#include "semble/semble_meta.h"


namespace radmat
{

  ThreePointDataTag::ThreePointDataTag()
  {
    qsq_label = -10000.;
    left_E.resize(1); 
    left_E = SEMBLE::toScalar(double(0.)); 
    right_E = left_E;  
  }


  ENSEM::EnsemReal 
    ThreePointDataTag::Q2() const
    {
      // use the energies to determine the size params of the 
      // ensem objects then zero them, then set them based on 
      // the value of the Real types 
      ENSEM::EnsemReal Q0, Q1, Q2, Q3, QQ; 
      ENSEM::Real zero,q1,q2,q3; 
      zero = ENSEM::toDouble(0.); 
      Q0 = (left_E - right_E) * (left_E - right_E); 
      Q1 = zero*Q0;
      Q2 = zero*Q0;
      Q3 = zero*Q0;

      q1 = ENSEM::toDouble( mom_fac*( left_mom[0] - right_mom[0] ) ); 
      q2 = ENSEM::toDouble( mom_fac*( left_mom[1] - right_mom[1] ) ); 
      q3 = ENSEM::toDouble( mom_fac*( left_mom[2] - right_mom[2] ) ); 

      Q1 = q1*q1; 
      Q2 = q2*q2; 
      Q3 = q3*q3; 

      // Q2 = - q \dot q 
      QQ = -Q0 + Q1 + Q2 + Q3;
      
    //  double qq = SEMBLE::toScalar(ENSEM::mean(QQ) ); 
    //  if( qq < 0 ) 
    //  {
    //    std::cout << __PRETTY_FUNCTION__ 
    //      << mom_string() << " is messed up right now??" 
    //      << " think q2 is " << qq << std::endl;  
    //    double a0,a1,a2,a3; 
    //    a0 = SEMBLE::toScalar(ENSEM::mean(Q0));
    //    a1 = SEMBLE::toScalar(ENSEM::mean(Q1));
    //    a2 = SEMBLE::toScalar(ENSEM::mean(Q2));
    //    a3 = SEMBLE::toScalar(ENSEM::mean(Q3));
    //    std::cout << " Q0 " << a0 << std::endl;
    //    std::cout << " Q1 " << a1 << std::endl;
    //    std::cout << " Q2 " << a2 << std::endl;
    //    std::cout << " Q3 " << a3 << std::endl;
    //  }

      return QQ; 
    }

  void write(ADATIO::BinaryWriter &bin, const ThreePointDataTag &d)
  {
    write(bin,d.origin_rep);
    write(bin,d.data_rep);
    ADATIO::writeDesc(bin,d.file_id);
    ADATIO::writeDesc(bin,d.mat_elem_id); 
    ADATIO::write(bin,d.left_row);
    ADATIO::write(bin,d.gamma_row);
    ADATIO::write(bin,d.right_row);
    ADATIO::write(bin,d.left_mom);
    ADATIO::write(bin,d.q); 
    ADATIO::write(bin,d.right_mom);
    ENSEM::write(bin,d.left_E);
    ENSEM::write(bin,d.right_E); 
    ADATIO::write(bin,d.qsq_label); 
    ADATIO::write(bin,d.mom_fac); 
  }


  void read(ADATIO::BinaryReader &bin, ThreePointDataTag &d)
  {
    read(bin,d.origin_rep);
    read(bin,d.data_rep);
    ADATIO::readDesc(bin,d.file_id);
    ADATIO::readDesc(bin,d.mat_elem_id); 
    ADATIO::read(bin,d.left_row);
    ADATIO::read(bin,d.gamma_row);
    ADATIO::read(bin,d.right_row);
    ADATIO::read(bin,d.left_mom);
    ADATIO::read(bin,d.q); 
    ADATIO::read(bin,d.right_mom);
    ENSEM::read(bin,d.left_E);
    ENSEM::read(bin,d.right_E); 
    ADATIO::read(bin,d.qsq_label); 
    ADATIO::read(bin,d.mom_fac); 
  }



} 
