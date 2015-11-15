/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

* File Name : data_representation_three_point.cc

* Purpose :

* Creation Date : 24-03-2014

* Last Modified : Tue 25 Mar 2014 02:25:42 PM EDT

* Created By : shultz

_._._._._._._._._._._._._._._._._._._._._.*/


#include "data_representation_three_point.h"


namespace radmat
{
  void write(ADATIO::BinaryWriter &bin, const DataRep3pt &d)
  {
    ADATIO::writeDesc(bin,d.l);
    ADATIO::writeDesc(bin,d.g);
    ADATIO::writeDesc(bin,d.r);
  }

  void read(ADATIO::BinaryReader &bin, DataRep3pt &d)
  {
    ADATIO::readDesc(bin,d.l); 
    ADATIO::readDesc(bin,d.g); 
    ADATIO::readDesc(bin,d.r); 
  }

}
