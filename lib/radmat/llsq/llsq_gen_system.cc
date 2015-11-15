/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

* File Name : llsq_gen_system.cc

* Purpose :

* Creation Date : 25-01-2013

* Last Modified : Fri Jan 25 14:51:35 2013

* Created By : shultz

_._._._._._._._._._._._._._._._._._._._._.*/



#include "llsq_gen_system.h"
#include <algorithm>


namespace radmat
{
  bool isActive(const LLSQDataPoint &d)
  { 
    return !!!(!d.zero.first && !d.one.first && !d.two.first && !d.three.first);   // true if active so see if we can accumulate a false and then if we do it had active
  }


  struct LLSQDataPoint_comparator
  {
    bool operator()(const LLSQDataPoint &a, const LLSQDataPoint &b)
    {
      if (isActive(a))
        return true;
      return !!!isActive(b);
    }
  } theLLSQDataPointComparator; 


  std::vector<LLSQDataPoint> left_sort(const std::vector<LLSQDataPoint> &indata)
  {
    std::vector<LLSQDataPoint> outdata(indata); 

    std::sort(outdata.begin(),outdata.end(),theLLSQDataPointComparator);

    return outdata;
  }

}
