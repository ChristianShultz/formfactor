/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : llsq_chisq.cc

 * Purpose :

 * Creation Date : 27-02-2013

 * Last Modified : Fri 06 Dec 2013 04:12:48 PM EST

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/


#include "llsq_chisq.h"
#include "semble/semble_semble.h"
#include "itpp/itbase.h"
#include <utility>

namespace radmat
{

  namespace
  {

    // NB: this is the covariance on the mean hence the extra (n-1)
    itpp::Mat<std::complex<double> > makeCov(const SEMBLE::SembleVector<std::complex<double> > &v)
    {
      const int dim = v.getN();
      const int ncfg = v.getB();
      itpp::Mat<std::complex<double> > out(dim,dim);
      itpp::Vec<std::complex<double> > work, mean; 
      mean = v.mean(); 
      out.zeros();

      for(int cfg = 0; cfg < ncfg; ++cfg)
      {
        work = v[cfg] - mean;

        for(int row = 0; row < dim; ++row)
          for(int col = row; col < dim; ++col) // hermitian matrix..
          {
            std::complex<double> val = std::conj(work[row])*work[col];
            out(row,col) += val;
            if(row != col)
              out(col,row) += std::conj(val);
          }
      }

      out /= double(ncfg*(ncfg - 1)); 

      //      std::cout << __func__ << std::endl;
      //      std::cout << out << std::endl;



      return out; 
    }

    std::pair<int,itpp::Mat<std::complex<double> > >
      invert_cmplx(const itpp::Vec<double> &in, const double tol)
      {
        const int sz = in.size();
        int nrset = sz;  
        itpp::Mat<std::complex<double> > out(sz,sz);
        out.zeros(); 

        for(int elem = 0; elem < sz; ++elem)
          if(in[elem] > tol)
          {
            --nrset; 
            out(elem,elem) = std::complex<double>(1./in[elem],0.); 
          }
          else
            break; // they come out ordered high to low and its already zero

        return std::make_pair(nrset,out); 
      }


    std::pair<int,itpp::Mat<std::complex<double> > > 
      makeInvCov(const SEMBLE::SembleVector<std::complex<double> > &v, double tol)
      {
        itpp::Mat<std::complex<double> > U,V;
        itpp::Vec<double> s;

        if(!!!itpp::svd(makeCov(v),U,s,V))
        {
          std::cerr << __func__ << ": error performing SVD" << std::endl;
          exit(1);
        }

        std::pair<int,itpp::Mat<std::complex<double> > >  inv; 

        inv = invert_cmplx(s,tol) ; 

        return std::make_pair(inv.first, V * inv.second * itpp::hermitian_transpose(U));
      }


  }




  std::pair<double,int> ChisqAndDoF(const SEMBLE::SembleVector<std::complex<double> > &data,
      const SEMBLE::SembleVector<std::complex<double> > &thy, double tol)
  {
    std::pair<int,itpp::Mat<std::complex<double> > > foofoo = makeInvCov(data,tol); 
    itpp::Mat<std::complex<double> > Cinv = foofoo.second;
    const int DoF = data.getN() - foofoo.first;  
    SEMBLE::SembleVector<std::complex<double> > v;
    std::complex<double> val(0.,0.);

    v = thy - data;
    itpp::Vec<std::complex<double> > meanv = v.mean(); 


    val = itpp::conj( meanv ) * ( Cinv * meanv ); 

    if(fabs(val.imag())/double(DoF) > tol)
    {
      std::cerr << __func__ << ": warning chisq with large imag component!!!" << std::endl;
      std::cerr << "chisq = " << val << " DoF = " << DoF << std::endl;
    }

    return std::make_pair(val.real(),DoF); 
  }




}
