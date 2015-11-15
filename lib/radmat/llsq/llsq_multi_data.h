#ifndef LLSQ_MULTI_TIME_DATA_H
#define LLSQ_MULTI_TIME_DATA_H


#include "semble/semble_semble.h"
#include "ensem/ensem.h"
#include <vector>
#include <complex>




namespace radmat
{



  // another instance of caveat emptor 

  template<typename TAG_TYPE, typename SEMBLE_TYPE>
    struct LLSQMultiData
    {
      typedef TAG_TYPE TT; 
      typedef SEMBLE_TYPE ST; 


      LLSQMultiData(void) : first_append(true) {}
      LLSQMultiData(const SEMBLE::SembleMatrix<ST> &, const int, const TT &);
      LLSQMultiData(const LLSQMultiData<TT,ST> &);


      LLSQMultiData<TT,ST>& operator=(const LLSQMultiData<TT,ST> &); 

      SEMBLE::SembleMatrix<ST>& data(void) {return m_data;}
      const SEMBLE::SembleMatrix<ST>& data(void) const {return m_data;}
      std::vector<TT> tags(void) {return m_tags;}
      const std::vector<TT>& tags(void) const {return m_tags;}

      void append_row_semble(const SEMBLE::SembleVector<ST> &, const TT &);
      SEMBLE::SembleVector<ST> get_col_semble(const int col) const;
      SEMBLE::SembleVector<ST> get_row_semble(const int row) const;

      void append_row_ensem(const typename SEMBLE::PromoteEnsemVec<ST>::Type &, const TT &);
      typename SEMBLE::PromoteEnsemVec<ST>::Type get_col_ensem(const int col) const; 
      typename SEMBLE::PromoteEnsemVec<ST>::Type get_row_ensem(const int row) const;

      int nrows(void) const {return m_data.getN();}
      int ncols(void) const {return m_data.getM();}
      bool empty(void) const {return ( m_data.getN() == 0 );}

      const TT& get_tag(const int idx) const; 
      TT& get_tag(const int idx); 
      void set_tag(const int idx, const TT &); 
      void splash_tags(void) const; 

      bool first_append; 
      SEMBLE::SembleMatrix<ST> m_data; 
      std::vector<TT> m_tags; 
    };



  template<typename TT, typename ST>
    LLSQMultiData<TT,ST>::LLSQMultiData(const typename SEMBLE::SembleMatrix<ST> &s,
        const int n , const TT &tt)
    :  first_append(false), m_data(s) 
    {
      m_tags.resize(n,tt);
    }

  template<typename TT, typename ST>
    LLSQMultiData<TT,ST>::LLSQMultiData(const LLSQMultiData<TT,ST> &o)
    : first_append(o.first_append) , m_data(o.m_data) , m_tags(o.m_tags)  {}

  template<typename TT, typename ST>
    LLSQMultiData<TT,ST>& LLSQMultiData<TT,ST>::operator=(const LLSQMultiData<TT,ST> &o)
    {
      if(this != &o)
      {
        first_append = o.first_append; 
        m_data = o.m_data;
        m_tags.clear();
        m_tags = o.m_tags; 
      }
      return *this;
    }

  template<typename TT, typename ST>
    void LLSQMultiData<TT,ST>::append_row_semble(const SEMBLE::SembleVector<ST> &s, const TT &tt)
    {
      if(first_append)  
      {

        SEMBLE::SembleMatrix<ST> foofoo(s.getB(),1,s.getN());
        for(int mm = 0; mm < s.getN(); ++mm)
          foofoo.loadEnsemElement(0,mm, s.getEnsemElement(mm)); 

        m_data = foofoo; 

        if(!!!m_tags.empty())
        {
          std::cerr << __func__ << ": error appending rows" << std::endl;
          exit(1);
        }
        /*
           std::cout << "s" << std::endl;
           std::cout << s.mean() << std::endl;
           std::cout << s.getB() << "x" << s.getN() << std::endl;
           std::cout << "first_append" << std::endl;
           std::cout << m_data.getB() << "x" << m_data.getN() << "x" << m_data.getM() << std::endl;
           std::cout << "foofoo" << std::endl;
           std::cout << foofoo.getB() << "x" << foofoo.getN() << "x" << foofoo.getM() << std::endl;
         */
        m_tags.push_back(tt);
        first_append = false;  
      }
      else
      {
        /*
           std::cout << "s" << std::endl;
           std::cout << s.getB() << "x" << s.getN() << std::endl;
           std::cout << m_data.getB() << "x" << m_data.getN() << "x" << m_data.getM() << std::endl;
         */
        m_data.append_row(s);
        m_tags.push_back(tt); 
      }
    }

  template<typename TT, typename ST>
    SEMBLE::SembleVector<ST> LLSQMultiData<TT,ST>::get_col_semble(const int col) const
    {
      if(col > m_data.getM())
      {
        std::cerr << __func__ << ": error out of bounds" << std::endl;
        exit(1);
      }

      return m_data.getCol(col); 
    }


  template<typename TT, typename ST>
    SEMBLE::SembleVector<ST> LLSQMultiData<TT,ST>::get_row_semble(const int row) const
    {
      if(row > m_data.getN())
      {
        std::cerr << __func__ << ": error out of bounds" << std::endl;
        exit(1); 
      }

      return m_data.getRow(row); 
    }

  template<typename TT, typename ST>
    void  LLSQMultiData<TT,ST>::append_row_ensem(const typename SEMBLE::PromoteEnsemVec<ST>::Type &s, const TT &tt)
    {
      SEMBLE::SembleVector<ST> foo;
      foo = s; 
      this->append_row_semble(foo,tt); 
    }



  template<typename TT, typename ST>
    typename SEMBLE::PromoteEnsemVec<ST>::Type  LLSQMultiData<TT,ST>::get_col_ensem(const int col) const
    {
      SEMBLE::SembleVector<ST> foo;
      foo = this->get_col_semble(col); 
      typename SEMBLE::PromoteEnsemVec<ST>::Type ret; 
      ret.resize(foo.getB());
      ret.rsizeObs(foo.getN());
      for(int nn = 0; nn < foo.getN(); ++nn)
        ENSEM::pokeObs(ret,foo.getEnsemElement(nn),nn);

      return ret; 
    }


  template<typename TT, typename ST>
    typename SEMBLE::PromoteEnsemVec<ST>::Type  LLSQMultiData<TT,ST>::get_row_ensem(const int row) const
    {
      SEMBLE::SembleVector<ST> foo;
      foo = this->get_row_semble(row); 
      typename SEMBLE::PromoteEnsemVec<ST>::Type ret; 
      ret.resize(foo.getB());
      ret.resizeObs(foo.getN());
      for(int nn = 0; nn < foo.getN(); ++nn)
        ENSEM::pokeObs(ret,foo.getEnsemElement(nn),nn);

      return ret; 

    }


  template<typename TT, typename ST>
    const TT&  LLSQMultiData<TT,ST>::get_tag(const int idx) const
    {
      return m_tags.at(idx); 
    }


  template<typename TT, typename ST>
    TT&  LLSQMultiData<TT,ST>::get_tag(const int idx)
    {
      return m_tags.at(idx);
    }


  template<typename TT, typename ST>
    void  LLSQMultiData<TT,ST>::set_tag(const int idx, const TT &tt)
    {
      if(idx > m_tags.size())
      {
        std::cerr << __func__ << ": error out of bounds" << std::endl;
        exit(1);
      }

      m_tags[idx] = tt; 
    } 

template<typename TT, typename ST>
void LLSQMultiData<TT,ST>::splash_tags(void) const
{
  typename std::vector<TT>::const_iterator tt; 
  for(tt = m_tags.begin(); tt != m_tags.end(); ++tt)
    std::cout << tt->file_id << std::endl;
}





} // namespace radmat




#endif /* LLSQ_MULTI_TIME_DATA_H */
