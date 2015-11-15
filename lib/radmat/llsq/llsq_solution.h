#ifndef LLSQ_SOLUTION_H
#define LLSQ_SOLUTION_H 

#include "semble/semble_meta.h"
#include "radmat/data_representation/data_representation.h"
#include "io/adat_xmlio.h"


namespace radmat
{

  template<typename T> 
    struct FormFacSolutions
    {
      FormFacSolutions() {}

      FormFacSolutions(const SEMBLE::SembleMatrix<T> &ff, 
          const std::vector<ThreePointDataTag> &i,
          const std::vector<std::string> &names)
        : FF_t(ff) , Ingredients(i) , Names(names)
      { }

      void append_ingredients(const ThreePointDataTag &t)
      {
        Ingredients.push_back(t); 
      }

      void append_name(const std::string &n)
      {
        Names.push_back(n); 
      }

      SEMBLE::SembleMatrix<T> FF_t; 
      std::vector<ThreePointDataTag> Ingredients; 
      std::vector<std::string> Names; 
    }; 

  template<typename T> 
    void 
    write(ADATIO::BinaryWriter &bin, const FormFacSolutions<T> &f)
    {
      int nr = f.FF_t.getN(); 

      write(bin,nr); 
      for(int ff = 0; ff < nr; ++ff)
      {
        typename SEMBLE::PromoteEnsemVec<T>::Type e; 
        SEMBLE::SembleVector<T> foo = f.FF_t.getRow(ff); 
        e.resize(foo.getB()); 
        e.resizeObs(foo.getN()); 

        for(int n = 0; n < foo.getN(); ++n)
          ENSEM::pokeObs(e,foo.getEnsemElement(n),n); 

        ENSEM::write(bin,e); 

        writeDesc(bin,f.Names[ff]); 

      }

      int nt = f.Ingredients.size(); 

      write(bin,nt); 
      for(int i =0; i < nt; ++i)
        write(bin,f.Ingredients[i]); 
    } 

  template<typename T> 
    void 
    read(ADATIO::BinaryReader &bin, FormFacSolutions<T> &out)
    {
      FormFacSolutions<T> f; 
      int nr;
      read(bin,nr); 

      for(int ff = 0; ff < nr; ++ff)
      {
        typename SEMBLE::PromoteEnsemVec<T>::Type e; 
        SEMBLE::SembleVector<T> foo; 
        ENSEM::read(bin,e); 

        std::string n; 
        readDesc(bin,n); 

        f.append_name(n); 

        foo = e; 
        if( ff == 0 ) // first append 
          f.FF_t.reDim(foo.getB(),0,foo.getN()); 
        f.FF_t.append_row(foo); 
      }

      int nt; 
      read(bin,nt); 

      f.Ingredients.resize(nt);
      for(int i =0; i < nt; ++i)
        read(bin,f.Ingredients[i]); 

      out = f; 
    } 

}



#endif /* LLSQ_SOLUTION_H */
