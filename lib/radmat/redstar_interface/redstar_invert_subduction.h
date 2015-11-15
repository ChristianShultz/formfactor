#ifndef REDSTAR_INVERT_SUBDUCTION_H
#define REDSTAR_INVERT_SUBDUCTION_H


#include "radmat/utils/obj_expr_t.h"
#include "hadron/hadron_npart_npt_corr.h"
#include "adat/objfactory.h"
#include "adat/singleton.h"
#include "ensem/ensem.h"
#include <ostream>
#include <string>

namespace radmat
{

  //! a continuum expression for J^P , lambda and the lattice symmetry group we want to "undo"
  struct ContinuumBosonExprPrimitive
  {
    ContinuumBosonExprPrimitive(const int _J, const bool _parity, 
        const int _H, const std::string _group)
      : J(_J) , parity(_parity) , H(_H) , group(_group)
    {} 

    int J;           //! the total angular momentum -
    bool parity;        //! true is positive parity
    int H;           //! the helicity/z projection of spin at rest  (-J..0..J) 1 based is for fortran
    std::string group;  //! either "Oh", "D3" , "D2" , "D4" .. C when I get to it
  };

  //! write this into a string 
  std::string toString(const ContinuumBosonExprPrimitive &expr);

  //! stream a ContinuumExprPrimitive
  std::ostream& operator<<(std::ostream &os, const ContinuumBosonExprPrimitive &expr);


  //--------------------------------------------------------------------------------------

  //! a lattice expression for the irrep, row, and which symmetry group we are in
  struct LatticeExprPrimitive
  {
    LatticeExprPrimitive(const std::string _group , const std::string _irrep , const int _row) 
      : group(_group) , irrep(_irrep) , row(_row) 
    {}

    std::string group;  //! either "Oh" , "D3" , "D2" , "D4" .. 
    std::string irrep;  //! A1 B1 T1 etc 
    int row;            //! the row of the irrep 
  };

  //! write into a string
  std::string toString(const LatticeExprPrimitive &expr);

  //! stream a LatticeExprPrimitive
  std::ostream& operator<<(std::ostream &os, const LatticeExprPrimitive &p);


  //--------------------------------------------------------------------------------------

  //! shorthand
  typedef ObjExpr_t<ENSEM::Complex,LatticeExprPrimitive> LatticeIrrepExpr_t;

  //! shorthand
  typedef ListObjExpr_t<ENSEM::Complex,LatticeExprPrimitive> ListLatticeIrrepExpr_t; 

  //! determine the coefficients for inverting subduction
  ListLatticeIrrepExpr_t invertSubduction(const ContinuumBosonExprPrimitive &);

  //! conjugate an obj
  LatticeIrrepExpr_t conj(const LatticeIrrepExpr_t &); 

  //! conjugate the coefficients 
  ListLatticeIrrepExpr_t conj(const ListLatticeIrrepExpr_t &); 

  namespace InvertSubductionEnv
  {
    bool registerAll(void);
  }

} // namespace radmat




#endif /* REDSTAR_INVERT_SUBDUCTION_H */
