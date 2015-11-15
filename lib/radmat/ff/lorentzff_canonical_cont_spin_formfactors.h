#ifndef LORENTZFF_CANONICAL_CONT_SPIN_CANONICAL_FORMFACTORS_H
#define LORENTZFF_CANONICAL_CONT_SPIN_CANONICAL_FORMFACTORS_H 

#include "lorentzff_canonical_PiPi.h"
#include "lorentzff_canonical_PiPiStar.h"
#include "lorentzff_canonical_PiRho.h"
#include "lorentzff_canonical_RhoPi.h"
#include "lorentzff_canonical_RhoRho.h"
#include "lorentzff_canonical_one.h"
#include "radmat/utils/stringify.h"


namespace radmat
{
  // 0 x 0 
  ///////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////

  // spin holders 
  // <0p | j | 0p >
  struct J0pJ0p_diag; 
  REGISTER_STRINGIFY_TYPE( J0pJ0p_diag ); 

  struct J0pJ0p_diag : public PiPi<0,0>
  {
    virtual ~J0pJ0p_diag() {};
    virtual LorentzFFAbsBase_t* clone() const { return new J0pJ0p_diag() ; }
    virtual std::string reg_id() const { return Stringify< J0pJ0p_diag >(); }
  };

  // spin holders 
  // <0m | j | 0m>
  struct J0mJ0m_diag; 
  REGISTER_STRINGIFY_TYPE( J0mJ0m_diag ); 

  struct J0mJ0m_diag : public PiPi<0,0>
  {
    virtual ~J0mJ0m_diag() {};
    virtual LorentzFFAbsBase_t* clone() const { return new J0mJ0m_diag() ; }
    virtual std::string reg_id() const { return Stringify< J0mJ0m_diag >(); }
  };

  // spin holders 
  // <0p | j | 0p >
  struct J0pJ0p_tran; 
  REGISTER_STRINGIFY_TYPE( J0pJ0p_tran ); 

  struct J0pJ0p_tran : public PiPiStar<0,0>
  {
    virtual ~J0pJ0p_tran() {};
    virtual LorentzFFAbsBase_t* clone() const { return new J0pJ0p_tran() ; }
    virtual std::string reg_id() const { return Stringify< J0pJ0p_tran >(); }
  };

  // spin holders 
  // <0m | j | 0m>
  struct J0mJ0m_tran; 
  REGISTER_STRINGIFY_TYPE( J0mJ0m_tran ); 

  struct J0mJ0m_tran : public PiPiStar<0,0>
  {
    virtual ~J0mJ0m_tran() {};
    virtual LorentzFFAbsBase_t* clone() const { return new J0mJ0m_tran() ; }
    virtual std::string reg_id() const { return Stringify< J0mJ0m_tran >(); }
  };


  // 1 x 0 
  ///////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////

  // spin holders 
  // <1m | j | 0m >
  struct J1mJ0m_tran; 
  REGISTER_STRINGIFY_TYPE( J1mJ0m_tran ); 

  struct J1mJ0m_tran : public RhoPi<1,0>
  {
    virtual ~J1mJ0m_tran() {};
    virtual LorentzFFAbsBase_t* clone() const { return new J1mJ0m_tran() ; }
    virtual std::string reg_id() const { return Stringify< J1mJ0m_tran >(); }
  };

  // spin holders 
  // <1p | j | 0p >
  struct J1pJ0p_tran; 
  REGISTER_STRINGIFY_TYPE( J1pJ0p_tran ); 

  struct J1pJ0p_tran : public RhoPi<1,0>
  {
    virtual ~J1pJ0p_tran() {};
    virtual LorentzFFAbsBase_t* clone() const { return new J1pJ0p_tran() ; }
    virtual std::string reg_id() const { return Stringify< J1pJ0p_tran >(); }
  };

  // 0 x 1 
  ///////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////

  // spin holders 
  // <0m | j | 1m >
  struct J0mJ1m_tran; 
  REGISTER_STRINGIFY_TYPE( J0mJ1m_tran ); 

  struct J0mJ1m_tran : public PiRho<0,1>
  {
    virtual ~J0mJ1m_tran() {};
    virtual LorentzFFAbsBase_t* clone() const { return new J0mJ1m_tran() ; }
    virtual std::string reg_id() const { return Stringify< J0mJ1m_tran >(); }
  };


  // spin holders 
  // <0p | j | 1p >
  struct J0pJ1p_tran; 
  REGISTER_STRINGIFY_TYPE( J0pJ1p_tran ); 

  struct J0pJ1p_tran : public PiRho<0,1>
  {
    virtual ~J0pJ1p_tran() {};
    virtual LorentzFFAbsBase_t* clone() const { return new J0pJ1p_tran() ; }
    virtual std::string reg_id() const { return Stringify< J0pJ1p_tran >(); }
  };

  // 1 x 1 
  ///////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////


#if 1

  // use G1,G2,G3 decomp

  // spin holders 
  // <1p | j | 1p >
  struct J1pJ1p_diag; 
  REGISTER_STRINGIFY_TYPE( J1pJ1p_diag ); 

  struct J1pJ1p_diag : public RhoRho<1,1>
  {
    virtual ~J1pJ1p_diag() {};
    virtual LorentzFFAbsBase_t* clone() const { return new J1pJ1p_diag() ; }
    virtual std::string reg_id() const { return Stringify< J1pJ1p_diag >(); }
  };

  // spin holders 
  // <1m | j | 1m>
  struct J1mJ1m_diag; 
  REGISTER_STRINGIFY_TYPE( J1mJ1m_diag ); 

  struct J1mJ1m_diag : public RhoRho<1,1>
  {
    virtual ~J1mJ1m_diag() {};
    virtual LorentzFFAbsBase_t* clone() const { return new J1mJ1m_diag() ; }
    virtual std::string reg_id() const { return Stringify< J1mJ1m_diag >(); }
  };

#else

  // use multipole decomp

  // spin holders 
  // <1p | j | 1p >
  struct J1pJ1p_diag; 
  REGISTER_STRINGIFY_TYPE( J1pJ1p_diag ); 

  struct J1pJ1p_diag : public RhoRhoMultipole<1,1>
  {
    virtual ~J1pJ1p_diag() {};
    virtual LorentzFFAbsBase_t* clone() const { return new J1pJ1p_diag() ; }
    virtual std::string reg_id() const { return Stringify< J1pJ1p_diag >(); }
  };

  // spin holders 
  // <1m | j | 1m>
  struct J1mJ1m_diag; 
  REGISTER_STRINGIFY_TYPE( J1mJ1m_diag ); 

  struct J1mJ1m_diag : public RhoRhoMultipole<1,1>
  {
    virtual ~J1mJ1m_diag() {};
    virtual LorentzFFAbsBase_t* clone() const { return new J1mJ1m_diag() ; }
    virtual std::string reg_id() const { return Stringify< J1mJ1m_diag >(); }
  };

#endif 


  //  test classes 
  ////////////////////////////
  struct J1mJ1m_test;
  REGISTER_STRINGIFY_TYPE(J1mJ1m_test); 
  struct J1mJ1m_test : public ONE<1,1>
  {
    virtual ~J1mJ1m_test() {}; 
    virtual LorentzFFAbsBase_t* clone() const {return new J1mJ1m_test();}
    virtual std::string reg_id() const { return Stringify< J1mJ1m_test >();}
  };



} // radmat 

#endif /* LORENTZFF_CANONICAL_CONT_SPIN_CANONICAL_FORMFACTORS_H */
