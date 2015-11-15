#ifndef SPECTRUM_STATE_H
#define SPECTRUM_STATE_H 

#include "radmat/data_representation/data_representation.h"
#include "radmat/utils/mink_qsq.h"
#include <string>

namespace radmat
{

  // what lattice 
  enum LATTICE
  {
    L_3F_743
  };

  // states we currently care about 
  enum STATE
  {
    PI,
    PIp,
    PIpp,
    RHO,
    RHOp,
    RHOpp,
    RHOppp,
    A1
  };

  // translate enums to strings 
  struct enumStateHandler
  {
    static std::string handle(const int enumVal)
    {
      switch( enumVal )
      {
        case STATE::PI:
          return std::string("pi");
        case STATE::PIp:
          return std::string("pip");
        case STATE::PIpp:
          return std::string("pipp");
        case STATE::RHO:
          return std::string("rho");
        case STATE::RHOp:
          return std::string("rhop");
        case STATE::RHOpp:
          return std::string("rhopp");
        case STATE::RHOppp:
          return std::string("rhoppp");
        case STATE::A1:
          return std::string("a1");

        default : 
          return std::string("unknown"); 
      }
    }
  };

  // primitive 
  struct SpectrumState_p
  {
    virtual ~SpectrumState_p() {}
    virtual rHandle<LorentzRep> rep() const = 0; 
    virtual double rest_mass() const = 0; 
    virtual double boost_mass(const ADATXML::Array<int> &p) const = 0;
    virtual double mom_fac() const = 0; 
    virtual std::string name() const = 0; 
    virtual SpectrumState_p * clone() const = 0; 
  };

  // template instantiation of parameters  
  template< class CONT_REP ,
    int ENUM_STATE, 
    int ENUM_LATTICE,
    int MOM_FAC,           // times 100
    int REST_MASS,         // times 10000
    class ENUM_HANDLER  > 
      struct SpectrumState
      : public SpectrumState_p 
      {
        virtual ~SpectrumState() {}
        virtual rHandle<LorentzRep> rep() const 
        { 
          return rHandle<LorentzRep>( new CONT_REP() ); 
        }

        virtual double rest_mass() const 
        {
          return double(REST_MASS)/10000.; 
        }

        virtual double mom_fac() const 
        {
          return double(MOM_FAC) / 100.;
        }

        virtual double boost_mass(const ADATXML::Array<int> &p) const
        { 
          double rm = this->rest_mass(); 
          double mf = this->mom_fac(); 
          return sqrt(rm*rm + mf*mf*( p[0]*p[0] + p[1]*p[1] + p[2]*p[2])); 
        }

        virtual std::string name() const
        {
          return ENUM_HANDLER::handle(ENUM_STATE); 
        }

        virtual SpectrumState_p* clone() const 
        {
          return new SpectrumState(*this); 
        }

      };



}

#endif /* SPECTRUM_STATE_H */
