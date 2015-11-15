#ifndef DATA_REPRESENTATION_CUBIC_GROUPS_H
#define DATA_REPRESENTATION_CUBIC_GROUPS_H 

#include "data_representation_primitive_rep.h"
#include <string>
#include <vector>
#include "radmat/utils/handle.h"


namespace radmat
{

  ////////////////////////////////////////////////////////
  //
  //   SYMMETRY GROUPS  
  //

  template<typename T> 
  struct CubicRepEmbedType
  : public CubicRep_p
  {
    virtual ~CubicRepEmbedType() {}
    virtual std::string rep_id() const {return Stringify<T>();}
  };

  struct Oh; 
  struct D2; 
  struct D3; 
  struct D4; 

  REGISTER_STRINGIFY_TYPE(Oh);
  REGISTER_STRINGIFY_TYPE(D2);
  REGISTER_STRINGIFY_TYPE(D3);
  REGISTER_STRINGIFY_TYPE(D4);

  struct Oh : public CubicRepEmbedType<Oh> {};
  struct D2 : public CubicRepEmbedType<D2> {};
  struct D3 : public CubicRepEmbedType<D3> {};
  struct D4 : public CubicRepEmbedType<D4> {};


  
  ////////////////////////////////////////////////////////
  //
  //   REST IRREPS 
  //


  // all the irreps -- Oh
  struct A1;
  struct A2;
  struct T1;
  struct T2;
  struct E;

  REGISTER_STRINGIFY_TYPE(A1);
  REGISTER_STRINGIFY_TYPE(A2);
  REGISTER_STRINGIFY_TYPE(T1);
  REGISTER_STRINGIFY_TYPE(T2);
  REGISTER_STRINGIFY_TYPE(E); 

  struct A1 : public CubicRepEmbedType<A1> { enum {DIM = 1}; };
  struct A2 : public CubicRepEmbedType<A2> { enum {DIM = 1}; };
  struct T1 : public CubicRepEmbedType<T1> { enum {DIM = 3}; };
  struct T2 : public CubicRepEmbedType<T2> { enum {DIM = 3}; };
  struct E  : public CubicRepEmbedType<E>  { enum {DIM = 2}; };

  ////////////////////////////////////////////////////////
  //
  //   FLIGHT IRREPS 
  //

  // lg variants 
  struct E2; 
  struct B1; 
  struct B2; 
  
  REGISTER_STRINGIFY_TYPE(E2);
  REGISTER_STRINGIFY_TYPE(B1);
  REGISTER_STRINGIFY_TYPE(B2);

  struct E2 : public CubicRepEmbedType<E2> { enum {DIM = 2}; }; 
  struct B1 : public CubicRepEmbedType<B1> { enum {DIM = 1}; }; 
  struct B2 : public CubicRepEmbedType<B2> { enum {DIM = 1}; }; 

  // helicities for tokens
  struct H0  { enum {LAMBDA = 0}; virtual ~H0() {} };
  struct H1  { enum {LAMBDA = 1}; virtual ~H1() {} };
  struct H2  { enum {LAMBDA = 2}; virtual ~H2() {} };
  struct H3  { enum {LAMBDA = 3}; virtual ~H3() {} };
  struct H4  { enum {LAMBDA = 4}; virtual ~H4() {} };

  REGISTER_STRINGIFY_TYPE(H0);
  REGISTER_STRINGIFY_TYPE(H1);
  REGISTER_STRINGIFY_TYPE(H2);
  REGISTER_STRINGIFY_TYPE(H3);
  REGISTER_STRINGIFY_TYPE(H4);


  ////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////

  
  // the common part of rest and flight
  struct CubicRep
    : public CubicRep_p
  {
    virtual ~CubicRep() {}
    virtual std::string rep_group() const = 0; 
    virtual std::string rep_irrep() const = 0; 
    virtual int dim() const = 0; 
    virtual CubicRep * clone() const = 0; 
  };

  
  // rows are 1 based
  template<class G, class IRREP>
  struct CubicRepRest
    : public CubicRep
  {
    virtual ~CubicRepRest() {}; 
    virtual std::string rep_group() const { return Stringify<G>(); }
    virtual std::string rep_irrep() const { return Stringify<IRREP>(); }
    virtual std::string rep_id() const { return rep_irrep(); }
    virtual int dim() const { return IRREP::DIM; } 
    virtual CubicRep* clone() const { return new CubicRepRest( *this ); }
  };


  // type is CubicRep_t , group is Oh, irrep is w/e, id is irrep
  struct A1Rep_t : public CubicRepRest<Oh,A1> { virtual ~A1Rep_t() {} };
  struct A2Rep_t : public CubicRepRest<Oh,A2> { virtual ~A2Rep_t() {} };
  struct T1Rep_t : public CubicRepRest<Oh,T1> { virtual ~T1Rep_t() {} };
  struct T2Rep_t : public CubicRepRest<Oh,T2> { virtual ~T2Rep_t() {} };
  struct ERep_t  : public CubicRepRest<Oh,E > { virtual ~ERep_t()  {} };



  ////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////

  // add information about helicity rep
  struct CubicRepFlight_p
    : public CubicRep
  {
    virtual ~CubicRepFlight_p() {}; 
    virtual std::string rep_h() const = 0; 
    virtual int helicity() const = 0;  
  };


  // rows are 1 based
  template<class HRep, class G, class IRREP>
    struct CubicRepFlight
    : public CubicRepFlight_p
  {
    virtual ~CubicRepFlight() {}; 
    virtual std::string rep_group() const { return Stringify<G>(); }
    virtual std::string rep_irrep() const { return Stringify<IRREP>(); }
    virtual std::string rep_h() const { return Stringify<HRep>(); }
    virtual int helicity() const { return HRep::LAMBDA; }
    virtual std::string rep_id() const { return rep_h() + rep_group() + rep_irrep(); }
    virtual int dim() const { return IRREP::DIM; } 
    virtual CubicRep* clone() const { return new CubicRepFlight(*this); }
  };


  // D4 
  struct H0D4A1Rep_t : public CubicRepFlight<H0,D4,A1> { virtual ~H0D4A1Rep_t() {} };
  struct H0D4A2Rep_t : public CubicRepFlight<H0,D4,A2> { virtual ~H0D4A2Rep_t() {} };
  struct H1D4E2Rep_t : public CubicRepFlight<H1,D4,E2> { virtual ~H1D4E2Rep_t() {} };
  struct H2D4B1Rep_t : public CubicRepFlight<H2,D4,B1> { virtual ~H2D4B1Rep_t() {} };
  struct H2D4B2Rep_t : public CubicRepFlight<H2,D4,B2> { virtual ~H2D4B2Rep_t() {} };
  struct H3D4E2Rep_t : public CubicRepFlight<H3,D4,E2> { virtual ~H3D4E2Rep_t() {} };
  struct H4D4A1Rep_t : public CubicRepFlight<H4,D4,A1> { virtual ~H4D4A1Rep_t() {} };
  struct H4D4A2Rep_t : public CubicRepFlight<H4,D4,A2> { virtual ~H4D4A2Rep_t() {} };


  // D2 
  struct H0D2A1Rep_t : public CubicRepFlight<H0,D2,A1> { virtual ~H0D2A1Rep_t() {} };
  struct H0D2A2Rep_t : public CubicRepFlight<H0,D2,A2> { virtual ~H0D2A2Rep_t() {} };
  struct H1D2B1Rep_t : public CubicRepFlight<H1,D2,B1> { virtual ~H1D2B1Rep_t() {} };
  struct H1D2B2Rep_t : public CubicRepFlight<H1,D2,B2> { virtual ~H1D2B2Rep_t() {} };
  struct H2D2A1Rep_t : public CubicRepFlight<H2,D2,A1> { virtual ~H2D2A1Rep_t() {} };
  struct H2D2A2Rep_t : public CubicRepFlight<H2,D2,A2> { virtual ~H2D2A2Rep_t() {} };
  struct H3D2B1Rep_t : public CubicRepFlight<H3,D2,B1> { virtual ~H3D2B1Rep_t() {} };
  struct H3D2B2Rep_t : public CubicRepFlight<H3,D2,B2> { virtual ~H3D2B2Rep_t() {} };
  struct H4D2A1Rep_t : public CubicRepFlight<H4,D2,A1> { virtual ~H4D2A1Rep_t() {} };
  struct H4D2A2Rep_t : public CubicRepFlight<H4,D2,A2> { virtual ~H4D2A2Rep_t() {} };


  // D3 
  struct H0D3A1Rep_t : public CubicRepFlight<H0,D3,A1> { virtual ~H0D3A1Rep_t() {} };
  struct H0D3A2Rep_t : public CubicRepFlight<H0,D3,A2> { virtual ~H0D3A2Rep_t() {} };
  struct H1D3E2Rep_t : public CubicRepFlight<H1,D3,E2> { virtual ~H1D3E2Rep_t() {} };
  struct H2D3E2Rep_t : public CubicRepFlight<H2,D3,E2> { virtual ~H2D3E2Rep_t() {} };
  struct H3D3A1Rep_t : public CubicRepFlight<H3,D3,A1> { virtual ~H3D3A1Rep_t() {} };
  struct H3D3A2Rep_t : public CubicRepFlight<H3,D3,A2> { virtual ~H3D3A2Rep_t() {} };
  struct H4D3E2Rep_t : public CubicRepFlight<H4,D3,E2> { virtual ~H4D3E2Rep_t() {} };



  namespace CubicRepresentationFactoryEnv
  {
    bool registerAll();
    std::vector<std::string> cubic_keys(); 
    rHandle<CubicRep> callFactory(const std::string &id); 
  }

} 


#endif /* DATA_REPRESENTATION_CUBIC_GROUPS_H */
