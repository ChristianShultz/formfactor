#ifndef OBJ_EXPR_T_H
#define OBJ_EXPR_T_H

#include <list>
#include <ostream>


namespace radmat
{

  template<typename Coefficient, typename Object>
    struct ObjExpr_t
    {
      ObjExpr_t() {}
      ObjExpr_t(const Coefficient &c, const Object &o)
        : m_coeff(c) , m_obj(o)
      {}

      Coefficient m_coeff;     //! some coefficient 
      Object m_obj;            //! some object
    };

  template<typename Coefficient, typename Object>
    struct ListObjExpr_t
    {
      typedef typename std::list<ObjExpr_t<Coefficient,Object> > List_t; //! shorthand
      typedef typename List_t::const_iterator const_iterator;            //! iterator type
      typedef ObjExpr_t<Coefficient,Object> ListObj_t; 
      typedef Coefficient Coeff_t;
      typedef Object Obj_t;

      ListObjExpr_t(){}
      ListObjExpr_t(const int &) {} // for itpp
      ListObjExpr_t(const ListObjExpr_t<Coefficient,Object> &o) : m_expr(o.m_expr) {}
      ListObjExpr_t(const ObjExpr_t<Coefficient,Object> &expr) {m_expr.push_back(expr);}
      ListObjExpr_t<Coefficient,Object>& operator=(const ListObjExpr_t<Coefficient,Object> &o)
      {
        if(this != &o) {m_expr = o.m_expr;} return *this;
      }
      void push_back(const ObjExpr_t<Coefficient,Object> &expr) {m_expr.push_back(expr);}
      void clear(void) {m_expr.clear();}
      const_iterator begin(void) const {return m_expr.begin();}
      const_iterator end(void) const {return m_expr.end();}
      int size(void) const {return m_expr.size();}

      List_t m_expr;
    };


  // a cast like function
  template<typename Coefficient, typename Object>
    ListObjExpr_t<Coefficient,Object> 
      convert_to_list(const ObjExpr_t<Coefficient,Object> &o)
      {
        return ListObjExpr_t<Coefficient,Object>(o); 
      }

  //-----------------------------------------------------------------------

  //! stream an ObjExpr_t
  template<typename C, typename O>
    std::ostream& operator<<(std::ostream &os , const ObjExpr_t<C,O> &o)
    {
      os << o.m_coeff << " X (" << o.m_obj << ")";
      return os;
    }

  //! stream a ListObjExpr_t
  template<typename C, typename O>
    std::ostream& operator<<(std::ostream &os , const ListObjExpr_t<C,O> &o)
    {
      typename ListObjExpr_t<C,O>::const_iterator it;
      for(it = o.m_expr.begin(); it != o.m_expr.end(); ++it)
        os << *it << "  ";

      return os;
    }


  //-----------------------------------------------------------------------

  //! negate an ObjExpr_t
  template<typename C, typename O>
    ObjExpr_t<C,O> operator-(const ObjExpr_t<C,O> &o)
    {
      return ObjExpr_t<C,O>(-o.m_coeff, o.m_obj);
    }

  //! negate a ListObjExpr_t 
  template<typename C, typename O>
    ListObjExpr_t<C,O> operator-(const ListObjExpr_t<C,O> &o)
    {
      ListObjExpr_t<C,O> dest;
      typename ListObjExpr_t<C,O>::const_iterator it;
      for(it = o.m_expr.begin(); it != o.m_expr.end(); ++it)
        dest.push_back(-*it); 
      return dest;
    }


  //-----------------------------------------------------------------------

  //! add two ObjExpr_t 
  template<typename C, typename O>
    ListObjExpr_t<C,O> operator+(const ObjExpr_t<C,O> &o1, const ObjExpr_t<C,O> &o2)
    {
      ListObjExpr_t<C,O> dest(o1);
      dest.push_back(o2);
      return dest;
    }

  //! add two ListObjExpr_t
  template<typename C, typename O>
    ListObjExpr_t<C,O> operator+(const ListObjExpr_t<C,O> &o1, const ListObjExpr_t<C,O> &o2)
    {
      ListObjExpr_t<C,O> dest;
      typename ListObjExpr_t<C,O>::const_iterator it;
      dest.m_expr = o1.m_expr;

      for(it = o2.m_expr.begin(); it != o2.m_expr.end(); ++it)
        dest.push_back(*it);

      return dest;
    }

  //! add a Obj_Expr_t to a ListObjExpr_t
  template<typename C, typename O>
  ListObjExpr_t<C,O> operator+(const ListObjExpr_t<C,O> &o1, const ObjExpr_t<C,O> &o2)
  {
    ListObjExpr_t<C,O> tmp(o2);
    return o1 + tmp; 
  }

  //! add a ListObj_Expr_t to a ObjExpr_t
  template<typename C, typename O>
  ListObjExpr_t<C,O> operator+(const ObjExpr_t<C,O> &o1, const ListObjExpr_t<C,O> &o2)
  {
    ListObjExpr_t<C,O> tmp(o1);
    return tmp + o2; 
  }



  //-----------------------------------------------------------------------

  //! subtract two ObjExpr_t 
  template<typename C, typename O>
    ListObjExpr_t<C,O> operator-(const ObjExpr_t<C,O> &o1, const ObjExpr_t<C,O> &o2)
    {
      return (o1 + (-o2));  
    }

  //! subtract two ListObjExpr_t
  template<typename C, typename O>
    ListObjExpr_t<C,O> operator-(const ListObjExpr_t<C,O> &o1, const ListObjExpr_t<C,O> &o2)
    {
      return (o1 + (-o2));
    }


  //-----------------------------------------------------------------------

  //! multiply ObjExpr_t 
  template<typename C, typename O>
    ObjExpr_t<C,O> operator*(const C &c, const ObjExpr_t<C,O> &o2)
    {
      return ObjExpr_t<C,O>(c*o2.m_coeff,o2.m_obj);  
    }

  //! multiply a ListObjExpr_t
  template<typename C, typename O>
    ListObjExpr_t<C,O> operator*(const C &c, const ListObjExpr_t<C,O> &o2)
    {
      ListObjExpr_t<C,O> dest;
      typename ListObjExpr_t<C,O>::const_iterator it;

      for(it = o2.m_expr.begin(); it != o2.m_expr.end(); ++it)
        dest.push_back(c**it);

      return dest;
    }


  //-----------------------------------------------------------------------

  template<typename C, typename O>
    ObjExpr_t<C,O> operator/(const ObjExpr_t<C,O> &o1, const C &c)
    {
      return ObjExpr_t<C,O>(c*o1.m_coeff, o1.m_obj);
    } 

  template<typename C, typename O>
    ListObjExpr_t<C,O> operator/(const ListObjExpr_t<C,O> &o1, const C &c)
    {
      ListObjExpr_t<C,O> dest;
      typename ListObjExpr_t<C,O>::const_iterator it;

      for(it = o1.m_expr.begin(); it != o1.m_expr.end(); ++it)
        dest.push_back(*it/c);

      return dest;
    } 



} // namespace radmat




#endif /* OBJ_EXPR_T_H */
