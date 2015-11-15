/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name :

 * Purpose :

 * Creation Date : 12-10-2012

 * Last Modified : Thu 03 Apr 2014 05:05:57 PM EDT

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/




#include "redstar_invert_subduction.h"
#include "radmat/utils/pow2assert.h"
#include "radmat/utils/handle.h"
#include "radmat/utils/polarization_tensors.h"
#include "radmat/utils/printer.h"
#include "radmat/utils/tokenize.h"
#include "hadron/subduce_tables_factory.h"
#include "semble/semble_meta.h"

#include <sstream>
#include <vector>
#include <algorithm>
#include <list>
#include <complex>
#include <exception>


namespace  // a bunch of local stuff to make my life easier
{

  Hadron::SubduceTable* callSubduceFactory(const std::string &id)
  {
    return Hadron::TheSubduceTableFactory::Instance().createObject(id);
  }

  std::vector<std::string> getSubTableKeys(void)
  {
    return Hadron::TheSubduceTableFactory::Instance().keys();
  }

  bool match_string(const std::string &pattern, const std::string &text)
  {
    return std::search(text.begin(), text.end(), pattern.begin(), pattern.end()) != text.end();
  }

  std::vector<std::string> tokenize(const std::string &s,
      const std::string &delim, 
      const bool keep_empty_tokens=false)
  {
    std::vector<std::string> tokens;

    if(delim.empty())
      return std::vector<std::string>(1,s);

    std::string::const_iterator front,back;
    front = s.begin();
    const std::string::const_iterator end = s.end();

    while(true)
    {

      back = std::search(front,end,delim.begin(),delim.end());
      std::string token(front,back);
      if(keep_empty_tokens || !!!token.empty())
        tokens.push_back(token);

      if(back == end)
        break;

      front = back + delim.size();
    }

    return tokens; 
  }



  std::vector<std::string> getSubTableBosonKeys(void)
  {
    std::vector<std::string> bosons;
    std::vector<std::string> all_keys(getSubTableKeys());
    std::vector<std::string>::const_iterator it;
    std::string dont_match("o"); // all of the fermion labels feature something like J1o2 for spin half
    for(it = all_keys.begin(); it != all_keys.end(); ++it)
      if(!!!match_string(dont_match,*it))
        bosons.push_back(*it);
    return bosons;
  }

  std::vector<std::string> getBosonKeysPattern(const std::string &pattern)
  {
    std::vector<std::string> boson_keys = getSubTableBosonKeys();
    std::vector<std::string> oct_keys;
    std::vector<std::string>::const_iterator it;
    for(it = boson_keys.begin(); it != boson_keys.end(); ++it)
      if(match_string(pattern,*it))
        oct_keys.push_back(*it);
    return oct_keys;
  }

  std::vector<std::string> getBosonKeys(const std::string &group)
  {
    if(group == "Oh")
      return getBosonKeysPattern(std::string("J"));
    if(group == "D4")
      return getBosonKeysPattern(group);
    if(group == "D3")
      return getBosonKeysPattern(group);
    if(group == "D2") 
      return getBosonKeysPattern(group);
    if(group == "C4nm0") 
      return getBosonKeysPattern(group);
    if(group == "C4nnm")
      return getBosonKeysPattern(group);

    std::cerr << __func__ <<": ERROR: the group " << group 
      << " is not a supported/available symmetry group" << std::endl;
    exit(1);
  }

  bool inFlight(const std::string &group)
  {
    if(group == "Oh")
      return false;
    if(group == "D4")
      return true;
    if(group == "D3")
      return true;
    if(group == "D2") 
      return true;
    if(group == "C4nm0") 
      return true;
    if(group == "C4nnm")
      return true;

    std::cerr << __func__ << ": ERROR: the group " << group 
      << " is not a supported/available symmetry group" << std::endl;
    exit(1);
  }

  //! move from (-j..0..j) to (1..(2j + 1)) indexing using the convention that the most positive helicity is mapped to 1
  int remapHelicity_1based(const radmat::ContinuumBosonExprPrimitive &expr)
  {
    const int H = expr.H; 
    const int J = expr.J;

    if( abs(H) > abs(J) )
    {
      std::cerr << __func__ << ": error: out of bounds" << std::endl;
      exit(1);
    }
                
    return ::radmat::remapHelicity_1based(H,J); 
  }

  struct irrepKey
  {
    irrepKey(const std::string rep , const std::string sub)
      : irrep(rep) , subduce_key(sub)
    {}

    std::string irrep;
    std::string subduce_key;
  };

  // get all of the irreps that the cont operator went into dependent on whatever group it was in
  std::vector<irrepKey> getIrreps(const radmat::ContinuumBosonExprPrimitive &expr)
  {
    const std::string delimiter1("->");
    const std::string delimiter2(",1");
    const std::vector<std::string> group_keys(getBosonKeys(expr.group));
    std::stringstream parse_id_stream; 

    if(inFlight(expr.group))
    {
      parse_id_stream << "H" << abs(expr.H);

      if(abs(expr.H) == 0)   // N.B. we are using the eta_p thing. eta_p = P(-1)^J
      {
        int p = expr.parity ? 1 : -1; 
        int eta_p = (expr.J % 2 == 0) ? p : -p;
        (eta_p > 0) ? parse_id_stream << "+" : parse_id_stream << "-";   
      }
    }
    else
      parse_id_stream << "J" << expr.J;     

    std::string parse_id = parse_id_stream.str();

    std::vector<irrepKey> irreps;  
    std::vector<std::string>::const_iterator it;

    for(it = group_keys.begin(); it != group_keys.end(); ++it)
    {
      if(match_string(parse_id,*it))
      {
        //   std::cout << __func__ << ": matched " << *it << std::endl;
        std::vector<std::string> tokens = tokenize(*it,delimiter1);
        irreps.push_back( irrepKey(tokenize(tokens[1],delimiter2)[0], *it) );
      }
    }
    return irreps;
  }


  radmat::LatticeIrrepExpr_t makeLatticeIrrepExpr(const ENSEM::Complex &coefficient,
      const std::string &group, 
      const std::string &irrep,
      const int row)
  {
    return radmat::LatticeIrrepExpr_t(coefficient, radmat::LatticeExprPrimitive(group,irrep,row));
  }

  // the patterns for this one are easy and self evident
  radmat::ListLatticeIrrepExpr_t invertSubduction_rest(const radmat::ContinuumBosonExprPrimitive & expr)
  {
    std::vector<irrepKey> irreps = getIrreps(expr);
    std::vector<irrepKey>::const_iterator it;
    radmat::ListLatticeIrrepExpr_t lattice_expr;

    const int Hbase1 = remapHelicity_1based(expr);

    for(it = irreps.begin(); it != irreps.end(); ++it)
    {
      radmat::rHandle<Hadron::SubduceTable> subduce(callSubduceFactory(it->subduce_key)); 
      int rep_bound = subduce->dimt();

      for(int irrep_row = 1; irrep_row <= rep_bound; ++irrep_row)
      {

#if 0        
           std::cout <<"group: " << expr.group << "\n" 
           << "irrep: " << it->irrep << "\n"
           << "rep_bound: " << rep_bound << "\n"
           << "Hbase1: " << Hbase1 << "\n"
           << "irrep_row: " << irrep_row << "\n"
           << "subduce = " << (*subduce)(irrep_row,Hbase1) << "\n"
           << std::endl;
#endif
          

        if(ENSEM::toBool(ENSEM::localNorm2( (*subduce)(irrep_row,Hbase1))  > ENSEM::Real(0.0)))
          lattice_expr.push_back(
              makeLatticeIrrepExpr(
                (*subduce)(irrep_row,Hbase1),
                expr.group,
                it->irrep ,
                irrep_row));
      } 
    }

    return lattice_expr; 
  }

  // the patterns for this one suck 
  radmat::ListLatticeIrrepExpr_t invertSubduction_flight(const radmat::ContinuumBosonExprPrimitive & expr)
  {

    int parity = (expr.parity) ? 1 : -1;
    int eta_p = (expr.J % 2 == 0) ? parity : -parity;  

    std::vector<irrepKey> irreps = getIrreps(expr);
    std::vector<irrepKey>::const_iterator it;
    radmat::ListLatticeIrrepExpr_t lattice_expr;

#if 0
    std::cout << __PRETTY_FUNCTION__ << std::endl; 
    std::cout << radmat::toString( expr ) << std::endl;
    for(it = irreps.begin(); it != irreps.end(); ++it)
    {
      std::cout << it->irrep << std::endl;
    }
#endif

    for(it = irreps.begin(); it != irreps.end(); ++it)
    {
      radmat::rHandle<Hadron::SubduceTable>  subduce( callSubduceFactory(it->subduce_key) );
      int rep_bound = subduce->dimt();

      if(abs(expr.H) == 0)
      {
        // paranoia 
        if(rep_bound != 1)
        {
          std::cerr << __func__ << ": Error, it looks like a Helicity zero piece is going into something"
            << " other than an A irrep which is wrong at the point of writing this code" << std::endl;
          exit(1);
        }
        if(irreps.size() != 1)
        {
          std::cerr << __func__ << ": Error, it looks like the Helicity zero piece is getting split across "
            << "different irreps which doesn't make sense." << std::endl; 
        }

        lattice_expr.push_back(
            makeLatticeIrrepExpr(
              (*subduce)(1,1),
              expr.group,
              it->irrep ,
              1));
      }
      else if (expr.H > 0) // positive helicity
      {
        for(int irrep_row = 1; irrep_row <= rep_bound; ++irrep_row)
        {
          lattice_expr.push_back(
              makeLatticeIrrepExpr(
                (*subduce)(irrep_row,1),
                expr.group,
                it->irrep ,
                irrep_row));

        }
      }
      else if (expr.H < 0) // negative helicity
      {
        for(int irrep_row = 1; irrep_row <= rep_bound; ++irrep_row)
        {
          lattice_expr.push_back(
              makeLatticeIrrepExpr(
                ENSEM::Real(double(eta_p))*(*subduce)(irrep_row,2),
                expr.group,
                it->irrep ,
                irrep_row));

        }
      }
      else
      {
        std::cerr << __func__ << " Something went wrong, probably negative zero helicity" << std::endl;
        exit(1);
      }

    } // loop over irrep strings 

    return lattice_expr;
  }




} // namespace anonymous







namespace radmat
{

  namespace InvertSubductionEnv
  {
    namespace
    {
      bool local_reg = false; 
    }

    bool registerAll(void)
    {
      printer_function<debug_print_reg_all>(__PRETTY_FUNCTION__); 
      bool success = true; 

      if(!!! local_reg )
        success &= Hadron::SubduceTableEnv::registerAll();

      if( !!! success )
      {
        throw std::string("failed to reg Hadron::SubduceTableEnv"); 
      }

      return success; 
    }
  } 


  //! write this into a string for the factory
  std::string toString(const ContinuumBosonExprPrimitive &expr)
  {
    std::stringstream ss;

    ss << "_J_" << expr.J << "_H_" << expr.H << "_parity_";
    if(expr.parity)
      ss << "p";
    else
      ss << "m";

    ss << "_" << expr.group;

    return ss.str();
  }

  //! stream a ContinuumExprPrimitive
  std::ostream& operator<<(std::ostream &os, const ContinuumBosonExprPrimitive &expr)
  {
    os << toString(expr);
    return os;
  }

  //! write into a string
  std::string toString(const LatticeExprPrimitive &p)
  {
    std::stringstream ss;
    ss << "_" << p.group << "_" << p.irrep << "_row_" << p.row;
    return ss.str();
  }

  //! stream a LatticeExprPrimitive
  std::ostream& operator<<(std::ostream &os, const LatticeExprPrimitive &p)
  {
    os << toString(p);
    return os;
  }



  ListLatticeIrrepExpr_t invertSubduction(const ContinuumBosonExprPrimitive & expr)
  {
    if(inFlight(expr.group))
      return invertSubduction_flight(expr);
    else
      return invertSubduction_rest(expr);
  }


  LatticeIrrepExpr_t conj(const LatticeIrrepExpr_t &expr)
  {
    return LatticeIrrepExpr_t(ENSEM::conj(expr.m_coeff) , expr.m_obj);
  }


  ListLatticeIrrepExpr_t conj(const ListLatticeIrrepExpr_t &expr)
  {
    ListLatticeIrrepExpr_t dest; 
    ListLatticeIrrepExpr_t::const_iterator it; 
    for(it = expr.begin(); it != expr.end(); ++it)
      dest = dest + conj(*it); 

    return dest; 
  }

} // namespace radmat
