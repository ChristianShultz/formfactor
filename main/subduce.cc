/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

* File Name : cubic_to_cart.cc

* Purpose :

* Creation Date : 29-04-2014

* Last Modified : Tue 29 Apr 2014 11:06:33 AM EDT

* Created By : shultz

_._._._._._._._._._._._._._._._._._._._._.*/

#include "radmat/register_all/register_all.h"
#include "radmat/data_representation/data_representation.h"
#include "radmat/redstar_interface/redstar_interface.h"
#include "radmat/utils/utils.h"
#include <string>
#include <sstream> 
#include <iostream>

using namespace radmat; 

const SubduceTableMap::irrep_sub_table*
  get_subduce_table(const rHandle<LorentzRep> &cont , const rHandle<CubicRep> &cub)
{
  std::string id = make_subduce_table_map_id(cont,cub); 
  std::cout << "subduce key " << id << std::endl;
  return TheSmarterSubduceTableMap::Instance().get_table(id); 
}

void print_subduce_table(int argc, char *argv[])
{

  if(argc != 4) 
  {
    std::cout << __func__ << " needs <cub> <cont> " << std::endl;
    exit(1); 
  }

  std::string contr,cubr; 
  { std::istringstream val(argv[2]); val >> cubr; }
  { std::istringstream val(argv[3]); val >> contr; }

  rHandle<LorentzRep> cont = LorentzRepresentationFactoryEnv::callFactory(contr); 
  rHandle<CubicRep> cub = CubicRepresentationFactoryEnv::callFactory(cubr); 

  const SubduceTableMap::irrep_sub_table *tab = get_subduce_table(cont,cub); 

  for(int i = 0; i < tab->size(); ++i)
  {
    int row = i +1; 
    std::cout << "row " << row << std::endl;
    SubduceTableMap::sub_list::const_iterator it; 
    for(it = tab->begin(row); it != tab->end(row); ++it)
     std::cout << "[" << it->first << "] x <" << it->second << ">   "; 
    std::cout << "\n\n";   
  }
}

// utilitiy 
typedef ListObjExpr_t<std::complex<double> , std::string> objs_t; 
itpp::Vec<objs_t> operator*(const itpp::Mat<std::complex<double> > &m, const itpp::Vec<objs_t> &v)
{
  itpp::Vec<objs_t> ret(m.rows()); 

  for(int row = 0; row < m.rows(); ++row)
    for(int col = 0; col < m.cols(); ++col)
    {
      if(std::norm(m(row,col)) > 0.0001)
      {
        // init 
        if(ret[row].m_expr.size() == 0)
          ret[row] =  m(row,col)*v(col);
        else
          ret[row] = ret[row] + m(row,col)*v(col);
      }
    }

  return ret; 
}


void J1m_subduce_components(int argc, char *argv[])
{

  if(argc != 6) 
  {
    std::cout << __func__ << " needs <cub> <px py pz> " << std::endl;
    exit(1); 
  }

  std::string cubr;
  ADATXML::Array<int> mom(3); 
  { std::istringstream val(argv[2]); val >> cubr; }
  { std::istringstream val(argv[3]); val >> mom[0]; }
  { std::istringstream val(argv[4]); val >> mom[1]; }
  { std::istringstream val(argv[5]); val >> mom[2]; }


  rHandle<LorentzRep> cont = LorentzRepresentationFactoryEnv::callFactory("J1m"); 
  rHandle<CubicRep> cub = CubicRepresentationFactoryEnv::callFactory(cubr); 

  const SubduceTableMap::irrep_sub_table *tab = get_subduce_table(cont,cub); 

  itpp::Vec<objs_t> cart(3),hel;
  cart[0] = objs_t::ListObj_t(std::complex<double>(1.,0.),"x");
  cart[1] = objs_t::ListObj_t(std::complex<double>(1.,0.),"y");
  cart[2] = objs_t::ListObj_t(std::complex<double>(1.,0.),"z");

  itpp::Mat<std::complex<double> > eps = eps3d(mom,true); 
  hel = eps * cart; 

  std::map<int,int> row_map; 
  row_map[1] = 0; 
  row_map[0] = 1; 
  row_map[-1] = 2; 

  for(int i = 0; i < tab->size(); ++i)
  {
    int row = i +1; 

    objs_t res; 
    

    std::cout << "row " << row << std::endl;
    SubduceTableMap::sub_list::const_iterator it; 
    for(it = tab->begin(row); it != tab->end(row); ++it)
      res = res + it->first * hel[row_map[it->second]]; 

    std::cout << res << "\n\n";   
  }

}


// typedef 
typedef void (*fptr)(int argc, char *argv[]);


// a map of operation names and function pointers
std::map<std::string , fptr> options; 

// init the map 
void init_options(void)
{
  options.insert(std::pair<std::string,fptr>("print_subduce_table",&print_subduce_table));
  options.insert(std::pair<std::string,fptr>("J1m_subduce_components",&J1m_subduce_components));
}

void handle_work(const std::string &op, int argc, char *argv[])
{

  init_options(); 

  if(options.find(op) == options.end())
  {
    std::cerr << " unrecognized op " << op 
      << " more intelligent choices are " << std::endl; 
    std::map<std::string , fptr>::const_iterator it; 
    for(it = options.begin(); it != options.end(); ++it)
      std::cerr << it->first << std::endl; 
    exit(1); 
  }

  fptr foo = options[op];

  foo(argc,argv); 
}


int main(int argc, char *argv[])
{
  AllFactoryEnv::registerAll(); 

  if(argc < 2)
  {
    std::cout << "error: <mode> stuff " << std::endl; 
    exit(1); 
  }

  std::string mode; 
  { std::istringstream val(argv[1]); val >> mode; }

  handle_work(mode,argc,argv); 

  return 0;
}

