/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : build_correlators_bad_data_repository.cc

 * Purpose :

 * Creation Date : 01-11-2013

 * Last Modified : Wed 30 Apr 2014 01:26:10 PM EDT

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/


#include "radmat/database/database.h"
#include "construct_correlators_bad_data_repository.h"
#include "semble/semble_file_management.h"
#include "hadron/ensem_filenames.h"
#include "io/adat_xmlio.h"
#include <omp.h>
#include <vector>
#include <map>

namespace radmat
{
  void 
    BuildCorrsLocalBadDataRepo_t::insert(const int tid, 
        const std::vector<RadmatExtendedKeyHadronNPartIrrep_t> &v)
    {
#pragma omp critical
      {
        if(norms.find(tid) != norms.end())
          norms.find(tid)->second.insert(norms.find(tid)->second.begin(),v.begin(),v.end());
        else
          norms.insert(
              std::map<int,std::vector<RadmatExtendedKeyHadronNPartIrrep_t> >::value_type(tid,v));
      }
    }

  void BuildCorrsLocalBadDataRepo_t::insert(const int tid,
      const std::vector<Hadron::KeyHadronNPartNPtCorr_t> &v)
  {
#pragma omp critical
    {
      if(norms.find(tid) != norms.end())
        tpc.find(tid)->second.insert(tpc.find(tid)->second.begin(),v.begin(),v.end());
      else
        tpc.insert(
            std::map<int,std::vector<Hadron::KeyHadronNPartNPtCorr_t> >::value_type(tid,v));
    }        
  }


  std::vector<Hadron::KeyHadronNPartNPtCorr_t> 
    BuildCorrsLocalBadDataRepo_t::merge_bad_data_npt(void) const
  {
    typedef std::map<std::string,Hadron::KeyHadronNPartNPtCorr_t>::value_type value_type; 
    std::map<std::string,Hadron::KeyHadronNPartNPtCorr_t> uniq; 
    std::map<std::string,Hadron::KeyHadronNPartNPtCorr_t>::const_iterator it; 
    std::map<int,std::vector<Hadron::KeyHadronNPartNPtCorr_t > >::const_iterator it1; 
    std::vector<Hadron::KeyHadronNPartNPtCorr_t>::const_iterator it2; 
    std::vector<Hadron::KeyHadronNPartNPtCorr_t> ret; 

    for(it1 = tpc.begin(); it1 != tpc.end(); ++it1)
      for(it2 = it1->second.begin(); it2 != it1->second.end(); ++it2)
        uniq.insert(value_type(Hadron::ensemFileName(*it2),*it2));

    for(it = uniq.begin(); it != uniq.end(); ++it)
      ret.push_back(it->second); 

    return ret; 
  }

  std::vector<RadmatExtendedKeyHadronNPartIrrep_t> 
    BuildCorrsLocalBadDataRepo_t::merge_bad_data_norm(void) const
    {
      typedef std::map<std::string,RadmatExtendedKeyHadronNPartIrrep_t>::value_type value_type; 
      std::map<std::string,RadmatExtendedKeyHadronNPartIrrep_t> uniq; 
      std::map<std::string,RadmatExtendedKeyHadronNPartIrrep_t>::const_iterator it; 
      std::map<int,std::vector<RadmatExtendedKeyHadronNPartIrrep_t > >::const_iterator it1; 
      std::vector<RadmatExtendedKeyHadronNPartIrrep_t>::const_iterator it2; 
      std::vector<RadmatExtendedKeyHadronNPartIrrep_t> ret; 

      for(it1 = norms.begin(); it1 != norms.end(); ++it1)
        for(it2 = it1->second.begin(); it2 != it1->second.end(); ++it2)
          uniq.insert(value_type(fileName(*it2),*it2));

      for(it = uniq.begin(); it != uniq.end(); ++it)
        ret.push_back(it->second); 

      return ret; 

    }


  void BuildCorrsLocalBadDataRepo_t::dump_baddies(void) const
  {

    std::vector<Hadron::KeyHadronNPartNPtCorr_t> bcv = merge_bad_data_npt();
    std::vector<RadmatExtendedKeyHadronNPartIrrep_t> bnv =  merge_bad_data_norm();

    if(!!!bcv.empty())
    {
      std::cout << __func__ << ": printing missing three-point xml " << std::endl;
      ADATXML::XMLBufferWriter corrs;
      ADATXML::Array<Hadron::KeyHadronNPartNPtCorr_t> bc;

      bc.resize(bcv.size()); 
      for(unsigned int i = 0; i < bcv.size(); ++i)
        bc[i] = bcv[i];

      write(corrs,"NPointList",bc);

      std::string pth = SEMBLE::SEMBLEIO::getPath();
      SEMBLE::SEMBLEIO::makeDirectoryPath(pth + std::string("missing")); 

      std::ofstream out("missing/npt.list.xml");
      corrs.print(out);
      out.close();

      std::vector<Hadron::KeyHadronNPartNPtCorr_t>::const_iterator it; 

      out.open("missing/npt.ensemFileNames.list"); 
      for(it = bcv.begin(); it != bcv.end(); ++it)
        out << Hadron::ensemFileName(*it) << "\n";
      out.close(); 

    }

    if(!!!bnv.empty())
    {
      std::cout << __func__ << ": printing missing state info xml " << std::endl;
      ADATXML::XMLBufferWriter norms;
      ADATXML::Array<RadmatExtendedKeyHadronNPartIrrep_t> bn;

      std::vector<RadmatExtendedKeyHadronNPartIrrep_t>::const_iterator it;


      bn.resize(bnv.size());
      for(unsigned int i = 0; i < bnv.size(); ++i)
        bn[i] = bnv[i];

      write(norms,"BadNorms",bn);

      std::string pth = SEMBLE::SEMBLEIO::getPath();
      SEMBLE::SEMBLEIO::makeDirectoryPath(pth + std::string("missing")); 

      std::ofstream out("missing/normalizations.list.xml");
      norms.print(out);
      out.close(); 
    }

  }



} // radmat


