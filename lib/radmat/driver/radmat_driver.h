#ifndef RADMAT_DRIVER_H
#define RADMAT_DRIVER_H


#include <string>
#include "radmat/utils/handle.h"
#include "radmat/construct_data/construct_correlators.h"
#include "radmat/llsq/llsq_multi_data.h"
#include "radmat_driver_props.h"
#include "radmat_single_q2_driver.h"

namespace radmat
{

  struct RadmatDriver
  {
    RadmatDriver(void) {}

    //! do a three point analysis and extract form factors    
    void run_program(const std::string &inifile);

    //! split work b/c qdp isn't playing nice 
    void xml_handler(const std::string &inifile, const std::string &mode); 

    //! just run up to the point of building xml 
    void build_xml(const std::string &inifile); 

    //! split up xml on p^2
    void build_xml_split_p2(const std::string &inifile);  

    //! only build canonical momenta at the insertion
    void build_xml_canon_p2(const std::string &inifile);  

    //! split up xml on p^2 then split into N sets 
    //    -- N is set in the .cc file 
    void build_xml_split_p2_N(const std::string &inifile);

    //! split up xml on octant 
    void build_xml_split(const std::string &inifile); 

    //! two point xml hack
    void build_xml_twopoint(const std::string &inifile);

    //    //! just figure out what disconnected graphs we have to nuke for redstar
    //    void nuke_graph(const std::string &inifile ,
    //                    const std::string &graph_db,
    //                    const std::string &nuke_xml_out);   
    //
    //    //! stick in stubs for gen_prop generation
    //    void build_stub_xml(const std::string &inifile); 

    private:

    bool read_ini, built_correlators, init_llsq, solved_llsq, fit_formfacs;
    bool chisq_analysis; 

    void init_false(void); 
    void check_exit_ini(void) const {check_exit(read_ini,__func__);}
    void check_exit_corrs(void) const {check_exit(built_correlators,__func__);}
    void check_exit_init_llsq(void) const {check_exit(init_llsq,__func__);}
    void check_exit_llsq(void) const {check_exit(solved_llsq,__func__);}
    void check_exit_fit(void) const {check_exit(fit_formfacs,__func__);}
    void check_exit_chisq(void) const {check_exit(chisq_analysis,__func__);}
    void check_exit(const bool &, const char *) const;

    // for running the full analysis 
    bool read_xmlini(const std::string &ini); 
    bool build_correlators(void);
    bool solve_llsq(void);
    bool fit_ffs(void); 
    bool do_chisq_analysis(void);  
    bool make_FF_of_Q2_plots(void);
    bool print_Q2_list(void);  


    RDriverProps_t m_ini;  
    ConstructCorrelators m_correlators; 
    std::vector<bool> good_qs;
    std::vector<rHandle<LLSQLatticeMultiData> > multi_lattice_data; 
    std::vector<RadmatSingleQ2Driver> linear_systems_of_Q2; 
  };


}


#endif /* RADMAT_DRIVER_H */
