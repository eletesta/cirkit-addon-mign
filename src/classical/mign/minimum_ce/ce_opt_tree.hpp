#ifndef CE_OPT_TREE_HPP
#define CE_OPT_TREE_HPP

#include <functional>
#include <map>

#include "core/properties.hpp"
#include "classical/mig/mig.hpp"
#include "classical/mig/mig_utils.hpp"
#include "classical/mign/minimum_ce/mig_strash.hpp"



namespace cirkit
{

//mig_function rewrite_default_maj( mig_graph& mig_new, const mig_function& a, const mig_function& b, const mig_function& c );

mig_graph ce_opt_tree( const mig_graph mig,const properties::ptr& settings = properties::ptr(),
                                const properties::ptr& statistics = properties::ptr() );
    
std::map<mig_node, std::vector<std::pair<mig_function,unsigned>>> init_visited_table_ce( const mig_graph& mig, mig_graph& mig_new );

std::vector<std::pair<mig_function,unsigned>> comple (std::vector<std::pair<mig_function,unsigned>> ciao, bool complemented); 
    
std::vector<std::pair<mig_function,unsigned>>  change_ce (mig_graph& mig_new, std::vector<std::pair<mig_function,unsigned>> a, std::vector<std::pair<mig_function,unsigned>> b, std::vector<std::pair<mig_function,unsigned>> c, bool outcompl, mig_graph mig, mig_node last);
    
bool e_output (mig_graph old_mig, mig_node node); 

}

#endif

// Local Variables:
// c-basic-offset: 2
// eval: (c-set-offset 'substatement-open 0)
// eval: (c-set-offset 'innamespace 0)
// End: