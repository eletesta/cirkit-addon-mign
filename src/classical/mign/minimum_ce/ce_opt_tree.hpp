#ifndef CE_OPT_TREE_HPP
#define CE_OPT_TREE_HPP

#include <functional>
#include <map>

#include "core/properties.hpp"
#include "classical/mign/mign.hpp"
#include "classical/mign/mign_utils.hpp"


namespace cirkit
{
    
std::map<mign_node, std::vector<std::pair<mign_function,unsigned>>> init_visited_table_ce( const mign_graph& mig, mign_graph& mig_new );

std::vector<std::pair<mign_function,unsigned>> comple (std::vector<std::pair<mign_function,unsigned>> ciao, bool complemented); 
    
std::vector<std::pair<mign_function,unsigned>> change_ce (mign_graph& mig_new, std::vector<std::pair<mign_function,unsigned>> a, std::vector<std::pair<mign_function,unsigned>> b, std::vector<std::pair<mign_function,unsigned>> c, bool outcompl, mign_graph mig, mign_node last);
    
}

#endif

// Local Variables:
// c-basic-offset: 2
// eval: (c-set-offset 'substatement-open 0)
// eval: (c-set-offset 'innamespace 0)
// End: