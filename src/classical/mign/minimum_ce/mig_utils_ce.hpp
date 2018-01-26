#ifndef MIG_UTILS_CE_HPP
#define MIG_UTILS_CE_HPP

#include <iostream>
#include <map>

#include "classical/mig/mig.hpp"

namespace cirkit
{


unsigned mig_print_stats_depth( const mig_graph& mig, std::ostream& os = std::cout );
unsigned depth_mig (const mig_graph mig);
std::vector<mig_function> get_children( const mig_graph& mig, const mig_node& node );

unsigned compute_level( const mig_graph& mig, const mig_node& node );

}

#endif