#ifndef MIGN_REWRITE_HPP
#define MIGN_REWRITE_HPP

#include <functional>
#include <map>


#include<classical/mign/mign.hpp>
#include<core/properties.hpp>


namespace cirkit
{

using mign_substitutes_map_t = boost::optional<std::unordered_map<mign_node, mign_function>>;

mign_graph mign_rewrite_top_down( const mign_graph& mign,
                                const properties::ptr& settings,
                                const properties::ptr& statistics ); 
								
mign_graph mign_rewrite_top_down_sub( const mign_graph& mign, const properties::ptr& settings,
								                                const properties::ptr& statistics );

}

#endif