#ifndef MIGN_REWRITE_HPP
#define MIGN_REWRITE_HPP

#include <functional>
#include <map>


#include<classical/mign/mign.hpp>
#include<core/properties.hpp>


namespace cirkit
{

//using mig_maj_rewrite_func_t = std::function<mig_function(mig_graph&, const mig_function&, const mig_function&, const mig_function&)>;

mign_graph mign_rewrite_top_down( const mign_graph& mign,
                                const properties::ptr& settings,
                                const properties::ptr& statistics ); 


}

#endif