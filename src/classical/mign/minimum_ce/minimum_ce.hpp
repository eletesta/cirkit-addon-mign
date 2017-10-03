#ifndef MINIMUM_CE_HPP
#define MINIMUM_CE_HPP

#include <core/properties.hpp>
#include <classical/mign/mign.hpp>
#include <classical/utils/truth_table_utils.hpp>

namespace cirkit
{
    
mign_graph min_ce_with_sat (mign_graph mig, const tt& spec,
                              const properties::ptr& settings = properties::ptr(),
                              const properties::ptr& statistics = properties::ptr());
							  
}

#endif
