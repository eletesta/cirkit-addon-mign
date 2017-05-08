/* CirKit: A circuit toolkit
 * Copyright (C) 2009-2015  University of Bremen
 * Copyright (C) 2015-2016  EPFL
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * @file threshold_synthesis.hpp
 *
 * @brief Threshold's synthesis algorithm by A. Neutzling
 *
 * @author Eleonora Testa 
 * @since  2.3
 */

#ifndef THRESHOLD_SYNTHESIS_HPP
#define THRESHOLD_SYNTHESIS_HPP

#include <core/properties.hpp>
#include <classical/utils/truth_table_utils.hpp>
#include <classical/mign/mign.hpp>

namespace cirkit
{

	
class resulting 
{
public:
		std::pair<unsigned, std::vector<unsigned>> t_and_w; 
		std::vector<bool> nega_un; 
};

mign_graph threshold_synthesis( const tt& func, 
                           const properties::ptr& settings = properties::ptr(),
                           const properties::ptr& statistics = properties::ptr() );

resulting compute_T_w ( tt& func); 

}

#endif

// Local Variables:
// c-basic-offset: 2
// eval: (c-set-offset 'substatement-open 0)
// eval: (c-set-offset 'innamespace 0)
// End:
