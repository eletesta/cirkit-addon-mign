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
 * @file xmg_io.hpp
 *
 * @brief Utilities for the flow map majority and almost majority 
 *
 * @author Eleonora Testa 
 * @since  2.3
 */

#ifndef MIGN_UTILS_MAJTT_HPP
#define MIGN_UTILS_MAJTT_HPP

#include <string>
#include <ostream>

#include <core/properties.hpp>
#include <classical/mign/mign.hpp>
#include <classical/utils/truth_table_utils.hpp>

namespace cirkit
{
tt calculate_tt_reminder_larger (tt func, tt tt_p); 
tt calculate_tt_reminder_smaller (tt func, tt tt_p); 
bool all_smaller(tt func, tt tt_p); 
bool all_larger(tt func, tt tt_p);
std::pair<unsigned,boost::dynamic_bitset<>> is_maj (tt func); 
std::tuple<unsigned,boost::dynamic_bitset<>,tt> is_almost_maj (tt func); 
}

#endif

// Local Variables:
// c-basic-offset: 2
// eval: (c-set-offset 'substatement-open 0)
// eval: (c-set-offset 'innamespace 0)
// End:
