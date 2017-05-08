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
 * @file mign.hpp
 *
 * @brief Create the bitcount with nmig. poly depth
 *
 * @author Eleonora Testa
 * @since  2.3
 */

#ifndef BITCOUNT_HPP
#define BITCOUNT_HPP


#include <classical/mign/mign.hpp>
#include <boost/optional/optional.hpp>

namespace cirkit 
{

mign_function find_mign(const std::unordered_map<std::string, mign_function>& name_to_function, const std::string& name );
mign_graph bitcount (unsigned Nin); 
mign_graph exact_count (unsigned Nin, unsigned count); 
mign_graph threshold (unsigned Nin, unsigned threshold, unsigned polarity,std::vector<unsigned> weigths); 
mign_graph add2_luca (unsigned input);
mign_graph addm (unsigned m, unsigned n); 
mign_graph mult2 (unsigned int Nin); 
}

#endif

// Local Variables:
// c-basic-offset: 2
// eval: (c-set-offset 'substatement-open 0)
// eval: (c-set-offset 'innamespace 0)
// End:
