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
 * @brief MIG-n graph
 *
 * @author Eleonora Testa
 * @since  2.3
 */

#ifndef MIG_TO_MIGN_HPP
#define MIG_TO_MIGN_HPP


#include <classical/mign/mign.hpp>
#include <classical/mig/mig.hpp>


namespace cirkit 
{

mign_graph mig_to_mign_not_strash (const mig_graph& mig); 
mign_graph mig_to_mign (const mig_graph& mig); 
mig_graph  mign_to_mig (const mign_graph& mign); 

}

#endif

// Local Variables:
// c-basic-offset: 2
// eval: (c-set-offset 'substatement-open 0)
// eval: (c-set-offset 'innamespace 0)
// End:
