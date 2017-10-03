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
 * @file mign_flow_almost_maj.hpp
 *
 * @brief FlowMap algorithm for MIGns which considers MAJ-n and "Almost" MAJ-n
 *
 * @author Eleonora Testa 
 * @since  2.3
 */

#ifndef MIGN_FLOW_ALMOST_MAJ_HPP
#define MIGN_FLOW_ALMOST_MAJ_HPP

#include <core/properties.hpp>
#include <classical/mign/mign.hpp>

namespace cirkit
{

void mign_flow_almost_maj ( mign_graph& mign, const properties::ptr& settings = properties::ptr(), const properties::ptr& statistics = properties::ptr() );

}

#endif

// Local Variables:
// c-basic-offset: 2
// eval: (c-set-offset 'substatement-open 0)
// eval: (c-set-offset 'innamespace 0)
// End:
