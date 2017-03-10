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

#include "mign_dfs.hpp"

#include <boost/range/algorithm.hpp>

namespace cirkit
{

/******************************************************************************
 * Types                                                                      *
 ******************************************************************************/

/******************************************************************************
 * Private functions                                                          *
 ******************************************************************************/

/******************************************************************************
 * Public functions                                                           *
 ******************************************************************************/

mign_dfs_visitor::mign_dfs_visitor( const mign_graph& mign )
  : mign( mign )
{
}

void mign_dfs_visitor::finish_vertex( const mign_node& node, const mign_graph::graph_t& g )
{
  boost::default_dfs_visitor::finish_vertex( node, g );
  if ( mign.is_input( node ) )
  {
    if ( node == 0u )
    {
      finish_constant( node, mign );
    }
    else
    {
      finish_input( node, mign );
    }
  }
  else //if ( mign.is_maj( node ) )
  {
    const auto cciaoo = mign.children( node );
    finish_maj_node( node, cciaoo, mign );
  }
}

}

// Local Variables:
// c-basic-offset: 2
// eval: (c-set-offset 'substatement-open 0)
// eval: (c-set-offset 'innamespace 0)
// End:
