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
 * @file mign_dfs.hpp
 *
 * @brief MIGN DFS functions
 *
 * @author Mathias Soeken
 * @since  2.3
 */

#ifndef MIGN_DFS_HPP
#define MIGN_DFS_HPP

#include <boost/graph/depth_first_search.hpp>

#include <string>
#include <ostream>

#include <classical/mign/mign.hpp>

namespace cirkit
{

class mign_dfs_visitor : public boost::default_dfs_visitor
{
public:
  mign_dfs_visitor( const mign_graph& mign );

  virtual void finish_vertex( const mign_node& node, const mign_graph::graph_t& g );

  virtual void finish_input( const mign_node& node, const mign_graph& mign ) = 0;
  virtual void finish_constant( const mign_node& node, const mign_graph& mign ) = 0;
 // virtual void finish_xor_node( const mign_node& node, const mign_function& a, const mign_function& b, const mign_graph& mign ) = 0;
  virtual void finish_maj_node( const mign_node& node, const std::vector<mign_function>& cciaoo, const mign_graph& mign ) = 0;

private:
  const mign_graph& mign;
};

}

#endif

// Local Variables:
// c-basic-offset: 2
// eval: (c-set-offset 'substatement-open 0)
// eval: (c-set-offset 'innamespace 0)
// End:
