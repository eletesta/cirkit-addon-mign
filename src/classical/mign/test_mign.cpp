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

#include "test_mign.hpp"

#include <fstream>
#include <iostream>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/topological_sort.hpp>
#include <boost/range/iterator_range.hpp>

#include <core/utils/graph_utils.hpp>

namespace cirkit
{

/******************************************************************************
 * Types                                                                      *
 ******************************************************************************/

struct my_vertex_info
{
  my_vertex_info() {}
  my_vertex_info( const std::string& name ) : name( name ) {}

  std::string name;
  bool        value = false;
};

std::ostream& operator<<( std::ostream& os, const my_vertex_info& v )
{
  os << v.name;
  return os;
}

struct my_edge_info
{
  bool complemented = false;
};

// using my_graph_t = boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS,
//                                          boost::property<boost::vertex_name_t, my_vertex_info>,
//                                          boost::property<boost::edge_name_t, my_edge_info>>;
// using my_vertex_t = boost::graph_traits<my_graph_t>::vertex_descriptor;
// using my_edge_t = boost::graph_traits<my_graph_t>::edge_descriptor;

using my_graph_t = digraph_t<boost::property<boost::vertex_name_t, my_vertex_info>, boost::property<boost::edge_name_t, my_edge_info>>;
using my_vertex_t = vertex_t<my_graph_t>;
using my_edge_t = edge_t<my_graph_t>;

struct my_edge_writer
{
  my_edge_writer( const my_graph_t& g ) : g( g ) {}

  void operator()( std::ostream& os, const my_edge_t& edge ) const
  {
    auto edge_info = boost::get( boost::edge_name, g );

    if ( edge_info[edge].complemented )
    {
      os << "[style=dashed,color=red]";
    }
  }

private:
  const my_graph_t& g;
};

/******************************************************************************
 * Private functions                                                          *
 ******************************************************************************/

/******************************************************************************
 * Public functions                                                           *
 ******************************************************************************/

void test_boost_graph_library()
{
  my_graph_t g;

  auto a = add_vertex( g );
  auto b = add_vertex( g );
  auto c = add_vertex( g );
  auto d = add_vertex( g );

  auto vertex_info = boost::get( boost::vertex_name, g );
  vertex_info[a] = my_vertex_info( "A" );
  vertex_info[b] = my_vertex_info( "B" );
  vertex_info[c] = my_vertex_info( "C" );
  vertex_info[d] = my_vertex_info( "D" );

  std::vector<double> weights( boost::num_vertices( g ) );
  weights[a] = 1.2;
  weights[b] = 3.2;
  weights[c] = 6.7;

  auto e1 = add_edge( c, a, g ).first;
  auto e2 = add_edge( c, b, g ).first;
  auto e3 = add_edge( b, d, g ).first;
  add_edge( a, d, g ).first;

  auto edge_info = boost::get( boost::edge_name, g );
  edge_info[e1].complemented = true;
  edge_info[e2].complemented = false;
  edge_info[e3].complemented = true;

  //std::cout << boost::num_vertices( g ) << std::endl;
  //std::cout << boost::num_edges( g ) << std::endl;

  std::ofstream os( "/tmp/test.dot", std::ostream::out );
  boost::write_graphviz( os, g, boost::make_label_writer( vertex_info ), my_edge_writer( g ) );
  os.close();

  /* precompute parents */
  auto inedges = precompute_ingoing_edges( g );
  auto indegrees = precompute_in_degrees( g );

  /* traditional way (before C++11) */
  boost::graph_traits<my_graph_t>::vertex_iterator it, itEnd;
  boost::tie( it, itEnd ) = boost::vertices( g );

  for ( ; it != itEnd; ++it )
  {
    std::cout << "vertex " << *it << " has name " << vertex_info[*it].name << std::endl;
  }

  for ( const auto& v : boost::make_iterator_range( boost::vertices( g ) ) )
  {
    std::cout << "vertex " << v << " has name " << vertex_info[v].name << " (C++11 way)" << std::endl;

    std::cout << "- out degree: " << boost::out_degree( v, g ) << std::endl;
    std::cout << "- in degree: " << indegrees[v] << std::endl;

    for ( const auto& w : boost::make_iterator_range( boost::adjacent_vertices( v, g ) ) )
    {
      std::cout << "- child is " << vertex_info[w].name << std::endl;
    }

    auto itIn = inedges.find( v );
    if ( itIn != inedges.end() )
    {
      for ( const auto& e : itIn->second )
      {
        std::cout << "- parent is " << vertex_info[boost::source( e, g )].name << std::endl;
      }
    }
  }

  for ( const auto& e : boost::make_iterator_range( boost::edges( g ) ) )
  {
    auto v = boost::source( e, g );
    auto w = boost::target( e, g );
    std::cout << vertex_info[v].name << " is connected to " << vertex_info[w].name;
    if ( edge_info[e].complemented )
    {
      std::cout << " and is complemented";
    }
    std::cout << std::endl;
  }

  std::vector<my_vertex_t> topsort( boost::num_vertices( g ) );
  boost::topological_sort( g, topsort.begin() );

  std::cout << "topological order:";
  for ( const auto& v : topsort )
  {
    std::cout << " " << vertex_info[v].name;
  }
  std::cout << std::endl;
}

}

// Local Variables:
// c-basic-offset: 2
// eval: (c-set-offset 'substatement-open 0)
// eval: (c-set-offset 'innamespace 0)
// End:
