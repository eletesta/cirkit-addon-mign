/* CirKit: A circuit toolkit
 * Copyright (C) 2009-2015  University of Bremen
 * Copyright (C) 2015-2017  EPFL
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

#include "mign_mffc.hpp"

#include <map>

#include <boost/graph/depth_first_search.hpp>
#include <boost/range/algorithm.hpp>

#include <classical/mign/mign_bitmarks.hpp>
#include <classical/mign/mign_dfs.hpp>
#include <classical/mign/mign_utils.hpp>

namespace cirkit
{

/******************************************************************************
 * Types                                                                      *
 ******************************************************************************/
struct mign_size_visitor : public mign_dfs_visitor
{
public:
  mign_size_visitor( const mign_graph& mign, const std::vector<mign_node>& support, unsigned& size )
    : mign_dfs_visitor( mign ),
      support( support ),
      size( size )
  {
  }

  void finish_input( const mign_node& node, const mign_graph& mign )
  {
  }

  void finish_constant( const mign_node& node, const mign_graph& mign )
  {
  }

  void finish_maj_node( const mign_node& node, const std::vector<mign_function>& a, const mign_graph& mign )
  {
    if ( boost::find( support, node ) == support.end() )
    {
      ++size;
    }
  }

private:
  const std::vector<mign_node>& support;
  unsigned& size;
};

struct mign_cone_visitor : public mign_dfs_visitor
{
public:
  mign_cone_visitor( const mign_graph& mign, const std::vector<mign_node>& support, std::vector<mign_node>& cone )
    : mign_dfs_visitor( mign ),
      support( support ),
      cone( cone )
  {
  }

  void finish_input( const mign_node& node, const mign_graph& mign )
  {
  }

  void finish_constant( const mign_node& node, const mign_graph& mign )
  {
  }

  void finish_maj_node( const mign_node& node, const std::vector<mign_function>& a, const mign_graph& mign )
  {
    if ( boost::find( support, node ) == support.end() )
    {
      cone.push_back( node );
    }
  }

private:
  const std::vector<mign_node>& support;
  std::vector<mign_node>& cone;
};

/******************************************************************************
 * Private functions                                                          *
 ******************************************************************************/

/* inspired by: https://bitbucket.org/alanmi/abc/src/9f0e7e81524337aeebec196723c23915c7354982/src/aig/gia/giaMffc.c?at=default&fileviewer=file-view-default */

unsigned mign_mffc_node_deref( mign_graph& mign, mign_node n )
{
  auto counter = 0u;

  if ( mign.is_input( n ) )
  {
    return 0;
  }

  for ( auto child : mign.children( n ) )
  {
    assert( mign.get_ref( child.node ) > 0u );
    if ( mign.dec_ref( child.node ) == 0u )
    {
      counter += mign_mffc_node_deref( mign, child.node );
    }
  }

  return counter + 1u;
}

unsigned mign_mffc_node_ref( mign_graph& mign, mign_node n )
{
  auto counter = 0u;

  if ( mign.is_input( n ) )
  {
    return 0;
  }

  for ( auto child : mign.children( n ) )
  {
    if ( mign.inc_ref( child.node ) == 0u )
    {
      counter += mign_mffc_node_ref( mign, child.node );
    }
  }

  return counter + 1u;
}

void mign_mffc_node_collect( mign_graph& mign, mign_node n, std::vector<mign_node>& support )
{
  if ( mign.bitmarks().is_marked( n ) ) return;
  mign.bitmarks().mark( n );

  if ( mign.get_ref( n ) > 0u || mign.is_input( n ) )
  {
    support.push_back( n );
    return;
  }

  for ( const auto& child : mign.children( n ) )
  {
    mign_mffc_node_collect( mign, child.node, support );
  }
}

void mign_mffc_mark_recurse( mign_graph& mign, mign_node curr, const std::vector<mign_node>& support )
{
  if ( mign.bitmarks().is_marked( curr ) ) return;
  mign.bitmarks().mark( curr );
  for (const auto& c : mign.children(curr))
  {
    if ( std::find( support.begin(), support.end(), c.node ) == support.end() )
    {
      mign_mffc_mark_recurse( mign, c.node, support );
    }
    else
    {
      /* last node to mark on this path */
      mign.bitmarks().mark( c.node );
    }
  }
}

/******************************************************************************
 * Public functions                                                           *
 ******************************************************************************/

unsigned mign_compute_mffc( mign_graph& mign, mign_node n, std::vector<mign_node>& support )
{
  assert( !mign.is_input( n ) );

  mign.init_refs();
  mign.inc_output_refs();
  mign.bitmarks().init_marks( mign.size() );

  const auto size1 = mign_mffc_node_deref( mign, n );
  for ( const auto& child : mign.children( n ) )
  {
    mign_mffc_node_collect( mign, child.node, support );
  }
  const auto size2 = mign_mffc_node_ref( mign, n );

  assert( size1 == size2 );

  return size1;
}

std::map<mign_node, std::vector<mign_node>> mign_mffcs( mign_graph& mign )
{
  std::map<mign_node, std::vector<mign_node>> map;

  auto nodes = mign_output_deque( mign );

  while ( !nodes.empty() )
  {
    const auto n = nodes.front();
    nodes.pop_front();

    if ( map.find( n ) != map.end() || mign.is_input( n ) ) { continue; }

    std::vector<mign_node> support;
    mign_compute_mffc( mign, n, support );
    map[n] = support;

    for ( auto s : support )
    {
      nodes.push_back( s );
    }
  }

  return map;
}

unsigned mign_mffc_size( const mign_graph& mign, mign_node n, const std::vector<mign_node>& support )
{
  std::map<mign_node, boost::default_color_type> colors;
  auto size = 0u;

  boost::depth_first_visit( mign.graph(), n,
                            mign_size_visitor( mign, support, size ),
                            boost::make_assoc_property_map( colors ),
                            [&support]( const mign_node& node, const mign_graph::graph_t& g ) { return boost::find( support, node ) != support.end(); } );

  return size;
}

std::vector<mign_node> mign_mffc_cone( const mign_graph& mign, mign_node n, const std::vector<mign_node>& support )
{
  std::map<mign_node, boost::default_color_type> colors;
  std::vector<mign_node> cone;

  boost::depth_first_visit( mign.graph(), n,
                            mign_cone_visitor( mign, support, cone ),
                            boost::make_assoc_property_map( colors ),
                            [&support]( const mign_node& node, const mign_graph::graph_t& g ) { return boost::find( support, node ) != support.end(); } );

  return cone;
}

bool mign_mffc_contains( const mign_graph& mign, mign_node root, const std::vector<mign_node>& support, mign_node curr )
{
  if ( root == curr ) return true;
  for ( const auto& c : mign.children(root) )
  {
    if (std::find(support.begin(), support.end(), c.node) != support.end())
    {
      continue;
    }
    if (mign_mffc_contains(mign, c.node, support, curr))
    {
      return true;
    }
  }
  return false;
}

unsigned mign_mffc_tipsize( const mign_graph& mign, const std::map<mign_node, std::vector<mign_node>> mffcs, mign_node curr )
{
  for ( const auto& mffc : mffcs )
  {
    if ( mign_mffc_contains(mign, mffc.first, mffc.second, curr ) )
    {
      return mign_mffc_size(mign, curr, mffc.second );
    }
  }
  return 0u; /* not contained in any mffc */
}

void mign_mffc_mark( mign_graph& mign, mign_node root, const std::vector<mign_node>& support )
{
  mign.bitmarks().init_marks( mign.size() );
  mign_mffc_mark_recurse( mign, root, support );
}

}

// Local Variables:
// c-basic-offset: 2
// eval: (c-set-offset 'substatement-open 0)
// eval: (c-set-offset 'innamespace 0)
// End:
