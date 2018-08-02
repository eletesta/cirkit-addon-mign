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

#include "mign_bitmarks.hpp"

#include <classical/mign/mign.hpp>

namespace cirkit
{

/******************************************************************************
 * mign_bitmarks                                                               *
 ******************************************************************************/

void mign_bitmarks::init_marks( unsigned size, unsigned num_colors )
{
  assert( num_colors > 0u );

  marks.resize( num_colors );
  for ( auto& layer : marks )
  {
    layer.resize( size );
    layer.reset();
  }

  used.resize( num_colors );
  used = ~used;
}

bool mign_bitmarks::is_marked( mign_node n, unsigned color ) const
{
  assert( marks.size() > color );
  return n < marks[color].size() && marks[color][n];
}

void mign_bitmarks::mark( mign_node n, unsigned color )
{
  assert( marks.size() > color );
  if ( n < marks[color].size() )
  {
    marks[color].set( n );
  }
}

void mign_bitmarks::resize_marks( mign_node n )
{
  assert( marks.size() > 0u );
  if ( n >= marks[0u].size() )
  {
    for ( auto i = 0u; i < num_layers(); ++i )
    {
      marks[i].resize( n+1 );
    }
  }
}

void mign_bitmarks::unmark( mign_node n, unsigned color )
{
  assert( marks.size() > color );
  if ( n < marks[color].size() )
  {
    marks[color].reset( n );
  }
}

void mign_bitmarks::invert( unsigned color )
{
  assert( marks.size() > color );
  marks[color] = ~marks[color];
}

unsigned mign_bitmarks::alloc()
{
  if ( used.count() == used.size() )
  {
    const auto nlayers = num_layers();
    marks.resize( nlayers + 1u );
    used.resize( nlayers + 1u );
    used.set( nlayers );
    return nlayers;
  }

  for ( auto i = 0u; i < used.size(); ++i )
  {
    if ( !used[i] )
    {
      used.set(i);
      return i;
    }
  }

  assert( false && "unreachable" );
  return 1; 
}

void mign_bitmarks::free( unsigned color )
{
  assert( marks.size() > color );
  assert( used[color] );
  marks[color].reset();
  used.reset(color);
}

boost::dynamic_bitset<> mign_bitmarks::get_used() const
{
  return used;
}

void mign_bitmarks::set_used( const boost::dynamic_bitset<>& used_ )
{
  used = used_;
}

unsigned mign_bitmarks::num_layers() const
{
  return marks.size();
}

unsigned mign_bitmarks::num_used_layers() const
{
  return used.count();
}

void mign_bitmarks::reset( unsigned color )
{
  assert( marks.size() > color );
  marks[color].reset();
}

unsigned mign_bitmarks::count( unsigned color ) const
{
  assert( marks.size() > color );
  return marks[color].count();
}

unsigned mign_bitmarks::size() const
{
  /* all layers are supposed to have the same size */
  assert( marks.size() > 0u );
  return marks[0u].size();
}

void mign_bitmarks::move( unsigned dst, unsigned src )
{
  assert( marks.size() > dst );
  assert( marks.size() > src );
  assert( marks[dst].size() == marks[src].size() );
  marks[dst] |= marks[src];
}

boost::dynamic_bitset<> mign_bitmarks::get( unsigned color ) const
{
  return marks[color];
}

/******************************************************************************
 * private functions                                                          *
 ******************************************************************************/

/******************************************************************************
 * public functions                                                           *
 ******************************************************************************/

void mign_mark_inner_nodes( mign_graph& mign, const std::vector<mign_node>& nodes, unsigned color )
{
  for ( const auto& node : nodes )
  {
    if ( mign.is_input(node) || mign.fanout_count(node) == 0u ) continue;
    mign.bitmarks().mark( node, color );
  }
}

}

// Local Variables:
// c-basic-offset: 2
// eval: (c-set-offset 'substatement-open 0)
// eval: (c-set-offset 'innamespace 0)
// End:
