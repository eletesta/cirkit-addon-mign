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

/**
 * @file mign_bitmarks.hpp
 *
 * @brief Store MIGN bitmarks
 *
 * @author Heinz Riener
 * @since  2.3
 */

#ifndef MIGN_BITMARKS_HPP
#define MIGN_BITMARKS_HPP

#include <classical/mign/mign.hpp>

#include <boost/dynamic_bitset.hpp>

#include <vector>

namespace cirkit
{

class mign_bitmarks
{
public:
  void init_marks( unsigned size, unsigned num_colors = 1u );
  bool is_marked( mign_node n, unsigned color = 0u ) const;
  void mark( mign_node n, unsigned color = 0u );
  void resize_marks( mign_node n );
  void unmark( mign_node n, unsigned color = 0u );
  void invert( unsigned color = 0u );
  void reset( unsigned color = 0u );
  unsigned count( unsigned color = 0u ) const;
  unsigned size() const;
  void move( unsigned dst, unsigned src );
  
  unsigned alloc();
  void free( unsigned color );

  unsigned num_layers() const;
  unsigned num_used_layers() const;

  boost::dynamic_bitset<> get( unsigned color = 0u ) const;

  boost::dynamic_bitset<> get_used() const;
  void set_used( const boost::dynamic_bitset<>& used_ );
  
private:
  std::vector<boost::dynamic_bitset<>> marks;
  boost::dynamic_bitset<> used;
}; /* mign_bitmarks */


void mign_mark_inner_nodes( mign_graph& mign, const std::vector<mign_node>& nodes, unsigned color );

}

#endif

// Local Variables:
// c-basic-offset: 2
// eval: (c-set-offset 'substatement-open 0)
// eval: (c-set-offset 'innamespace 0)
// End:
