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
 * @file mign_mffc.hpp
 *
 * @brief Compute MFFCs
 *
 * @author Mathias Soeken
 * @since  2.3
 */

#ifndef MIGN_MFFC_HPP
#define MIGN_MFFC_HPP

#include <map>
#include <vector>

#include <classical/mign/mign.hpp>

namespace cirkit
{

unsigned mign_compute_mffc( mign_graph& mign, mign_node n, std::vector<mign_node>& support );
std::map<mign_node, std::vector<mign_node>> mign_mffcs( mign_graph& mign );

/* counts nodes including the root, but excluding the leafs */
unsigned mign_mffc_size( const mign_graph& mign, mign_node n, const std::vector<mign_node>& support );
/* returns nodes including the root, but excluding the leafs */
std::vector<mign_node> mign_mffc_cone( const mign_graph& mign, mign_node n, const std::vector<mign_node>& support );

/* check if curr is contained in the mffc defined by the pair (root,support) */
bool mign_mffc_contains( const mign_graph& mign, mign_node root, const std::vector<mign_node>& support, mign_node curr );
/* compute the number of nodes of the subcone (within mffcs) rooted by curr */
unsigned mign_mffc_tipsize( const mign_graph& mign, const std::map<mign_node, std::vector<mign_node>> mffcs, mign_node curr );

/* mark all nodes of a mffc defined by the pair (root,support) on mign */
void mign_mffc_mark( mign_graph& mign, mign_node root, const std::vector<mign_node>& support );

}

#endif

// Local Variables:
// c-basic-offset: 2
// eval: (c-set-offset 'substatement-open 0)
// eval: (c-set-offset 'innamespace 0)
// End:
