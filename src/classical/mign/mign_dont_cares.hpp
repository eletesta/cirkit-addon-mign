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
 * @file xmg_dont_cares.hpp
 *
 * @brief Compute various don't cares
 *
 * @author Mathias Soeken
 * @since  2.4
 */

#ifndef MIGN_DONT_CARES_HPP
#define MIGN_DONT_CARES_HPP

#include <boost/dynamic_bitset.hpp>

#include <classical/mign/mign.hpp>
#include <classical/xmg/xmg.hpp>

namespace cirkit
{

/* checks whether `pattern` is a don't care at `node`, meaning that by inverting node, one does not see a difference at the output -- uses XMGs -- only works with mign made of 3 inputs */
bool mign_is_observable_at_node( const mign_graph& mign, mign_node node, const boost::dynamic_bitset<>& pattern , const std::map<xmg_node, bool>& assignment);

}

#endif

// Local Variables:
// c-basic-offset: 2
// eval: (c-set-offset 'substatement-open 0)
// eval: (c-set-offset 'innamespace 0)
// End:
