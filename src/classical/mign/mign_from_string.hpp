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
 * @file mig_from_string.hpp
 *
 * @brief Creates an MIG from a string
 *
 * @author Mathias Soeken
 * @since  2.3
 */

#ifndef MIGN_FROM_STRING_HPP
#define MIGN_FROM_STRING_HPP

#include <string>
#include <vector>

#include <core/properties.hpp>
#include <classical/mign/mign.hpp>

namespace cirkit
{

mign_graph mign_from_string( const std::string& expr,
                           const properties::ptr& settings = properties::ptr(),
                           const properties::ptr& statistics = properties::ptr() );

mign_function mign_from_string( mign_graph& mig, const std::string& expr,
                              const properties::ptr& settings = properties::ptr(),
                              const properties::ptr& statistics = properties::ptr() );

std::string mign_to_string( const mign_graph& mig, const mign_function& f,
                           const properties::ptr& settings = properties::ptr(),
                           const properties::ptr& statistics = properties::ptr() );

}

#endif

// Local Variables:
// c-basic-offset: 2
// eval: (c-set-offset 'substatement-open 0)
// eval: (c-set-offset 'innamespace 0)
// End:
