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
 * @file xmg_io.hpp
 *
 * @brief I/O routines for MIGN _ Verilog Files and Write for BLIF 
 *
 * @author Eleonora Testa 
 * @since  2.3
 */

#ifndef MIGN_IO_HPP
#define MIGN_IO_HPP

#include <string>

#include <core/properties.hpp>
#include <classical/mign/mign.hpp>

namespace cirkit
{

std::string get_operand_name( const mign_graph& mign, const mign_function& f );
mign_function find_function( const std::unordered_map<std::string, mign_function>& name_to_function, const std::string& name ); 
mign_graph read_verilog_compact( const std::string& filename );

void write_verilog_compact( const mign_graph& mign, const std::string& filename, const properties::ptr& settings = properties::ptr() );
void write_verilog_perm ( const mign_graph& mign); /*const std::string& filename, const properties::ptr& settings = properties::ptr()*/ 

void write_blif ( const mign_graph& mign,const std::string& filename , const properties::ptr& settings);
}

#endif

// Local Variables:
// c-basic-offset: 2
// eval: (c-set-offset 'substatement-open 0)
// eval: (c-set-offset 'innamespace 0)
// End:
