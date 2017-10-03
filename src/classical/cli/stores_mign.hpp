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
 * @file stores_experimental.hpp
 *
 * @brief Meta-data for stores
 *
 * @author Mathias Soeken
 * @since  2.3
 */

#ifndef MIGN_STORES_HPP
#define MIGN_STORES_HPP

#include <alice/command.hpp>

#include <string>
#include <vector>

#include <cli/stores.hpp>
#include <core/utils/bdd_utils.hpp>

#include <classical/aig.hpp>
#include <classical/mig/mig.hpp>
#include <classical/mign/mign.hpp>
#include <classical/utils/expression_parser.hpp>
#include <classical/mig/mig_utils.hpp>

namespace alice
{

using namespace cirkit;

/******************************************************************************
 * mign_graph                                                                 *
 ******************************************************************************/

template<>
struct store_info<mign_graph>
{
  static constexpr const char* key         = "migns";
  static constexpr const char* option      = "mign";
  static constexpr const char* mnemonic    = "M";
  static constexpr const char* name        = "MIG-n";
  static constexpr const char* name_plural = "MIg-n's";
};

template<>
std::string store_entry_to_string<mign_graph>( const mign_graph& mign );

template<>
inline bool store_can_convert<mig_graph, mign_graph>() { return true; }

template<>
mign_graph store_convert<mig_graph, mign_graph>( const mig_graph& mig );

template<>
inline bool store_can_read_io_type<mign_graph, io_verilog_tag_t>( command& cmd ) { return true; }

template<>
mign_graph store_read_io_type<mign_graph, io_verilog_tag_t>( const std::string& filename, const command& cmd );

template<>
inline bool store_can_write_io_type<mign_graph, io_verilog_tag_t>( command& cmd ) { return true; }

template<>
void store_write_io_type<mign_graph, io_verilog_tag_t>( const mign_graph& mign, const std::string& filename, const command& cmd );
    
template<>
void print_store_entry_statistics<mign_graph>( std::ostream& os, const mign_graph& mign );

template<>
inline bool store_can_convert<mign_graph, mig_graph>() { return true; }

template<>
mig_graph store_convert<mign_graph, mig_graph>( const mign_graph& mign );

/*template<>
inline bool store_can_convert<bdd_function_t, mign_graph>() { return true; }

template<>
mign_graph store_convert<bdd_function_t, mign_graph>( const bdd_function_t& bdd );*/


/*template<>
inline bool store_can_write_io_type<mign_graph, io_blif_tag_t>( command& cmd ) { return true; }

template<>
void store_write_io_type<mign_graph, io_blif_tag_t>( const mign_graph& mign, const std::string& filename, const command& cmd );*/
    

}

#endif

// Local Variables:
// c-basic-offset: 2
// eval: (c-set-offset 'substatement-open 0)
// eval: (c-set-offset 'innamespace 0)
// End:
