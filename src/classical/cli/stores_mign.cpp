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

#include "stores_mign.hpp"

#include <sstream>

#include <boost/format.hpp>

#include <core/properties.hpp>
#include <core/graph/depth.hpp>
#include <core/utils/range_utils.hpp>

#include <classical/io/write_bench.hpp>
#include <classical/mign/mign_io.hpp>
#include <classical/mign/mign_utils.hpp>
#include <classical/mign/mig_to_mign.hpp>
#include <classical/mign/mign_lut_based_synthesis.hpp> 

namespace alice
{

using namespace cirkit;

/******************************************************************************
 * mign_graph                                                                 *
 ******************************************************************************/

template<>
std::string store_entry_to_string<mign_graph>( const mign_graph& mign )
{
  const auto name = mign.name();
  return boost::str( boost::format( "%s i/o = %d/%d" ) % ( name.empty() ? "(unnamed)" : name ) % mign.inputs().size() % mign.outputs().size() );
}

template<>
mign_graph store_convert<mig_graph, mign_graph>( const mig_graph& mig )
{
  return mig_to_mign( mig );
}

template<>
mig_graph store_convert<mign_graph, mig_graph>( const mign_graph& mign )
{
  return mign_to_mig( mign );
}

template<> 
void print_store_entry_statistics<mign_graph>( std::ostream& os, const mign_graph& mign )
{
  mign_print_stats( mign, os );
}

template<>
mign_graph store_read_io_type<mign_graph, io_verilog_tag_t>( const std::string& filename, const command& cmd )
{
  return read_verilog_compact( filename );
}

template<>
void store_write_io_type<mign_graph, io_verilog_tag_t>( const mign_graph& mign, const std::string& filename, const command& cmd )
{
  write_verilog_compact( mign, filename );
}


/*template<>
void store_write_io_type<mign_graph, io_blif_tag_t>( const mign_graph& mign, const std::string& filename, const command& cmd )
{
	auto settings = std::make_shared<properties>();
    write_blif( mign, filename, settings );
}*/

}
// Local Variables:
// c-basic-offset: 2
// eval: (c-set-offset 'substatement-open 0)
// eval: (c-set-offset 'innamespace 0)
// End:
