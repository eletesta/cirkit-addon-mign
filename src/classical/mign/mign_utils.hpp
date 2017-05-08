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
 * @brief I/O routines for MIGN _ Verilog Files 
 *
 * @author Eleonora Testa 
 * @since  2.3
 */

#ifndef MIGN_UTILS_HPP
#define MIGN_UTILS_HPP

#include <string>
#include <ostream>

#include <core/properties.hpp>
#include <classical/mign/mign.hpp>
#include <classical/utils/truth_table_utils.hpp>

namespace cirkit
{
mign_function create_mign_function(mign_node n, bool complemented); 
std::vector<unsigned> mign_compute_levels( const mign_graph& mign );
unsigned evaluate_depth (const mign_graph& mign);
float evaluate_depth_th (const mign_graph& mign);
float evaluate_energy (mign_graph& mign); 
float leakage_energy (mign_graph& mign, float depth); 
void mign_print_stats(const mign_graph& mign, std::ostream& os);
 
std::deque<mign_node> mign_output_deque( const mign_graph& mign );
tt mign_simulate_cut( const mign_graph& mign, mign_node root, const std::vector<mign_node>& leafs ); 
tt mign_simulate_cut_func( const mign_graph& mign, mign_node root, const std::vector<mign_node>& leafs );

std::vector<mign_graph> simplify_vector_mign(std::vector<mign_graph>& mign); 
unsigned input_threshold (unsigned tt_size, unsigned T, unsigned polarity, std::vector<unsigned> weigths); 

void analysis (mign_graph& mign);
}

#endif

// Local Variables:
// c-basic-offset: 2
// eval: (c-set-offset 'substatement-open 0)
// eval: (c-set-offset 'innamespace 0)
// End:
