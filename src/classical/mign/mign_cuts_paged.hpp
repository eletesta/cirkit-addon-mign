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
 * @file mign_cuts_paged.hpp
 *
 * @brief MiGN cut enumeration
 *
 * @author Eleonora Testa 
 * @since  2.3
 */

#ifndef MIGN_CUTS_PAGED_HPP
#define MIGN_CUTS_PAGED_HPP

#include <map>
#include <vector>

#include <boost/dynamic_bitset.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <boost/range/iterator_range.hpp>
#include <boost/range/combine.hpp>
#include <boost/iterator/zip_iterator.hpp>

#include <core/properties.hpp>
#include <core/utils/paged_memory.hpp>
#include <classical/mign/mign.hpp>
#include <classical/utils/truth_table_utils.hpp>

namespace cirkit
{


class mign_cuts_paged final
{
public:
  using cut = paged_memory::set;
  using cone = paged_memory::set;

  mign_cuts_paged( mign_graph& mign, unsigned k, const properties::ptr& settings = properties::ptr() );
  mign_cuts_paged( mign_graph& mign, unsigned k,
                  const std::vector<mign_node>& start,
                  const std::vector<mign_node>& boundary,
                  const std::vector<std::pair<unsigned, unsigned>>& levels,
                  const properties::ptr& settings = properties::ptr() );

  const mign_graph& mign() const;

  unsigned total_cut_count() const;
  double enumeration_time() const;

  unsigned memory() const;
  unsigned count( mign_node node ) const;
  boost::iterator_range<paged_memory::iterator> cuts( mign_node node );
  boost::iterator_range<paged_memory::iterator> cut_cones( mign_node node );

  tt simulate( mign_node node, const cut& c ) const;
  tt simulate_bis( mign_node node, const cut& c ) const;
  unsigned depth( mign_node node, const cut& c ) const;
  unsigned size( mign_node node, const cut& c ) const;

  unsigned index( const cut& c ) const;
  cut      from_address( unsigned address );

  void foreach_cut( const std::function<void(mign_node, cut&)>& func );

private:
  void enumerate();
  //void enumerate_with_xor_blocks( const std::unordered_map<mign_node, mign_xor_block_t>& blocks );
  //void enumerate_partial( const std::vector<mign_node>& start, const std::vector<mign_node>& boundary );

  void enumerate_node_with_bitsets( mign_node n, const std::vector<mign_node>& ns );

  using local_cut_vec_t = std::vector<std::tuple<boost::dynamic_bitset<>, unsigned, boost::dynamic_bitset<>>>;
  local_cut_vec_t enumerate_local_cuts( mign_node n1, mign_node n2, unsigned max_cut_size );
  local_cut_vec_t enumerate_local_cuts( mign_node n1, mign_node n2, mign_node n3, unsigned max_cut_size );
  local_cut_vec_t enumerate_local_cuts( mign_node n1, mign_node n2, mign_node n3, mign_node n4, mign_node n5, unsigned max_cut_size );
  local_cut_vec_t enumerate_local_cuts( const std::vector<mign_node>& ns, unsigned max_cut_size ); 
  void merge_cut( local_cut_vec_t& local_cuts, const boost::dynamic_bitset<>& new_cut, unsigned min_level, const boost::dynamic_bitset<>& new_cone ) const;

  std::vector<unsigned> get_extra( unsigned depth, unsigned size ) const;

private:
  const mign_graph& _mign;
  unsigned         _k;
  unsigned         _priority = 10u;
  unsigned         _extra    = 0u;
  paged_memory     data;
  paged_memory     cones;

  double           _enumeration_time = 0.0;

  unsigned         _top_index = 0u; /* index when doing topo traversal */

  std::vector<std::pair<unsigned, unsigned>> _levels;
};

}

#endif

// Local Variables:
// c-basic-offset: 2
// eval: (c-set-offset 'substatement-open 0)
// eval: (c-set-offset 'innamespace 0)
// End:
