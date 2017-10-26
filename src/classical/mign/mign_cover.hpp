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
 * @file mign_cover.hpp
 *
 * @brief Store MIGN covers
 *
 * @author Eleonora Testa 
 * @since  2.3
 */

#ifndef MIGN_COVER_HPP
#define MIGN_COVER_HPP

#include <string>
#include <vector>

#include <boost/range/iterator_range.hpp>

#include <classical/mign/mign.hpp>
#include <classical/mign/mign_cuts_paged.hpp>

namespace cirkit
{

class mign_cover
{
public:
  using index_range = boost::iterator_range<std::vector<unsigned>::const_iterator>;

  mign_cover( unsigned cut_size, const mign_graph& mign );
  mign_cover(); 

  void add_cut( mign_node n, const mign_cuts_paged::cut& cut , const unsigned T, const std::vector<unsigned> w, const std::vector<bool> neg_un);
  void add_cut_maj( mign_node n, const mign_cuts_paged::cut& cut, const boost::dynamic_bitset<> tt_maj, const unsigned tt_maj_and_or, tt reminder, std::string exact); 
  bool has_cut( mign_node n ) const;
  unsigned has_threshold (mign_node n) const; 
  std::vector<unsigned> has_weights (mign_node n) const; 
  std::vector<bool> neg_una(mign_node n) const; 
  index_range cut( mign_node n ) const;
  boost::dynamic_bitset<> has_tt_maj(mign_node n) const; 
  unsigned has_tt_almostmaj_or_maj(mign_node n) const;
  tt has_tt_reminder(mign_node n) const;
  std::string has_tt_exact (mign_node n) const; 

  inline unsigned cut_size() const { return _cut_size; }
  inline unsigned lut_count() const { return count; }
  
  
private:
  unsigned                              _cut_size; /* remember cut_size */

  std::vector<unsigned>                 threshold; // T = 0 means no "big" threshold. Normal majority node. 
  std::vector<std::vector<unsigned>>    weights; 
  std::vector<std::vector<bool>>        neg_un; 
  std::vector<boost::dynamic_bitset<>>  tt_is_maj;
  std::vector<unsigned>                 tt_is_almostmaj_or_maj;
  std::vector<tt>                       tt_is_reminder;
  std::vector<std::string>              tt_is_exact; 
  std::vector<unsigned>                 offset; /* address from node index to leafs, 0 if unused */
  std::vector<unsigned>                 leafs;  /* first element is unused, then | #leafs | l_1 | l_2 | ... | l_k | */
  unsigned                              count = 0u;
};

 mign_graph mign_cover_write (const mign_graph& mign); 
 mign_graph mign_cover_write_maj (const mign_graph& mign); 
 std::vector<mign_graph> mign_cover_write_multi (mign_graph& mign); 
 
}

#endif

// Local Variables:
// c-basic-offset: 2
// eval: (c-set-offset 'substatement-open 0)
// eval: (c-set-offset 'innamespace 0)
// End:
