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

#include "mign_simulate.hpp"

#include <boost/range/algorithm_ext/iota.hpp>

#include <classical/mign/math_utils.hpp> 

#include <core/utils/combinations.hpp>

namespace cirkit
{

/******************************************************************************
 * Truth table simulation                                                     *
 ******************************************************************************/

tt mign_tt_simulator::get_input( const mign_node& node, const mign_graph& mign ) const
{
  return tt_nth_var( mign.input_index( node ) );
}

tt mign_tt_simulator::get_constant() const
{ 
  return tt_const0();
}

tt mign_tt_simulator::invert( const tt& v ) const
{
  return ~v;
}


tt mign_tt_simulator::maj_op( const mign_node& node, const std::vector<tt>& v ) const
{
	
	//std::cout << "simulation of majority node " << node << " che ha come figli: "<< std::endl;  
  const auto n = v.size(); 
  //std::cout << " simulate con figli = " << n << std::endl; 
  for ( auto x = 0; x <n; ++x)
  {
	  
  }
  std::vector<unsigned> numbers(n);
  boost::iota( numbers, 0u );

  auto fact_size = fact(n); 
  auto fact_nk = fact(n-2); 
  auto number_c = fact_size/(2*fact_nk); 
  std::vector<std::vector<unsigned>> combinations(number_c);
  
  auto count = 0, counting = 0; 
  
  std::vector<tt> _v;
  
  //std::cout << " after each combination " << std::endl; 
  for ( auto t = 0; t < v.size(); ++t)
  {
	   _v.push_back(v[t]); 
	  //if (tt_num_vars(v[t]) < n)
		 // tt_extend(_v[t],n); 	
	 
  }
  
  auto value = 0; 
  for (auto x = 0; x < v.size(); ++x)
  {
	  value = std::max<int>( value, tt_num_vars( _v[x] ) );
 
  }
  
  for (auto x = 0; x < v.size(); ++x)
  {
	  //if (tt_num_vars(_v[x]) < value)
		  tt_extend(_v[x], value); 
	  //std::cout << " tt = " << _v[x] << std::endl; 
  }
  
   
   //for (auto & comb:combinations)
  //{
	//  assert(comb.size()== 2); 
	//  tt_align(_v[comb[0]],_v[comb[1]]);
  //}
  
  
   auto pot = _v[0].size();
  tt result(pot);  
  for ( auto x = 0; x < pot; ++x)
  {
	  auto num_ones = 0; 
	  for ( auto y = 0; y <_v.size(); ++y)
	  {
		  if (_v[y][x] == 1)
			  ++num_ones; 
		  else 
			  num_ones = num_ones; 
	  }
	  result[x] = (num_ones > _v.size() - num_ones); 
		 
  }

  return result;
}

}

// Local Variables:
// c-basic-offset: 2
// eval: (c-set-offset 'substatement-open 0)
// eval: (c-set-offset 'innamespace 0)
// End:
