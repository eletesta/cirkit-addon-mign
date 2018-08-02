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

#include "mign_utils.hpp"

#include <fstream>
#include <unordered_map>
#include <vector>

#include <boost/format.hpp>
#include <classical/mign/mign_simulate.hpp>
#include <classical/mign/math_utils.hpp>

using boost::format;

namespace cirkit
{

/******************************************************************************
 * Types                                                                      *
 ******************************************************************************/

tt tt_const0_leaf( unsigned leaf_size)
{
	auto c = int_pow2(leaf_size); 
	return boost::dynamic_bitset<>( c, 0u );
}

tt tt_const1_leaf(unsigned leaf_size)
{
	  return ~tt_const0_leaf(leaf_size);
}

/******************************************************************************
 * Private functions                                                          *
 ******************************************************************************/
	
int compute_depth_rec ( const mign_graph& mign, mign_node node, std::vector<int>& depths)
{

	if ( depths[node] >= 0 ) return depths[node];
	
	  if ( mign.is_input( node ) )
	  {
	    return depths[node] = 0;
	  }
	  else
	  {
	    int depth = 0;
	    for ( const auto& c : mign.children( node ) )
	    {
	      depth = std::max( compute_depth_rec( mign, c.node, depths ), depth );
	    }
	    return depths[node] = ( depth + 1 );
	  }		
}

int compute_depth_inv_rec ( const mign_graph& mign, mign_node node, std::vector<int>& depths)
{

	if ( depths[node] >= 0 ) return depths[node];
	
	  if ( mign.is_input( node ) )
	  {
	    return depths[node] = 0;
	  }
	  else
	  {
	    int depth = 0;
	    for ( const auto& c : mign.children( node ) )
	    {
		  auto d = compute_depth_inv_rec( mign, c.node, depths ); 
  		  if ((c.node != 0) && (c.complemented == 1))
  		  {
  			  d++; 
  		  }
	      depth = std::max( d, depth );
		  
	    }
		
	    return depths[node] = ( depth + 1 );
	  }		
}

float compute_depth_rec_th ( const mign_graph& mign, mign_node node, std::vector<float>& depths)
{

	if ( depths[node] >= 0 ) return depths[node];
	
	  if ( mign.is_input( node ) )
	  {
	    return depths[node] = 0;
	  }
	  else
	  {
	    float depth = 0;
	    
		const auto children = mign.children( node ); 
	    for ( const auto& c : children )
	    {
		  auto d = compute_depth_rec_th( mign, c.node, depths ); 
  		  if ((c.node != 0) && (c.complemented == 1))
  		  {
  			  auto c_c =  mign.children( c.node ); 
			  if (c_c.size() == 3 )
				  d = d + 7.92; 
			  else if (c_c.size() == 5 )
				  d = d + 8.27; 
			  else if (c_c.size() == 7 )
				  d = d + 8.56; 
			  else if (c_c.size() == 9 )
				  d = d + 8.81; 
			  else 
				  d++; 
  		  }
	      depth = std::max( d, depth );
	    }
		
		if ((children.size() == 3) || (children.size() > 9)) // to avoid problems 
	    {	
			depths[node] = ( depth + 29.52 );
		}
			
		else if (children.size() == 5)
		{
			
			 depths[node] = ( depth + 30.60 );
		}
	   
		else if (children.size() == 7)
	    {
			
			depths[node] = ( depth + 32.42 ); 	
		}
		
		else 
		{
			depths[node] = ( depth + 34.47 );
		}
		return depths[node]; 
	  }		
}


bool check_one (std::vector<mign_function> children)
{
	for (auto & c :children)
	{
		if ((c.node == 0) && (c.complemented == 1))
			return true; 
	}
	return false; 
}

unsigned compute_ce (const mign_graph& mign)
{

	unsigned ce = 0u; 
	for ( auto& node : mign.topological_nodes())
	{
		 //std::cout << " nodo " << node << std::endl; 
				
	     if ( mign.is_input( node) ) { continue; }

	    const auto c = mign.children( node);
	
	  	for (auto x :c)
	  	{
			//std::cout << " figlio = " << x.node << " con complemented = " << x.complemented << std::endl; 
	  		if (( x.node != 0) && (x.complemented == 1))
			{
				++ce; 
			}
	  	}	
	}
	
	for ( const auto& output : mign.outputs())
	{
		//std::cout << " output " << output.first.node << std::endl; 
		if ( output.first.complemented == 1)
			++ce; 
	}
	
	//std::cout << " ce = " << std::endl; 
	
	return ce; 
}

unsigned compute_dangling (const mign_graph& mign)
{

	unsigned dan = 0u; 
	for ( auto& node : mign.topological_nodes())
	{		
	     if ( mign.is_input( node) ) { continue; }

	    const auto c = mign.children( node);
	
	  	for (auto x :c)
	  	{
			//std::cout << " figlio = " << x.node << " con complemented = " << x.complemented << std::endl; 
	  		if ((mign.is_input( x.node)) || (x.node == 0))
			{
				++dan; 
			}
	  	}	
	}
		
	return dan; 
}

unsigned compute_inv (const mign_graph& mign)
{

	unsigned inv = 0u; 
	std::vector<mign_node> count; 
	
	for ( auto& node : mign.topological_nodes())
	{		
	    if ( mign.is_input( node) ) { continue; }

	    const auto c = mign.children( node);
	
	  	for (auto x :c)
	  	{

	  		if (( x.node != 0) && (x.complemented == 1))
			{
				
				if (std::find(count.begin(), count.end(), x.node) != count.end())
					continue; 
				else 
				{
					++inv;
					count.push_back(x.node); 
				}
				 
			}
	  	}	
	}
	
	for ( const auto& output : mign.outputs())
	{
		if (( output.first.complemented == 1) && ( output.first.node != 0))
		{
			if (std::find(count.begin(), count.end(), output.first.node) != count.end())
				continue; 
			else 
			{
				inv++; 
				count.push_back(output.first.node); 
			}
			
		}
		
	}
		
	return inv; 
}


unsigned compute_larger (const mign_graph& mign)
{
	unsigned larger = 0u; 
	for ( auto& node : mign.topological_nodes())
	{
	     if ( mign.is_input( node) ) { continue; }

	      const auto c = mign.children( node);
		  if (c.size() > larger)
			  larger = c.size(); 
	} 
	return larger; 
}

unsigned compute_smaller (const mign_graph& mign)
{
	unsigned smaller = 101u; 
	for ( auto& node : mign.topological_nodes())
	{
	     if ( mign.is_input( node) ) { continue; }

	      const auto c = mign.children( node);
		  if (c.size() < smaller)
			  smaller = c.size(); 
	} 
	return smaller; 
}

void print_big (const mign_graph& mign, std::ostream& os = std::cout)
{
	std::vector<std::pair<unsigned,unsigned>> store; 
	auto flag = 0u; 
	for ( auto& v : mign.topological_nodes())
	{
	flag = 0; 
     if ( mign.is_input( v) ) { continue; }
	 
	 const auto c = mign.children(v);
	 for ( auto& curr : store )
		 {
			 if (c.size() == curr.first)
			 {
				 flag = 1; 
				 ++curr.second; 
			 }
			 else continue; 
		 } 
		 if (flag == 0)
		 {
			 store.push_back({c.size(),1}); 
		 }
	}
	
	for (auto& node : store)
	{
		os << " There are " << node.second << " MAJs with " << node.first << " inputs" << std::endl; 
	}
}

void print_depth_statistics(mign_graph mign)
{
	
    std::vector<int> depths (mign.size(), -1);
    int depth = -1;
    
    for ( const auto& i : mign.outputs())
    {
        depth = std::max (compute_depth_rec(mign, i.first.node,depths ), depth);
    }
	
	std::vector<std::pair<int,unsigned>> count_depth; 
	
	for (auto x = 0; x <depths.size(); ++x)
	{
		//std::cout << " depth del nodo " << x << " = "<< depths[x] << std::endl; 
		auto flag = 1; 
		for (auto & pair : count_depth)
		{
			if (depths[x] == pair.first)
				{
					flag = 0; 
					++pair.second;  
				}
		}
		if (flag == 1)
			count_depth.push_back({depths[x],1});
		
	}
	
	for (auto & pair : count_depth)
	{
		//std::cout << " There are " << pair.second << " nodes in level" << pair.first << std::endl; 
	}
}
/******************************************************************************
 * Public functions                                                           *
 ******************************************************************************/

mign_function create_mign_function(mign_node n, bool complemented)
{
	mign_function f; 
	f.node = n; 
	f.complemented = complemented; 
	return f; 
}

float evaluate_energy (mign_graph& mign)
{
	
	float energy = 0; 
	mign.compute_fanout(); 
	for ( auto& node : mign.topological_nodes())
	{
	     if ( mign.is_input( node) ) { continue; }

	      const auto c = mign.children( node);
		  if ((c.size() <= 3 ) || (c.size() > 9))
			  energy = energy + 0.59 + 0.013*mign.fanout_count(node); 
		  else if (c.size() == 5 ) 
			  energy = energy + 0.66 + 0.013*mign.fanout_count(node);
		 else  if (c.size() == 7 ) 
			  energy = energy + 0.71 + 0.013*mign.fanout_count(node);
		  else if (c.size() == 9 ) 
			  energy = energy + 0.77 + 0.013*mign.fanout_count(node);
	}
	
	return energy; 
}

float leakage_energy (mign_graph& mign, float depth)
{
	
	float energy = 0; 
	mign.compute_fanout(); 
	for ( auto& node : mign.topological_nodes())
	{
	     if ( mign.is_input( node) ) { continue; }

	      const auto c = mign.children( node);
		  if ((c.size() <= 3 ) || (c.size() > 9))
			  energy = energy + 0.6*0.28*depth; 
		  else if (c.size() == 5 ) 
			  energy = energy + 0.6*0.26*depth; 
		 else  if (c.size() == 7 ) 
			  energy = energy + 0.6*0.26*depth; 
		  else if (c.size() == 9 ) 
			  energy = energy + 0.6*0.28*depth; 
	}
	
	return energy; 
}

std::vector<unsigned> mign_compute_levels( const mign_graph& mign )
{
  std::vector<unsigned> levels( mign.size() );
  for ( const auto& n : mign.topological_nodes() )
  {
    if ( mign.is_input( n ) )
    {
      levels[n] = 0u;
    }
    else
    {
      auto level = 0u;
      for ( const auto& c : mign.children( n ) )
      {
        level = std::max( level, levels[c.node] );
      }
      levels[n] = level + 1u;
    }
  }
  return levels;
}

 unsigned evaluate_depth (const mign_graph& mign)
   {
        std::vector<int> depths (mign.size(), -1);
        int depth = -1;
        
        for ( const auto& i : mign.outputs())
        {
            depth = std::max (compute_depth_rec(mign, i.first.node,depths ), depth);
        }
		
        return static_cast<unsigned>( depth );
    }
	
    unsigned evaluate_depth_inv (const mign_graph& mign)
      {
           std::vector<int> depths (mign.size(), -1);
           int depth = -1;
        
           for ( const auto& i : mign.outputs())
           {
			   auto d = compute_depth_inv_rec(mign, i.first.node,depths ); 
			   if (i.first.complemented == 1)
				   d++;
               depth = std::max (d, depth);
			   
           }
		
           return static_cast<unsigned>( depth );
       }
	
float evaluate_depth_th (const mign_graph& mign)
      {
           std::vector<float> depths (mign.size(), -1);
           float depth = -1;
        
           for ( const auto& i : mign.outputs())
           {
               depth = std::max (compute_depth_rec_th(mign, i.first.node,depths ), depth);
           }
		
           return depth;
       }
	
std::deque<mign_node> mign_output_deque( const mign_graph& mign )
{
	 
	 std::deque<mign_node> deque;

	  for ( const auto& output : mign.outputs() )
	  {
	    deque.push_back( output.first.node );
	  }

	  return deque;
}


tt mign_simulate_cut( const mign_graph& mign, mign_node root, const std::vector<mign_node>& leafs )
{
	  std::map<mign_node, tt> inputs;
	  auto size_l = leafs.size(); 
	  
	  auto i = 0u;
	 
	  for ( auto child : leafs )
	  {
		  if (child == 0)
		  {
			  size_l--; 
			  inputs[child] = tt_const0(); 
		  }
		  else
	    inputs[child] = tt_nth_var( i++ );
	  }

	  mign_tt_simulator tt_sim;
	  mign_partial_node_assignment_simulator<tt> sim( tt_sim, inputs, tt_const0() );
	  
	  auto tt = simulate_mign_function( mign, root, sim );

	  tt_shrink( tt, size_l );
	  

	  return tt;
}

bool belong_input(mign_node n, std::vector<mign_node> done_inputs) 
{
	for ( auto & x : done_inputs)
	{
		if (x == n)
			return true; 
	}
	return false; 
}

    
void mign_print_stats (const mign_graph& mign, std::ostream& os = std::cout)
{
	const auto name = mign.name();
	auto depth = evaluate_depth (mign); 
	//auto depth_th = evaluate_depth_th (mign);
	auto mign_nc = mign; 
	//auto energy = evaluate_energy (mign_nc);
	//auto energy_leak = leakage_energy(mign_nc,depth_th); 
	auto ce = compute_ce (mign); 
	//auto inv = compute_inv(mign); 
	//auto depth_inv = evaluate_depth_inv(mign); 
	//auto dangling = compute_dangling(mign); 
	os << boost::str( boost::format( "%s i/o = %d/%d , size = %d, depth = %d, ce = %d" ) %( name.empty() ? "(unnamed)" : name ) % mign.inputs().size() % mign.outputs().size() % mign.num_gates() % depth %ce ) << std::endl; 
	os << "****************************************" << std::endl;
	//os << boost::str( boost::format( " Depth threshold = %.2f" ) % depth_th) << std::endl;
	//os << boost::str( boost::format( " Energy threshold = %.2f" ) % energy) << std::endl;
	//os << boost::str( boost::format( " Energy threshold = %.2f" ) % energy_leak) << std::endl;
	//os << "****************************************" << std::endl;
	
	print_depth_statistics(mign); 
	
	//auto larger = compute_larger (mign); 
	//auto smaller = compute_smaller (mign); 
	//os << boost::str( boost::format( "Mig-n: smaller n = %d, larger n = %d") %smaller %larger ) << std::endl; 
	
	print_big (mign, os);
	
}

unsigned input_threshold (unsigned size_tt, unsigned T, unsigned polarity, std::vector<unsigned> weigths)
{

	auto number = size_tt; 

		for (auto& w :weigths)
		{
			if ( w != 1)
			{
				for ( auto x = 2; x<=w; ++x)
				{

					++number; 
					
				}
			}
		
		}
		const auto num = number; 
	
		auto th = 0u; 
		if (polarity == 1)
		{
			th = num - T; 
		}
		else 
			th = T; 
		
	

	if (th > num/2)
	{
		auto narity = 2* th -1; 
		auto filling = narity - num; 
		if (filling > 0)
		{
			for (auto i = 0; i < filling ; ++ i)
			{
				number++; 
			}
		}
		return number; 
		
	}
	else {
		auto narity = 2*(num-th) + 1; 
		auto filling = narity - num; 
		
		if (filling > 0)
		{
			for (auto i = 0; i < filling ; ++ i)
			{
				number++; 
			}
		}
		
		return number; 
}

}
// When you have a vector of mign, with same ENergy and Delay for IMEC, simply the ones that are useless ( same detpth, energy etc)
std::vector<mign_graph> simplify_vector_mign(std::vector<mign_graph>& mign)
{
	if (mign.size() == 1)
		return mign; 
	
	std::vector<std::pair<float, unsigned>> couple, couple_p; 
	
	/*
	for ( auto m : mign)
	{ 
		auto count = 0; 
		auto flag = 1; 
		couple = couple_p; 
		auto depth = evaluate_depth_th(m); 
		auto energy = evaluate_energy(m); 
		for ( auto k : couple)
		{
			if ((k.first <= depth) && (k.second <= energy))
			{
				//count++; 
				flag = 0; 
				break; 
			}
			else if ((k.first == depth) && (k.second > energy))
			{
				//flag = 1; 
				couple_p.erase(couple_p.begin() + count); 
				count--; 
				//couple_p.push_back(std::make_pair(evaluate_depth_th(m),evaluate_energy(m)));
			}
			else if ((k.first > depth) && (k.second == energy))
			{
				//flag = 1; 
				
				couple_p.erase(couple_p.begin() + count); 
				count--; 
				//couple_p.push_back(std::make_pair(evaluate_depth_th(m),evaluate_energy(m)));
			}
			count ++; 
			 
		}
		if ( flag == 1)
			couple_p.push_back(std::make_pair(depth,energy));
	}*/
	
	std::cout << " There are " << mign.size() << " at the beginning. Simplification with Pareto :" << std::endl; 
	
	std::tuple<float, float, unsigned> best_level_pair; // last unsigned is the position in info 
	std::tuple<float, float, unsigned> best_energy_pair; 
	
    auto best_level = std::numeric_limits<float>::max();
    auto best_energy = std::numeric_limits<float>::max();
	
	auto count = 0; 
	
	couple_p = couple; 
	auto mign_p = mign; 
	
	for ( auto x : couple)
	{
		if ( x.first < best_level)
		{
			best_level_pair = std::make_tuple(x.first, x.second, count); 
			best_level = x.first; 
		}
		else if (x.first == best_level)
		{
			
			if (x.second < get<0>(best_level_pair))
			{
				best_level_pair = std::make_tuple(x.first, x.second, count); 
			}
			
		}
		if (x.second < best_energy)
		{
				best_energy_pair = std::make_tuple(x.first, x.second, count); 
				best_energy = x.second; 
		}
		
		else if (x.second == best_energy)
		{
			
			if (x.second < get<1>(best_energy_pair))
			{
				best_energy_pair = std::make_tuple(x.first, x.second, count); 
			}
			
		}
		count++; 
	}
	
	//std::cout << " best level couple " << get<0>(best_level_pair) << " " << get<1>(best_level_pair)<< std::endl; 
	//std::cout << " best energy couple " << get<0>(best_energy_pair) << " " << get<1>(best_energy_pair)<< std::endl; 
	
	if (get<0>(best_level_pair) == get<0>(best_energy_pair) && get<1>(best_level_pair) <= get<1>(best_energy_pair))
	{
		mign.resize(1); 
		mign[0] = mign_p[get<2>(best_level_pair)];
		return mign; 
	}
	
	else if ( get<0>(best_level_pair) == get<0>(best_energy_pair) && get<1>(best_level_pair) > get<1>(best_energy_pair))
	{
		mign.resize(1); 
		mign[0] = mign_p[get<2>(best_energy_pair)];
		return mign; 
	}
	
	else if ( get<0>(best_level_pair) <= get<0>(best_energy_pair) && get<1>(best_level_pair) == get<1>(best_energy_pair))
	{
		mign.resize(1); 
		mign[0] = mign_p[get<2>(best_level_pair)];
		return mign; 
	}
	
	else if ( get<0>(best_level_pair) > get<0>(best_energy_pair) && get<1>(best_level_pair) == get<1>(best_energy_pair))
	{
		mign.resize(1); 
		mign[0] = mign_p[get<2>(best_energy_pair)];
		return mign; 
	}
	
	
	auto flag_time = 0, flag_en = 0 ; 
	count = 0; 
	
	//std::cout << " info size = " << info.size() << std::endl; 
	for ( auto x : couple )
	{
		//std::cout << " couple x = " << x.node_to_level << "   " << x.node_to_energy<< std::endl; 
		
		if ((x.first > get<0>(best_level_pair)) && (x.second >= get<1>(best_level_pair)))
		{
			mign_p.erase(mign_p.begin() + count); 
			//std::cout << " caso 1" << std::endl; 
			count--;
		}
		
		else if ((x.first == get<0>(best_level_pair)) && (x.second > get<1>(best_level_pair)))
		{
			mign_p.erase(mign_p.begin() + count); 
			//std::cout << " caso 2" << std::endl;
			count--;
		}
		else if ((x.first == get<0>(best_level_pair)) && (x.second == get<1>(best_level_pair)))
		{
			if (flag_time == 1)
			{
				mign_p.erase(mign_p.begin() + count); 
				//std::cout << " caso 3" << std::endl;
				count--;
			} 
		    else flag_time = 1; 
		}
		else if ((x.second > get<1>(best_energy_pair)) && (x.first >= get<0>(best_energy_pair)))
		{
			mign_p.erase(mign_p.begin() + count); 
			count--;
			//std::cout << " caso 4" << std::endl;
		}
		else if ((x.second == get<1>(best_energy_pair)) && (x.first > get<0>(best_energy_pair)))
		{
			mign_p.erase(mign_p.begin() + count); 
			count--;
			//std::cout << " caso 5" << std::endl;
		}
		else if ((x.second == get<1>(best_energy_pair)) && (x.first == get<0>(best_energy_pair)))
		{
			if (flag_en == 1)
			{
				mign_p.erase(mign_p.begin() + count); 
				count--;
				//std::cout << " caso 6" << std::endl;
			}
		    else flag_en = 1; 
		}
		//std::cout << " primadi count ==ß" << count << std::endl; 
		count++; 
		//xstd::cout << " primadi count ==ß" << count << std::endl; 
	}
	
	mign.resize(mign_p.size());
	mign = mign_p; 
	
	//std::cout << " info size = " << info.size() << std::endl; 
	
	//std::cout << " best level couple " << get<0>(best_level_pair) << "   " << get<1>(best_level_pair)<< std::endl; 
	//std::cout << " best energy couple " << get<0>(best_energy_pair) << "   " << get<1>(best_energy_pair)<< std::endl; 
	
	std::cout << " There are " << mign.size() << " AFTER" << std::endl; 
	
	return mign; 
}

}
// Local Variables:
// c-basic-offset: 2
// eval: (c-set-offset 'substatement-open 0)
// eval: (c-set-offset 'innamespace 0)
// End:
