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

#include "mign_flow_map.hpp"

#include <vector>

#include <boost/dynamic_bitset.hpp>
#include <boost/format.hpp>
#include <boost/range/iterator_range.hpp>

#include <core/utils/timer.hpp>
#include <core/utils/program_options.hpp>
#include <core/utils/range_utils.hpp>
#include <classical/mign/mign_cover.hpp>
#include <classical/mign/mign_cuts_paged.hpp>
#include <classical/mign/mign_utils.hpp>
#include <classical/mign/threshold_synthesis.hpp>


#define L(x) if ( verbose ) { std::cout << x; }
#define LN(x) if ( verbose ) { std::cout << x << std::endl; }

using namespace boost::program_options;
using boost::format;

namespace cirkit
{

/******************************************************************************
 * Types                                                                      *
 ******************************************************************************/
struct save_info
{
public: 
	
void init (unsigned cut, float level, float energy, unsigned T, std::vector<unsigned> W, std::vector<bool> Neg, unsigned inputs, unsigned i, std::vector<std::pair<unsigned,unsigned>> common)
{
    node_to_cut = cut;
    node_to_level = level;
    node_to_energy = energy;
    threshold = T; 
    weights = W; 
    neg_un = Neg; 
	n_input = inputs; 
	combination_number = i; 	
	common_energy.resize(common.size() ); 
	common_energy = common;
}

float clear(std::vector<std::vector<save_info>>& info, mign_graph& mign)
{
	float energy_common = 0; 
	auto provv = common_energy; 
	
	for ( auto& x : common_energy)
		{
			auto flag = 0, count = 0; 
			for ( auto& y : provv)
			{
				if (mign.is_input(x.first))
					continue; 
				else 
				{
				if ((x.first == y.first) && (x.second == y.second) && (flag == 0))
					flag = 1; 
				else if ((x.first == y.first) && (x.second == y.second) && (flag == 1))
				{
					//std::cout << " COMMON! "<< std::endl; 
					if (info[x.first][x.second].n_input == 3) //|| (info[x.first][x.second].n_input == 5))
						energy_common = energy_common + 3.8; 
					else if ( info[x.first][x.second].n_input == 5)
						energy_common = energy_common + 2.4*2; 
					else if ( info[x.first][x.second].n_input == 7)
						energy_common = energy_common + 3*2; 
					else if (info[x.first][x.second].n_input == 9)
						energy_common = energy_common + 3.6*2; 
					
					provv.erase(provv.begin() + count);
						count--; 
				}
			    }
				count++; 
					
			}
		} 
		
		common_energy = provv; 
		return energy_common; 
}

    unsigned node_to_cut;
    float node_to_level;
    float node_to_energy;
    unsigned threshold; 
    std::vector<unsigned> weights; 
    std::vector<bool> neg_un; 
	unsigned n_input; 
	unsigned combination_number; 
	std::vector<std::pair<unsigned, unsigned>> common_energy;   // wichi leaf and in which combinations. Save for each node all the paths until PI. 
}; 


std::vector<save_info> simplify_info (std::vector<save_info>& info) // function to simply according to Pareto dominance ;) 
{
	if (info.size() == 1)
		return info; 
	
	std::tuple<float, float, unsigned> best_level_pair; // last unsigned is the position in info 
	std::tuple<float, float, unsigned> best_energy_pair; 
	
    auto best_level = std::numeric_limits<float>::max();
    auto best_energy = std::numeric_limits<float>::max();
	
	auto count = 0; 
	
	auto info_p = info; 
	
	for ( auto x : info)
	{
		if ( x.node_to_level < best_level)
		{
			best_level_pair = std::make_tuple(x.node_to_level, x.node_to_energy, count); 
			best_level = x.node_to_level; 
		}
		else if (x.node_to_level == best_level)
		{
			
			if (x.node_to_energy < get<0>(best_level_pair))
			{
				best_level_pair = std::make_tuple(x.node_to_level, x.node_to_energy, count); 
			}
			
		}
		if (x.node_to_energy < best_energy)
		{
				best_energy_pair = std::make_tuple(x.node_to_level, x.node_to_energy, count); 
				best_energy = x.node_to_energy; 
		}
		
		else if (x.node_to_energy == best_energy)
		{
			
			if (x.node_to_energy < get<1>(best_energy_pair))
			{
				best_energy_pair = std::make_tuple(x.node_to_level, x.node_to_energy, count); 
			}
			
		}
		count++; 
	}
	
	//std::cout << " best level couple " << get<0>(best_level_pair) << " " << get<1>(best_level_pair)<< std::endl; 
	//std::cout << " best energy couple " << get<0>(best_energy_pair) << " " << get<1>(best_energy_pair)<< std::endl; 
	
	
	// PER FARNE AVERE SOLAMENTE DUE : 
	

	if (get<0>(best_level_pair) == get<0>(best_energy_pair) && get<1>(best_level_pair) <= get<1>(best_energy_pair))
	{
		info.resize(1); 
		info[0] = info_p[get<2>(best_level_pair)];
		return info; 
	}
	
	else if ( get<0>(best_level_pair) == get<0>(best_energy_pair) && get<1>(best_level_pair) > get<1>(best_energy_pair))
	{
		info.resize(1); 
		info[0] = info_p[get<2>(best_energy_pair)];
		return info; 
	}
	
	else if ( get<0>(best_level_pair) <= get<0>(best_energy_pair) && get<1>(best_level_pair) == get<1>(best_energy_pair))
	{
		info.resize(1); 
		info[0] = info_p[get<2>(best_level_pair)];
		return info; 
	}
	
	else if ( get<0>(best_level_pair) > get<0>(best_energy_pair) && get<1>(best_level_pair) == get<1>(best_energy_pair))
	{
		info.resize(1); 
		info[0] = info_p[get<2>(best_energy_pair)];
		return info; 
	}
	
	
    info[0] = info_p[get<2>(best_level_pair)]; 
    info[1] = info_p[get<2>(best_energy_pair)]; 
	
	
	auto flag_time = 0, flag_en = 0 ; 
	count = 0; 
	
	//std::cout << " info size = " << info.size() << std::endl; 
	for ( auto x : info )
	{
		//std::cout << " couple x = " << x.node_to_level << "   " << x.node_to_energy<< std::endl; 
		
		if ((x.node_to_level > get<0>(best_level_pair)) && (x.node_to_energy >= get<1>(best_level_pair)))
		{
			info_p.erase(info_p.begin() + count); 
			//std::cout << " caso 1" << std::endl; 
			count--;
		}
		
		else if ((x.node_to_level == get<0>(best_level_pair)) && (x.node_to_energy > get<1>(best_level_pair)))
		{
			info_p.erase(info_p.begin() + count); 
			//std::cout << " caso 2" << std::endl;
			count--;
		}
		else if ((x.node_to_level == get<0>(best_level_pair)) && (x.node_to_energy == get<1>(best_level_pair)))
		{
			//if (flag_time == 1)
			//{
				info_p.erase(info_p.begin() + count); 
				//std::cout << " caso 3" << std::endl;
				count--;
				//} 
		   // else flag_time = 1; 
		}
		else if ((x.node_to_energy > get<1>(best_energy_pair)) && (x.node_to_level >= get<0>(best_energy_pair)))
		{
			info_p.erase(info_p.begin() + count); 
			count--;
			//std::cout << " caso 4" << std::endl;
		}
		else if ((x.node_to_energy == get<1>(best_energy_pair)) && (x.node_to_level > get<0>(best_energy_pair)))
		{
			info_p.erase(info_p.begin() + count); 
			count--;
			//std::cout << " caso 5" << std::endl;
		}
		else if ((x.node_to_energy == get<1>(best_energy_pair)) && (x.node_to_level == get<0>(best_energy_pair)))
		{
			//if (flag_en == 1)
			//{
				info_p.erase(info_p.begin() + count); 
				count--;
				//std::cout << " caso 6" << std::endl;
				//}
		    //else flag_en = 1; 
		}
		//std::cout << " primadi count ==ß" << count << std::endl; 
		count++; 
		//xstd::cout << " primadi count ==ß" << count << std::endl; 
	}
	
	auto h = 2; 
	for ( auto f = 0; f <info_p.size(); ++f)
	{
		info[h] = info_p[f]; 
		h++; 
	}
	
	if ( info.size() > 4)
		info.resize(4); 
	
	// SELECT ONLY 4 of the POSSIBILITIES TO CUT THE RUNNING TIME: 
	
		
	return info; 
} 

class mign_flow_map_multiobj_manager
{
public:
  mign_flow_map_multiobj_manager( mign_graph& mign, const properties::ptr& settings );

  void run();

private:
  void find_best_cuts();
  void extract_cover();

private:
  mign_graph&            mign;
  std::vector<std::vector<save_info>>             info; 

  std::shared_ptr<mign_cuts_paged> cuts;
  
  const properties::ptr& settings;

  /* settings */
  bool     verbose;
  unsigned cut_size;
};

mign_flow_map_multiobj_manager::mign_flow_map_multiobj_manager( mign_graph& mign, const properties::ptr& settings )
  : mign( mign ), 
   settings(settings),
   info( mign.size() )
   
{
  verbose  = get( settings, "verbose",  true );
  cut_size = get( settings, "cut_size", 6u);
}

void mign_flow_map_multiobj_manager::run()
{
  /* compute cuts */
  cuts = std::make_shared<mign_cuts_paged>( mign, cut_size, settings );
  LN( boost::format( "[i] enumerated %d cuts in %.2f secs" ) % cuts->total_cut_count() % cuts->enumeration_time() );

    find_best_cuts();
	extract_cover();
}

void mign_flow_map_multiobj_manager::find_best_cuts()
{
	std::cout << " [i] Find best cut" << std::endl;
	
	//std::vector<std::vector<save_info>> info; 
		
  for ( auto node : mign.topological_nodes() )
  {
    if ( mign.is_input( node ) )
    {
      assert( cuts->count( node ) == 1u );
	  std::vector<unsigned> w ( 1u, 0u); 
	  std::vector<bool> neg_un ( 1u, true); 
	  info[node].resize(1); 
	  std::vector<std::pair<unsigned, unsigned>> common; 
	  common.push_back(std::make_pair(node,0)); 
	  //save_info provv; 
	  //provv.init(cuts->cuts( node ).front().address(), 0u, 0u, 0u, w, neg_un, 1u,  0u); 
	  info[node][0].init(cuts->cuts( node ).front().address(), 0, 0, 0u, w, neg_un, 1u,  0u, common); 
	 
    }
    else
    {
		//std::cout << " NODO = " << node << std::endl; 
		
		std::vector<save_info> info_cuts; 
		
		auto number_of_com = 0; 
		  auto iter = 0u; 
		
      for ( const auto& cut : cuts->cuts( node ) )
      {
		 
        if ( cut.size() == 1u )  { continue; } /* ignore singleton cuts */   // to comprare with approach in abc max K 9 
			
        auto func = cuts->simulate(node,cut); 
		const auto num_vars = tt_num_vars( func );
		
		//std::cout << "[i] - {" << any_join( cut.range(), ", " ) << "}" << std::endl; 
		
		auto result = compute_T_w (func); 
		
		//std::cout << " DOpo threshold" << std::endl; 
		
		if (result.t_and_w.first == 0) 
			{ 
				continue; } 
		
		auto T = result.t_and_w.first; 
		auto weig = result.t_and_w.second; 
		//std::cout << "N of inputs" << std::endl; 
		auto n_of_input = input_threshold (num_vars,T, 0, weig) ; 
		auto neg_un = result.nega_un; 
		
		if (n_of_input >= 11)	{
				continue; }
				
        // for the threshold gates implementation by Arizona State University I need also the energy 

  	    const auto num_children = cut.size(); // numero di foglie --> figli ;) ogni leaf ha piu di una combinazione
		
  	    std::vector<unsigned> a(num_children+1,0u); 
  	    std::vector<unsigned> m(num_children+1,2u); 
		
  	  auto y = 1; 
	  for ( auto leaf : cut )
	  {
		 
	  	 m[y] = info[leaf].size(); 
		
		 ++y; 
		
	  }
	  	std::vector<std::vector<unsigned>> leaves; 
	  
	  //std::cout << "Mixed radix" << std::endl; 
	    mixed_radix( a, m, [&]( const std::vector<unsigned>& a ) {
        
    	  for (auto h = 0; h < a.size(); ++h)
    	  {
			  // std::cout << a[h]; 
	  
    	  }
		  // std::cout << std::endl; 
		   
			auto c = a; 
			c.erase (c.begin() + 0); 
	 
	    	  for (auto h = 0; h < c.size(); ++h)
	    	  {
				   //std::cout << c[h]; 
		  
	    	  }
			 //  std::cout << std::endl; 
	  	 leaves.push_back(c); 
	  
	          return false;
	        } );
	  
  	   
		//std::vector<std::vector<unsigned>> leaves; 
		
  	 /* while ( true )
  	  {
		  number_of_com++; 
		
		  leaves.push_back(a); 
		   
		  
  		  auto flag = 0; 
  		  for ( auto i = 1u; i <= num_children; ++i )
  		  {
  			 
  		    if (a[i] == m[i] -1) //|| (count(ns[c]) == 1))
  		    {
  		    	++flag; 
  				
  		    }
		
  			//std::cout << i << std::endl; 
  		  }
  		  if (flag == num_children)
  		  break; 
  	    // do something with as:
  	     auto j = a.size()-1;    
  	    while ( a[j] == m[j]-1 )
  	    {
			
  		 // std::cout << " j prima " << j << std::endl; 
  			auto x = j--;
  			
  		   a[x] = prova[0]; 
  		 
  	    }
  	 if (!j) {break;}
  	  a[j]++;  
  
  	   }*/
	   

        // prendo tutte le combinazioni delle leaf.  TODO  e per ogni combinazione i che va fino al numro delle combianzioni calcolo il local max. 
	 info_cuts.resize(leaves.size()); 
	 
	   /// to save common leaves for energy , becase some nodes could be common to different subgraphs 
	 iter = 0; 
        for ( auto x : leaves )
        {
	        float local_max_level = 0;
			float local_max_energy = 0;
			
			auto c = 0; // because a is longer 
			std::vector<std::pair<unsigned,unsigned>> common;
			
			for ( auto leaf : cut)
				
			{
				//std::cout << " x[c] = " << x[c] << std::endl; 
				//std::cout << " leaf " << leaf << std::endl; 
    			  local_max_level = std::max( local_max_level, info[leaf][x[c]].node_to_level );  // <-- ela combinazione i 
    		      local_max_energy = local_max_energy + info[leaf][x[c]].node_to_energy;
				  
				  for (auto h = 0; h <info[leaf][x[c]].common_energy.size(); h++)
				  {
					  common.push_back(info[leaf][x[c]].common_energy[h]); 
				  }
				  c++; 
			}
  
			  // some nodes could be shared between subgraphs ;) 
			
			info_cuts[iter].init(cut.address(),local_max_level, local_max_energy, T, weig, neg_un, n_of_input, iter, common); 
			auto energy_common = info_cuts[iter].clear(info,mign); 
			info_cuts[iter].node_to_energy = info_cuts[iter].node_to_energy - energy_common;
			//std::cout << " energy = " << info_cuts[iter].node_to_energy << " delay =" << info_cuts[iter].node_to_level << std::endl; 
			
  		  if (info_cuts[iter].n_input == 3)
  		  {
  			  info_cuts[iter].node_to_level = info_cuts[iter].node_to_level + 40;
  			 // std::cout << "node to level" << node_to_level[node] << std::endl;  
  			  info_cuts[iter].node_to_energy = info_cuts[iter].node_to_energy + 3.8 + 0.49;
  			  //std::cout << "node to energy" << node_to_energy[node] << std::endl;  
  			 // info_cuts[iter].common_energy.push_back(std::make_pair(node,iter)); 
			  
  		  }
  		  else if (info_cuts[iter].n_input == 5)
  		  {
  			  info_cuts[iter].node_to_level = info_cuts[iter].node_to_level + 44;
  			 // std::cout << "node to level" << node_to_level[node] << std::endl;  
  			  info_cuts[iter].node_to_energy = info_cuts[iter].node_to_energy + 4.8 + 0.49;
  			  //std::cout << "node to energy" << node_to_energy[node] << std::endl;  
  			 // info_cuts[iter].common_energy.push_back(std::make_pair(node,iter)); 
  		  }
  		  else if (info_cuts[iter].n_input == 7)
  		  {
  			  info_cuts[iter].node_to_level = info_cuts[iter].node_to_level + 50;
  			 // std::cout << "node to level" << node_to_level[node] << std::endl;  
  			  info_cuts[iter].node_to_energy = info_cuts[iter].node_to_energy + 6 + 0.49;
  			  //std::cout << "node to energy" << node_to_energy[node] << std::endl;  
  			  //info_cuts[iter].common_energy.push_back(std::make_pair(node,iter)); 
  		  }
  		  else if (info_cuts[iter].n_input == 9)
  		  {
  			  info_cuts[iter].node_to_level = info_cuts[iter].node_to_level + 60;
  			 // std::cout << "node to level" << node_to_level[node] << std::endl;  
  			  info_cuts[iter].node_to_energy = info_cuts[iter].node_to_energy + 7.2 + 0.49;
  			  //std::cout << "node to energy" << node_to_energy[node] << std::endl; 
			  //info_cuts[iter].common_energy.push_back(std::make_pair(node,iter));  
  		  }
		  
		   
			++iter; 
		}
			
      }
	  
	  //std::cout << " simplify info" << std::endl; 
	  simplify_info (info_cuts); 
	 
	  for ( auto x = 0; x < info_cuts.size(); ++x)
	  {
		  info_cuts[x].common_energy.push_back(std::make_pair(node,x)); 
	  }
	  
	  info[node].resize(info_cuts.size()); 
	  info[node] = info_cuts; 
    }
  }
}

unsigned find_in_comb(std::vector<mign_node> comb, mign_node node)
{
	auto count = 0; 
	for (auto& j : comb)
	{
		if (j == node)
			return count; 
		else 
			count++; 
	}
	
	return count; 
}

std::pair<std::vector<std::vector<unsigned>>, std::vector<mign_node>> node_cut_comb (std::vector<std::vector<save_info>> info)
{
	
	
	auto num_children = 0; 
	std::vector<mign_node> nodes; 
	
	for ( auto& y : info)
	{
		if (y.size() > 1)
			num_children++; 	
	}
 
    std::vector<unsigned> a(num_children+1,0u); 
    std::vector<unsigned> m(num_children+1,2u); 
    std::vector<unsigned> prova(num_children+1,0u);
	
	  auto y = 1u; 
	  auto i = 0u; 
	  for ( auto& x: info )
	  {
		  
		 if (x.size() > 1)
	  	{
	  	
	  	 m[y] = x.size(); 
		// std::cout << " m[y]" << m[y] << std::endl;  
		 nodes.push_back(i); 
		 ++y; 
	    }
		++i; 
		 //++y; 
	  }
	  
	std::pair<std::vector<std::vector<unsigned>>, std::vector<mign_node>> combinations; 
	combinations.second = nodes; 
	
	auto h = 0; 
  while ( true)
  {
	  auto c = a; 
	  c.erase (c.begin() + 0); 
  
	  combinations.first.push_back(c); 
	  
	  auto flag = 0; 
	  for ( auto i = 1u; i <= num_children; ++i )
	  {
		 
	    if (a[i] == m[i] -1) //|| (count(ns[c]) == 1))
	    {
	    	++flag; 	
	    }

	  }
	  
	  if (flag == num_children)
	  break; 
	  
    // do something with as:
     auto j = a.size()-1;    
    while ( a[j] == m[j]-1 )
    {
		
		auto x = j--;
		
	   a[x] = prova[0]; 
	 
    }
   if (!j) {break;}
  a[j]++;  

  h++; 
   }
   
   
	return combinations; 
}

std::pair<std::vector<std::vector<unsigned>>, std::vector<mign_node>> output_combinations (std::vector<std::vector<save_info>> info, mign_graph& mign)
{
	
	auto num_children = mign.outputs().size(); 
	//std::cout << " num children = "<< num_children << std::endl; 
	
    std::vector<unsigned> a( num_children+1, 0u );
    std::vector<unsigned> m (num_children+1,2u); 
 std::vector<mign_node> nodes; 
 
  auto y = 1u; 
  auto i = 0u; 
  for ( auto& x: info )
  {
	  if (mign.is_output(i))
  	 {
		 m[y] = x.size(); 
		
		 nodes.push_back(i); 
	     ++y; 
		 ++i; 
    }
	else 
	++i; 
  }
  
   
	std::vector<std::vector<unsigned>> combinations; 
	
    mixed_radix( a, m, [&]( const std::vector<unsigned>& a ) {
        
		
	
  	  auto c = a; 
  	  c.erase (c.begin() + 0); 
	 
  	  for (auto h = 0; h < c.size(); ++h)
  	  {
  		 // std::cout << c[h]; 
		  
  	  }
  	 // std::cout << std::endl; 
	 
  	  combinations.push_back(c); 
	  
          return false;
        } );
	
   
	return std::make_pair(combinations, nodes); 
}

unsigned get_comb_position(std::vector<mign_node> comb, mign_node node)
{
	auto count = 0; 
	for (auto& j : comb)
	{
		if (j == node)
			return count; 
		else 
			count++; 
	}
	
	return count; 
}

mign_node which_node (mign_node leaf, std::vector<std::pair<mign_node,std::vector<std::pair<unsigned, unsigned>>>> path_to_PI)
{
	for ( auto& g : path_to_PI)
	{
		for (auto& h : g.second)
		{
			if (h.first == leaf)
				return g.first; 
		}
	}
	return leaf; 
}

void mign_flow_map_multiobj_manager::extract_cover()
{
	// ci sono tant MIGS 
  std::cout << " [i] Extract Cover.." << std::endl;
  
  auto i = 0; 
  
    auto combinations = output_combinations (info, mign);   // ot spped up the algorithm, I will consider only putputs combinations since I already know the path until them
	
    //std::cout << "  avro " << combinations.first.size() << " covers" << std::endl; 
   
	const auto total = combinations.first.size(); 
    std::vector<mign_cover> covers; 
   
	covers.resize(total);  

   
    auto comb = 0; 
    for ( auto& c : combinations.first)
    {
  	 // std::cout << " c size " << c.size()  << std::endl; 
  	 // std::cout << " comb = " << comb << std::endl; 
  	  boost::dynamic_bitset<> visited( mign.size() );
  	  mign_cover covers_provv (cut_size, mign); 
  	  covers[comb] = covers_provv; 
  	  auto deque = mign_output_deque( mign );
	  std::vector<mign_node> leaves; 
	  //std::map<mign_node, std::vector<mign_node>> leaves; 
	  auto count = 0; 
	  auto flag = 0; 
	  
	 std::vector<std::pair<mign_node,std::vector<std::pair<unsigned, unsigned>>>> path_to_PI; 
	  //path_to_PI.resize(mign.outputs().size()); 
		  
  	  while ( !deque.empty() )
  	  {
		  
  	    auto node = deque.front();
		deque.pop_front();
		
		auto it = find (leaves.begin(), leaves.end(), node); 
		if (mign.is_output(node) && it == leaves.end())  // can happend that is a output but  it is also a leaves for some reason 
		{
			//std::cout << " output = " << node << std::endl; 
			
			auto position = get_comb_position (combinations.second,node); 
			//std::cout << " position " << position << std::endl; 
			auto cut = cuts->from_address( info[node][c[position]].node_to_cut );
			covers[comb].add_cut( node, cut, info[node][c[position]].threshold, info[node][c[position]].weights, info[node][c[position]].neg_un);
			visited[node] = true; 
			
					  	    for ( auto leaf : cut )
					  	    {
								leaves.push_back(leaf); 
					  	      deque.push_back( leaf );
					  	    }
							path_to_PI.push_back({node,info[node][c[position]].common_energy}); 
							continue; 
		}
		
  	    if ( mign.is_input( node ) || visited[node] ) { continue; }
  	    visited[node] = true;
	
		auto t = which_node (node, path_to_PI);   // e possibile che un nodo abbia due configurazioni? eh va beh prende la rpima in ogni caso e pace 
		auto position = get_comb_position (combinations.second,t); 
		//std::cout << " common energy node " << node << std::endl; 
			for ( auto& r : info[t][c[position]].common_energy) 
			{
				//std::cout << "r.first " << r.first << " n comb " << r.second << std::endl; 
				if (r.first != node)
					continue; 
				else 
				{
					//std::cout << " size =" << info[node].size() << std::endl; 
					auto cut = cuts->from_address( info[node][r.second].node_to_cut );
				covers[comb].add_cut( node, cut, info[node][r.second].threshold, info[node][r.second].weights, info[node][r.second].neg_un);
			
		  	    for ( auto leaf : cut )
		  	    {
					leaves.push_back(leaf); 
		  	      deque.push_back( leaf );
		  	    }
				break; 
			}
			}
  			
  	  }

  	  comb++; 
    }
	/*
  for ( auto& c : combinations.first)
  {
	 // std::cout << " c size " << c.size()  << std::endl; 
	 // std::cout << " comb = " << comb << std::endl; 
	  boost::dynamic_bitset<> visited( mign.size() );
	  mign_cover covers_provv (cut_size, mign); 
	  covers[comb] = covers_provv; 
	  auto deque = mign_output_deque( mign );
	  
	  while ( !deque.empty() )
	  {
	    auto node = deque.front();
		
	    deque.pop_front();

	    if ( mign.is_input( node ) || visited[node] ) { continue; }
	    visited[node] = true;
	
		if (info[node].size() == 1)
		{
			//std::cout << " nODo con siz = 1 " << node << std::endl; 
			auto cut = cuts->from_address( info[node][0].node_to_cut );
	    covers[comb].add_cut( node, cut, info[node][0].threshold, info[node][0].weights, info[node][0].neg_un, -1, -1);

	    for ( auto leaf : cut )
	    {
			
	      deque.push_back( leaf );
	    }
	    }
		else 
		{
			//std::cout << " nODo con size > 1 " << node << std::endl; 
			auto position = find_in_comb(combinations.second,node); 
			//std::cout << "position = " << position << std::endl; 
			auto cut = cuts->from_address( info[node][c[position]].node_to_cut );
			//std::cout << " c[position]" << c[position] << std::endl; 
	        covers[comb].add_cut( node, cut, info[node][c[position]].threshold, info[node][c[position]].weights, info[node][c[position]].neg_un, -1, -1);

	    for ( auto leaf : cut )
	    {
	      deque.push_back( leaf );
	    }
	}
	  
	  }

	  comb++; 
  }*/
   
 mign.set_multi_cover( covers );
 
}

/******************************************************************************
 * Private functions                                                          *
 ******************************************************************************/

/******************************************************************************
 * Public functions                                                           *
 ******************************************************************************/

void mign_flow_map_opt( mign_graph& mign, const properties::ptr& settings, const properties::ptr& statistics )
{
  mign_flow_map_multiobj_manager mgr( mign, settings );

  properties_timer t( statistics );
  mgr.run();
}


}

// Local Variables:
// c-basic-offset: 2
// eval: (c-set-offset 'substatement-open 0)
// eval: (c-set-offset 'innamespace 0)
// End:
