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
#include <ctime>

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

class mign_flow_map_manager
{
public:
  mign_flow_map_manager( mign_graph& mign, const properties::ptr& settings );

  void run();

private:
  void find_best_cuts();
  void extract_cover();

private:
  mign_graph&            mign;
  std::vector<unsigned> node_to_cut;
  std::vector<float> node_to_level;
  std::vector<unsigned> node_to_energy;
  std::vector<unsigned> threshold; 
  std::vector<std::vector<unsigned>> weights; 
  std::vector<std::vector<bool>> neg_un; 

  std::shared_ptr<mign_cuts_paged> cuts;
  
  const properties::ptr& settings;

  /* settings */
  bool     verbose;
  unsigned cut_size;
  unsigned priority_cut;
};

mign_flow_map_manager::mign_flow_map_manager( mign_graph& mign, const properties::ptr& settings )
  : mign( mign ), 
   node_to_cut( mign.size() ),
   threshold (mign.size()), 
   weights (mign.size()), 
   neg_un(mign.size()),
   node_to_energy( mign.size() ), 
   node_to_level( mign.size() ), 
   settings(settings)
{
  verbose  = get( settings, "verbose",  true );
  cut_size = get( settings, "cut_size", 6u);
  priority_cut = get( settings, "priority_cut", 0);
}

void mign_flow_map_manager::run()
{
  /* compute cuts */
	clock_t t1,t2,t3,t4;
	
  cuts = std::make_shared<mign_cuts_paged>( mign, cut_size, settings);
  LN( boost::format( "[i] enumerated %d cuts in %.2f secs" ) % cuts->total_cut_count() % cuts->enumeration_time() );

     t1 = clock();
     find_best_cuts();
	 t2 = clock();
	 float diff ((float)t2 - (float)t1);
	 float seconds = diff / CLOCKS_PER_SEC;
	 std::cout << format( "[i] run-time find best cut: %.2f seconds" ) % seconds << std::endl;
	 t3 = clock(); 
	extract_cover();
    t4 = clock();
    float difftwo ((float)t4 - (float)t3);
    seconds = difftwo / CLOCKS_PER_SEC;
    std::cout << format( "[i] run-time extract cover: %.2f seconds" ) % seconds << std::endl;
}

void mign_flow_map_manager::find_best_cuts()
{
	std::cout << " [i] Find best cut" << std::endl;
    for ( auto node : mign.topological_nodes() )
   {
    if ( mign.is_input( node ) )
    {
      assert( cuts->count( node ) == 1u );
      node_to_cut[node] = cuts->cuts( node ).front().address();
      node_to_level[node] = 0u;
	  node_to_energy[node] = 0u;
	  std::vector<unsigned> w ( 1u, 0u); 
	  threshold[node] = 0; 
	  weights[node] = w; 
    }
    else
    {
      auto best_level = std::numeric_limits<float>::max();
	  auto best_energy = std::numeric_limits<unsigned>::max();
      auto best_cut = 0u;
	  auto best_input = 0u; 
	  unsigned best_threshold = 0u; 
	  std::vector<unsigned> best_w; 
	  std::vector<bool> best_un; 
      //std::cout << "************* Node: " << node << std::endl; 
      for ( const auto& cut : cuts->cuts( node ) )
      {
		 
        if ( cut.size() == 1u )  { continue; } /* ignore singleton cuts */   // to comprare with approach in abc max K 9 
			//std::cout << " simulate" << std::endl;
        auto func = cuts->simulate(node,cut); 
		const auto num_vars = tt_num_vars( func );
		//std::cout << " tt " << func << std::endl; 
		//std::cout << "[i] - {" << any_join( cut.range(), ", " ) << "}" << std::endl; 
		
		//std::cout << " compute T w" << std::endl;
		auto result = compute_T_w ( func); 
		
		if (result.t_and_w.first == 0) 
			{ 
				continue; } 
		
		auto T = result.t_and_w.first; 
		//std::cout << " t = " << T << std::endl; 
		auto weig = result.t_and_w.second; 
		for ( auto& e : weig)
		{
			//std::cout << e; 
		}
		//std::cout << std::endl; 
		auto n_of_input = input_threshold (num_vars,T, 0, weig) ; 
		
		//std::cout << " it is a threshold " << std::endl; 
	
		if (n_of_input >= 11)	{
				continue; }
				
			
        float local_max_level = 0u;

        for ( auto leaf : cut )
        {
          local_max_level = std::max( local_max_level, node_to_level[leaf] );
        }

        if ( local_max_level < best_level )
        {
          best_level = local_max_level;
          best_cut = cut.address();
		  best_input = n_of_input; 
		  best_threshold = result.t_and_w.first; 
		  best_w = result.t_and_w.second; 
		  best_un = result.nega_un; 
        }
		
      }

      node_to_cut[node] = best_cut;
	  threshold[node] = best_threshold; 
	  weights[node] = best_w; 
	  neg_un[node] = best_un; 
	 
	  node_to_level[node] = best_level + 1u;
    }
  }
}

void mign_flow_map_manager::extract_cover()
{
  std::cout << " [i]Extract Cover.." << std::endl;
  boost::dynamic_bitset<> visited( mign.size() );

  mign_cover cover( cut_size, mign );
  auto deque = mign_output_deque( mign );

  while ( !deque.empty() )
  {
    auto node = deque.front();
    deque.pop_front();

    if ( mign.is_input( node ) || visited[node] ) { continue; }
    visited[node] = true;
	
    auto cut = cuts->from_address( node_to_cut[node] );
    cover.add_cut( node, cut, threshold[node], weights[node], neg_un[node]);

    for ( auto leaf : cut )
    {
      deque.push_back( leaf );
    }
  }

  mign.set_cover( cover );
}

/******************************************************************************
 * Private functions                                                          *
 ******************************************************************************/

/******************************************************************************
 * Public functions                                                           *
 ******************************************************************************/

void mign_flow_map( mign_graph& mign, const properties::ptr& settings, const properties::ptr& statistics )
{
  mign_flow_map_manager mgr( mign, settings );

  properties_timer t( statistics );
  mgr.run();
}


}

// Local Variables:
// c-basic-offset: 2
// eval: (c-set-offset 'substatement-open 0)
// eval: (c-set-offset 'innamespace 0)
// End:
