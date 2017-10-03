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

#include "mign_flow_almost_maj.hpp"

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

#define timer timer_class
#include <boost/progress.hpp>
#undef timer

#define L(x) if ( verbose ) { std::cout << x; }
#define LN(x) if ( verbose ) { std::cout << x << std::endl; }

//using namespace boost::program_options;
//using boost::format;

namespace cirkit
{

/******************************************************************************
 * Types                                                                      *
 ******************************************************************************/

	bool is_almost_maj (tt func, int& almost_node_is_maj, tt& new_tt)
	{
		// check wheter the tt func is a almost a majority. And also check what is the truth table that we have to add to have a majority ;)  
		/*auto maj = tt_num_vars( func ); 
	    tt tt_maj( 1u << maj );
	    boost::dynamic_bitset<> it( maj, 0 );
		int flag = 0;

	    do
	    {
	      m[it.to_ulong()] = it.count() > ( maj >> 1u );
	      inc( it );
	    } while ( it.any() );
	
		tt_allign(func,tt_maj); 
		auto counter = 0u; 
		new_tt.resize(tt_num_vars( func )); 
	
		for (auto i = 0u; i < tt_maj.sie(); i++)
		{
			if (func[i] == tt_maj[i])
			{
				counter++; 
				tt_new[i] = 0; 
			}
			else 
			{
				if (func[i] > tt_maj[i])
				{
					if ((flag == 0) || (flag == 1))
					{
						flag = 1; 
						tt_new[i] = 1;
						almost_node_is_maj = 1; 
					}
					else if (flag == -1)
					{
						return false; 
					}
				}
				else if (func[i] < tt_maj[i])
				{
					if ((flag == 0) || (flag == -1))
					{
						flag = -1; 
						tt_new[i] = 1;
						almost_node_is_maj = 0; 
					}
					else if (flag == 1)
					{
						return false; 
					}
				}
			}
		}
		if (counter > tt_maj - 10)
		{
			return true; 
		}
		else return false;*/ 
	
		return false; 
	}

	int find_tt_in_cut (tt new_tt, std::vector<tt> tt_to_cut, mign_node node)
	{
	
		/*for (auto n = 0; n < node; n++)
		{
			auto cut_tt = tt_to_vut[n]; 
		
		}*/
		return -1; 
	}

	bool is_maj_with_dcs (tt func)  // non é cosi facile su un grafico in movimento perché cambierá tutto. 
	{
		return false; 
	}
	
class mign_flow_almost_maj_manager
{
public:
  mign_flow_almost_maj_manager( mign_graph& mign, const properties::ptr& settings );

  void run();

private:
  void find_best_cuts();
  void extract_cover();

private:
  mign_graph&                         mign;
  std::vector<unsigned>               node_to_cut;
  std::vector<float>                  node_to_level;
  std::vector<unsigned>               threshold; 
  std::vector<std::vector<unsigned>>  weights; 
  std::vector<std::vector<bool>>      neg_un; 
  std::vector<int>                    almost_tt; 
  std::vector<int>                    almost_tt_is_maj; 
  std::shared_ptr<mign_cuts_paged>    cuts;
  std::vector<tt>                     tt_to_cut; 
  
  const properties::ptr& settings;

  /* settings */
  bool     verbose;
  unsigned cut_size;
  unsigned priority_cut;
  unsigned allow_almost;
  bool progress; 
};

mign_flow_almost_maj_manager::mign_flow_almost_maj_manager( mign_graph& mign, const properties::ptr& settings )
  : mign( mign ), 
   node_to_cut( mign.size() ),
   threshold (mign.size()), 
   weights (mign.size()), 
   neg_un(mign.size()),
   almost_tt(mign.size()),
   almost_tt_is_maj(mign.size()),
   node_to_level( mign.size() ), 
   settings(settings)
{
  verbose  = get( settings, "verbose",  true );
  cut_size = get( settings, "cut_size", 6u);
  progress = get( settings, "progress", true);
  priority_cut = get( settings, "priority_cut", 0u);
  allow_almost = get( settings, "allow_almost", 0u);
}

void mign_flow_almost_maj_manager::run()
{
  /* compute cuts */
  auto cuts_settings = std::make_shared<properties>();
  cuts_settings->set( "progress", progress );
  cuts = std::make_shared<mign_cuts_paged>( mign, cut_size, cuts_settings);
  LN( boost::format( "[i] enumerated %d cuts in %.2f secs" ) % cuts->total_cut_count() % cuts->enumeration_time()) ;

  find_best_cuts();
  extract_cover();
  
}

void mign_flow_almost_maj_manager::find_best_cuts()
{
	std::cout << " [i] Find best cut" << std::endl;

    for ( auto node : mign.topological_nodes() )
   {
    if ( mign.is_input( node ) )
    {
      assert( cuts->count( node ) == 1u );
      node_to_cut[node] = cuts->cuts( node ).front().address();
      node_to_level[node] = 0u;
    }
    else
    {
      auto best_level = std::numeric_limits<float>::max();
      auto best_cut = 0u;
	  auto best_input = 0u; 
	  auto best_threshold = 0u; 
	  std::vector<unsigned> best_w; 
	  std::vector<bool> best_un; 
	  auto best_almost_node = -1; 
	  auto best_almost_is_maj = -1; 
	  
	  tt tt_to_func; 
	  
      for ( const auto& cut : cuts->cuts( node ) )
      {
		 
        if ( cut.size() == 1u )  { continue; } /* ignore singleton cuts */   
        
		auto func = cuts->simulate(node,cut); 
		const auto num_vars = tt_num_vars( func );
		tt_to_func.resize(num_vars); 
		tt_to_func = func; 
		
		if (num_vars >= 11) continue; // Maj-5, MAJ-7 and MAJ-9
		
		auto result = compute_T_w (func); 
		auto T = result.t_and_w.first; 
		auto weig = result.t_and_w.second; 
		auto n_of_input = input_threshold (num_vars,T,0, weig);
		auto almost_node = -1; 
		auto almost_node_is_maj = -1; 
		
		if ((result.t_and_w.first == 0) || (n_of_input >= 11))
	    { 
			if (allow_almost == 0)
			{
				/*if (is_maj_with_dcs (func)) // check if it is a majority with dcs
				{
					T = num_vars/2 + 1; 
					for ( auto i = 0u; i < num_vars; i++)
					{
						weig[i] = 1; 
					}
				}*/
			}
			if (allow_almost == 1)
			{
				tt new_tt; 
				if (is_almost_maj (func, almost_node_is_maj, new_tt)) // means it is bigger
				{
						almost_node = find_tt_in_cut (new_tt, tt_to_cut, node);  // funzione che restituisce il nodo da usare :) 
				}
				else 
				{
					almost_node_is_maj = -1;  
					almost_node = -1;
				}
				/*  Algorithm: 
				
				1. Cerco se la truth table é "quasi " un MAJ -- decidere cosa vuol dire QUASI un maj 
				2. Cerco di capire se é maggiore o minore...
				3. Calcolo la trth table che dovrei sommarci - sottrarci 
				4. Cerco se tale tt o tt con input negati giá esise tra quelle che ci sono. 
				
				  */
				
			}
			else 
			{
				almost_node_is_maj = -1;  
				almost_node = -1; 
			}
			 
			if (( almost_node_is_maj >= 0) && (almost_node >= 0)) 
			{
				
		        float local_max_level = 0u;
				for ( auto leaf : cut )
		        {
		          local_max_level = std::max( local_max_level, node_to_level[leaf] + 1 );
		        }
				local_max_level = std::max( local_max_level, node_to_level[almost_node] + 1 );
		        
				if ( local_max_level < best_level )
		        {
		          best_level = local_max_level;
		          best_cut = cut.address();
				  best_input = num_vars; 
				  best_threshold = 0; 
				  std::vector<unsigned> best_w; 
				  std::vector<bool> best_un; 
				  best_almost_node = almost_node;  
				  best_almost_is_maj = almost_node_is_maj;
		        }
			}
			else 
				continue; 
	    } 
		
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
		  best_almost_node = -1; 
		  best_almost_is_maj = -1; 
        }
		
      }
	  
      tt_to_cut.push_back(tt_to_func); 
      node_to_cut[node] = best_cut;
	  threshold[node] = best_threshold; 
	  weights[node] = best_w; 
	  neg_un[node] = best_un; 
	  almost_tt[node] = best_almost_node; 
	  almost_tt_is_maj[node] = best_almost_is_maj;
	  
	  node_to_level[node] = best_level + 1u;
    }
  }
}

void mign_flow_almost_maj_manager::extract_cover()
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
    cover.add_cut( node, cut, threshold[node], weights[node], neg_un[node], almost_tt[node], almost_tt_is_maj[node]);

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

void mign_flow_almost_maj ( mign_graph& mign, const properties::ptr& settings, const properties::ptr& statistics )
{
  mign_flow_almost_maj_manager mgr( mign, settings );

  properties_timer t( statistics );
  mgr.run();
}


}

// Local Variables:
// c-basic-offset: 2
// eval: (c-set-offset 'substatement-open 0)
// eval: (c-set-offset 'innamespace 0)
// End:
