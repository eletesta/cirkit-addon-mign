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
#include <classical/mign/mign_utils_majtt.hpp>
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
	
class mign_flow_almost_maj_manager
{
public:
  mign_flow_almost_maj_manager( mign_graph& mign, const properties::ptr& settings );

  void run();

private:
  void find_best_cuts();
  void extract_cover();

private:
  mign_graph&                          mign;
  std::vector<unsigned>                node_to_cut;
  std::vector<float>                   node_to_level;
  std::shared_ptr<mign_cuts_paged>     cuts;
  std::vector<boost::dynamic_bitset<>> tt_is_maj; 
  std::vector<unsigned>                tt_is_almostmaj_or_maj; 
  std::vector<tt>                      reminder; 
  std::vector<std::string>             exact;
  
  const properties::ptr& settings;

  /* settings */
  bool     verbose;
  unsigned cut_size;
  unsigned extra; 
  unsigned priority_cut;
  unsigned allow_almost;
  bool progress; 
};

mign_flow_almost_maj_manager::mign_flow_almost_maj_manager( mign_graph& mign, const properties::ptr& settings )
  : mign( mign ), 
   node_to_cut( mign.size() ),
   node_to_level( mign.size() ), 
   tt_is_maj(mign.size()),
   tt_is_almostmaj_or_maj (mign.size()),
   reminder (mign.size()), 
   exact (mign.size()),
   settings(settings)
{
  verbose      = get( settings, "verbose",  true );
  cut_size     = get( settings, "cut_size", 8u);
  progress     = get( settings, "progress", true);
  extra        = get( settings, "extra", 0u);
  priority_cut = get( settings, "priority_cut", 0u);
  allow_almost = get( settings, "allow_almost", 0u);
}

void mign_flow_almost_maj_manager::run()
{
  /* compute cuts */
  auto cuts_settings = std::make_shared<properties>();
  cuts_settings->set( "progress", progress );
  cuts_settings->set( "extra", extra );
  cuts_settings->set( "almost", allow_almost );
  cuts = std::make_shared<mign_cuts_paged>( mign, cut_size, cuts_settings);
  LN( boost::format( "[i] enumerated %d cuts in %.2f secs" ) % cuts->total_cut_count() % cuts->enumeration_time()) ;

  find_best_cuts();
  extract_cover();
  
}

void mign_flow_almost_maj_manager::find_best_cuts()
{
	std::cout << " [i] Find best cut ... " << std::endl;

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
	  boost::dynamic_bitset<> best_tt_is_maj; 
	  unsigned best_tt_is_almostmaj_or_maj; 
	  tt best_tt_reminder;
	  std::string best_exact;  
	  
	  tt tt_to_func; 

	  auto cut_c = -1; 
	  for ( const auto& cut : cuts->cuts( node ) )
      {
		cut_c++; 
        auto flag = 0u; 
        if ( cut.size() == 1u )  { continue; } /* ignore singleton cuts */   
		/*auto func = cuts->simulate(node,cut); 
		const auto num_vars = tt_num_vars( func );
		tt_to_func.resize(num_vars); 
		tt_to_func = func; */
		
		//if (num_vars >= 9) continue; 
	    boost::dynamic_bitset<> local_tt_is_maj; 
		unsigned local_tt_almostmaj_or_maj;
		tt local_tt_reminder; 
		std::string local_tt_exact; 
	  
		if (cuts->is_maj_or_almost(node,cut_c) >= 0)
		{
			local_tt_is_maj = cuts->tt_comp_in(node,cut_c);
			local_tt_almostmaj_or_maj = (unsigned) cuts->is_maj_or_almost(node,cut_c);
			local_tt_reminder = cuts->tt_reminder(node,cut_c); 
			local_tt_exact = cuts->tt_exact(node,cut_c);
		}
		else continue; 
			
	    float local_max_level = 0u;

        for ( auto leaf : cut )
        {
          local_max_level = std::max( local_max_level, node_to_level[leaf] );
        }

        if ( local_max_level < best_level )  // MANCA IL CONTrollo suLLA SIZE DELL"EXACT "
        {
          best_level = local_max_level;
          best_cut = cut.address(); 
		 
		  best_tt_is_maj = local_tt_is_maj; 
		  best_tt_is_almostmaj_or_maj = local_tt_almostmaj_or_maj; 
		  best_tt_reminder = local_tt_reminder; 
		  best_exact = local_tt_exact; 
        }
		
      }
      node_to_cut[node] = best_cut; 
	  tt_is_maj[node] = best_tt_is_maj; 
	  tt_is_almostmaj_or_maj[node] = best_tt_is_almostmaj_or_maj; 
	  reminder[node] =  best_tt_reminder; 
	  exact[node] = best_exact; 
	  
	  if ((tt_is_almostmaj_or_maj[node] == 0u ) || (tt_is_almostmaj_or_maj[node] == 1u))
	  {
		node_to_level[node] = best_level + 1u;	
	  }
	  else 
	  {
		node_to_level[node] = best_level + 2u;
	  }
    }
  }
}

void mign_flow_almost_maj_manager::extract_cover()
{
  std::cout << " [i]Extract Cover ..." << std::endl;
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
    cover.add_cut_maj( node, cut, tt_is_maj[node], tt_is_almostmaj_or_maj[node], reminder[node], exact[node]);

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
