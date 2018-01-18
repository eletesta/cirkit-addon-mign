#ifndef FANOUT_FREE_REGIONS_HPP
#define FANOUT_FREE_REGIONS_HPP

#include <map>
#include <queue>
#include <vector>
#include <ios>
#include <string>
#include <fstream>


#include <boost/assign/std/vector.hpp>
#include <boost/bimap.hpp>
#include <boost/format.hpp>
#include <boost/function.hpp>
#include <boost/graph/adjacency_iterator.hpp>
#include <boost/graph/topological_sort.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>
#include <boost/range/algorithm.hpp>
#include <boost/range/iterator_range.hpp>
#include <boost/range/value_type.hpp>
#include <boost/regex.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/find_iterator.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/parameter.hpp> 


#include "core/utils/graph_utils.hpp"
#include "core/utils/range_utils.hpp"
#include "core/utils/timer.hpp"

#include "classical/mign/mign.hpp"

using namespace boost::assign;

namespace cirkit
{

/******************************************************************************
 * Topologically sorted output list                                           *
 ******************************************************************************/

struct topsort_compare_t
{
  explicit topsort_compare_t( mign_graph g )
    : topsortinv( g.num_gates() )
  {
	  auto f = g.topological_nodes(); 
      for ( const auto& v :  f)   {
		  auto pos = find(f.begin(), f.end(), v) - f.begin();
		  topsortinv[v] = pos; 
	   }
  }

  bool operator()( const mign_node& v1, const mign_node& v2 ) const
  {
    return topsortinv.at( v1 ) < topsortinv.at( v2 );
  }

private:
  std::vector<mign_node> topsortinv;
};

using ffr_output_queue_t = std::priority_queue<mign_node, std::vector<mign_node>, topsort_compare_t>;

ffr_output_queue_t make_output_queue( std::vector<mign_node>& outputs, mign_graph& g)
{
  
  /* dequeue for keeping track of FFR outputs */
  topsort_compare_t comp( g );
  return ffr_output_queue_t( outputs.begin(), outputs.end(), comp );
}

/******************************************************************************
 * Relabelling queue                                                          *
 ******************************************************************************/


using relabel_queue_t = std::priority_queue<mign_node>;

using relabel_map_t = std::map<mign_node, boost::bimap<unsigned, mign_node>>;


namespace detail
{

void compute_ffr_inputs_rec( const mign_node& v,
                             const mign_node& ffr_output,
                             ffr_output_queue_t& ffr_outputs,
                             std::vector<mign_node>& ffr_inputs,
							 bool relabel, relabel_queue_t& relabel_queue,
                             mign_graph& g )
{
  /* relabel? */
  if ( relabel )
  {
    relabel_queue.push( v );
  }

  /* primary input? */
  if ( g.is_input(v) )
  {
    if ( boost::find( ffr_inputs, v ) == ffr_inputs.end() )
    {
      ffr_inputs += v;
    }
  }
  /* if ffr output? */
  else if ( v != ffr_output && g.fanout_count(v) > 1u )
  {
    ffr_inputs += v;
    ffr_outputs.push( v );
  }
  else
  {
    for ( auto& adj : g.children(v))
    {
      compute_ffr_inputs_rec( adj.node, ffr_output, ffr_outputs, ffr_inputs,  relabel, relabel_queue, g );
    }
    for ( auto& adj : g.parents(v))
    {
      compute_ffr_inputs_rec( adj, ffr_output, ffr_outputs, ffr_inputs,  relabel, relabel_queue, g );
    }
  }
}

std::vector<mign_node> compute_ffr_inputs( const mign_node& output,
                                                 ffr_output_queue_t& ffr_outputs,
												 bool relabel, relabel_queue_t& relabel_queue,
                                                 mign_graph& g )
{
  std::vector<mign_node> ffr_inputs;
  compute_ffr_inputs_rec( output, output, ffr_outputs, ffr_inputs, relabel, relabel_queue, g );
  return ffr_inputs;
}

}

std::map<mign_node, std::vector<mign_node>> fanout_free_regions( mign_graph mign,
                                                                             const properties::ptr& settings = properties::ptr(),
                                                                             const properties::ptr& statistics = properties::ptr() )
{
  std::map<mign_node, std::vector<mign_node>> result;

  /* settings */
  const auto verbose      = get( settings, "verbose",      false );
        auto outputs      = get( settings, "outputs",      std::vector<mign_node>() );
  const auto relabel      = get( settings, "relabel",      true );
  const auto has_constant = get( settings, "has_constant", false );
  

  /* run-time */
  properties_timer t( statistics );

  /* relabeling */
  relabel_map_t relabel_map;

  /* pre-compute indegrees */
  mign.compute_fanout(); 
  mign.compute_parents(); 

  /* output queue */
  auto ffr_outputs = make_output_queue( outputs, mign);

  /* compute each FFR */
  while ( !ffr_outputs.empty() )
  {
    auto ffr_output = ffr_outputs.top();
    if ( result.find( ffr_output ) == result.end() )
    {
      relabel_queue_t relabel_queue;

      if ( relabel && has_constant )
      {
        relabel_queue.push( 0u );
      }

      auto ffr_inputs = detail::compute_ffr_inputs( ffr_output, ffr_outputs, relabel, relabel_queue, mign );

      if ( verbose )
      {
        std::cout << boost::format( "[i] found ffr region %d(%s)" ) % ffr_output % any_join( ffr_inputs, ", " ) << std::endl;
      }

      result.insert( {ffr_output, ffr_inputs} );

      if ( relabel )
      {
        typename relabel_map_t::mapped_type bm;

        auto pos  = 0u;
        auto last = -1;
        while ( !relabel_queue.empty() )
        {
          auto val = relabel_queue.top();

          if ( (int)val != last )
          {
            bm.insert( typename relabel_map_t::mapped_type::value_type( pos++, val ) );
            last = val;
          }

          relabel_queue.pop();
        }

        relabel_map.insert( {ffr_output, bm} );
      }
    }
    ffr_outputs.pop();
  }

  if ( relabel )
  {
    set( statistics, "relabel_map", relabel_map );
  }

  return result;
}

}

#endif