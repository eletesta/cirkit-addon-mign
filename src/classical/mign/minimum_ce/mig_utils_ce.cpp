#include "mig_utils_ce.hpp"

#include <boost/assign/std/vector.hpp>
#include <boost/format.hpp>
#include <boost/range/algorithm.hpp>
#include <boost/range/iterator_range.hpp>

#include "core/properties.hpp"
#include "core/graph/depth.hpp"
#include "core/utils/graph_utils.hpp"

#include "classical/mig/mig_utils.hpp"

//#include <simulate_mig.hpp>

using namespace boost::assign;

namespace cirkit
{

/******************************************************************************
 * Types                                                                      *
 ******************************************************************************/

/******************************************************************************
 * Private functions                                                          *
 ******************************************************************************/

/******************************************************************************
 * Public functions                                                           *
 ******************************************************************************/

unsigned mig_print_stats_depth( const mig_graph& mig, std::ostream& os )
{
  const auto& info = mig_info( mig );
  auto n = info.inputs.size();

  std::string name = info.model_name;
  if ( name.empty() )
  {
    name = "(unnamed)";
  }

  std::vector<mig_node> outputs;
  for ( const auto& output : info.outputs )
  {
    outputs += output.first.node;
  }

  std::vector<unsigned> depths;
  const auto depth = compute_depth( mig, outputs, depths );

  //os << boost::format( "[i] %20s: i/o = %7d / %7d  maj = %7d  lev = %4d" ) % name % n % info.outputs.size() % ( boost::num_vertices( mig ) - n - 1u ) % depth << std::endl;
    
    return depth;
}

unsigned depth_mig( const mig_graph& mig)
{
  const auto& info = mig_info( mig );
  //auto n = info.inputs.size();

  std::vector<mig_node> outputs;
  for ( const auto& output : info.outputs )
  {
    outputs += output.first.node;
  }

  std::vector<unsigned> depths;
  const auto depth = compute_depth( mig, outputs, depths );

//  os << boost::format( "[i] %20s: i/o = %7d / %7d  maj = %7d  lev = %4d" ) % name % n % info.outputs.size() % ( boost::num_vertices( mig ) - n - 1u ) % depth << std::endl;
    
    return depth;
}


std::vector<mig_function> get_children( const mig_graph& mig, const mig_node& node )
{
  std::vector<mig_function> children;

  for ( const auto& edge : boost::make_iterator_range( boost::out_edges( node, mig ) ) )
  {
    children += mig_to_function( mig, edge );
  }

  return children;
}
    
unsigned compute_level( const mig_graph& mig, const mig_node& node)
    {
        const auto& info = mig_info( mig );
       // auto n = info.inputs.size();
        
        std::string name = info.model_name;
        if ( name.empty() )
        {
            name = "(unnamed)";
        }
        
        std::vector<mig_node> outputs;
        for ( const auto& output : info.outputs )
        {
            outputs += output.first.node;
        }

        std::vector<unsigned> depths;
        compute_depth( mig, outputs, depths );
       
        
        return depths[node];
    }

    
}