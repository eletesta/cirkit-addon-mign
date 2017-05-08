#include "mign_rewrite.hpp"

#include <boost/graph/topological_sort.hpp>


namespace cirkit
{

/******************************************************************************
 * Types                                                                      *
 ******************************************************************************/

/******************************************************************************
 * Private functions                                                          *
 ******************************************************************************/

std::map<mign_node, mign_function> init_visited_table( const mign_graph& mign, mign_graph& mign_new )
{
  std::map<mign_node, mign_function> old_to_new;

  old_to_new[0] = mign_new.get_constant( false ); 
  for ( const auto& pi : mign.inputs() )
  {
    old_to_new[pi.first] = mign_new.create_pi( pi.second);
  }
 
  return old_to_new;
}


mign_function mign_rewrite_top_down_rec( const mign_graph& mign, mign_node node,
                                       mign_graph& mign_new,
                                       std::map<mign_node, mign_function>& old_to_new )
{
	std::vector<mign_function> operands; 
	
  /* visited */
  const auto it = old_to_new.find( node );
  if ( it != old_to_new.end() )
  {
    return it->second;
  }

  const auto c = mign.children(node); 
  for ( auto& x : c)
  {
  	operands.push_back(mign_rewrite_top_down_rec( mign, x.node, mign_new, old_to_new ) ^ x.complemented);
  }
  //operands.push_back(mign_rewrite_top_down_rec( mign, c[0].node, mign_new, old_to_new ) ^ c[0].complemented);
  //operands.push_back(mign_rewrite_top_down_rec( mign, c[1].node, mign_new, old_to_new ) ^ c[1].complemented); 
  //operands.push_back(mign_rewrite_top_down_rec( mign, c[2].node, mign_new, old_to_new ) ^ c[2].complemented);
  const auto f = mign_new.create_maj(operands); 

  old_to_new.insert( {node, f} );
  return f;
}

mign_function rewrite_max_fO_rec(mign_graph& mign, unsigned int fanout, mign_node node, mign_graph& mign_new, std::map<mign_node, mign_function>& old_to_new)
{
	  std::vector<mign_function> operands; 
	
	  /* visited */
	  const auto it = old_to_new.find( node );
	  if ( it != old_to_new.end() )
	  {
	    return it->second;
	  }

	  
	  const auto c = mign.children(node); 
	  for ( auto& x : c)
	  {
	  	operands.push_back(rewrite_max_fO_rec( mign, fanout, x.node, mign_new, old_to_new ) ^ x.complemented);
	  }
	  
	  const auto f = mign_new.create_maj_fo_restricted(operands, fanout); 

	  //old_to_new.insert( {node, f} );
	  return f;
}

/******************************************************************************
 * Public functions                                                           *
 ******************************************************************************/


mign_graph mign_rewrite_top_down( const mign_graph& mign,
                                const properties::ptr& settings,
                                const properties::ptr& statistics )
{
  /* settings */
 // const auto verbose = get( settings, "verbose", false );

  /* statistics */
  //properties_timer t( statistics );
	
  mign_graph mign_new;

  /* create constant and PIs */
  auto old_to_new = init_visited_table( mign, mign_new );

  /* map nodes */
  for ( const auto& output : mign.outputs())
  {
	 
    mign_new.create_po(mign_rewrite_top_down_rec( mign, output.first.node, mign_new, old_to_new ) ^ output.first.complemented, output.second );
  }

  return mign_new;
}

}
