#include "mign_rewrite.hpp"

#include <boost/graph/topological_sort.hpp>

#include <core/utils/range_utils.hpp>
#include <core/utils/timer.hpp>

#include <classical/mign/mign_bitmarks.hpp>

namespace cirkit
{

/******************************************************************************
 * Types                                                                      *
 ******************************************************************************/

/******************************************************************************
 * Private functions                                                          *
 ******************************************************************************/
void copy_bitmarks( const mign_graph& mign, mign_node node, mign_graph& mign_new, mign_node node_new )
	{
	  for ( auto i = 0u; i < mign.bitmarks().num_layers(); ++i )
	  {
	    if ( mign.bitmarks().is_marked(node, i) )
	    {
	      mign_new.bitmarks().mark(node_new, i);
	    }
	  }
}

std::map<mign_node, mign_function> init_visited_table_sub( const mign_graph& mign, mign_graph& mign_new, bool keep_bitmarks )
{
  std::map<mign_node, mign_function> old_to_new;
  old_to_new[0u] = mign_new.get_constant( false );
  for ( const auto& pi : mign.inputs() )
  {
    old_to_new[pi.first] = mign_new.create_pi( pi.second );
  }

  if ( keep_bitmarks && mign.bitmarks().num_layers() > 0u )
  {
    mign_new.bitmarks().resize_marks( mign.inputs().size() );
    copy_bitmarks( mign, 0u, mign_new, 0u );

    for ( const auto& pi : mign.inputs() )
    {
      copy_bitmarks( mign, pi.first, mign_new, old_to_new[pi.first].node );
    }
  }

  return old_to_new;
}
	
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

mign_function mign_rewrite_top_down_rec_sub( const mign_graph& mign, mign_node node,
                                       mign_graph& mign_new, std::map<mign_node, mign_function>& old_to_new,
                                       const mign_substitutes_map_t& substitutes,
                                       bool keep_bitmarks )
{
  /* reroute node if it is in substutitutes */
  auto complement = false;
  mign_substitutes_map_t::value_type::const_iterator it_s{};
  if ( substitutes && ( it_s = substitutes->find( node ) ) != substitutes->end() )
  {
    node = it_s->second.node;
    complement = it_s->second.complemented;
  }

  /* visited */
  const auto it = old_to_new.find( node );
  if ( it != old_to_new.end() )
  {
    return it->second ^ complement;
  }

  mign_function f;
 
  const auto c = mign.children( node );
  std::vector<mign_function> operands; 
  for (auto x = 0; x < c.size(); x++)
  {
	  operands.push_back(mign_rewrite_top_down_rec_sub( mign, c[x].node, mign_new, old_to_new, substitutes, keep_bitmarks ) ^ c[x].complemented); 
  }
  
  f = mign_new.create_maj(operands); 

  f.complemented = ( f.complemented != complement ); /* Boolean XOR */
  old_to_new.insert( {node, f} );

  if ( keep_bitmarks && mign.bitmarks().num_layers() > 0u )
  {
    mign_new.bitmarks().resize_marks(f.node);
    copy_bitmarks( mign, node, mign_new, f.node );
  }

  return f;
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

mign_graph mign_rewrite_top_down_sub( const mign_graph& mign, const properties::ptr& settings,
                                const properties::ptr& statistics )
{
  /* settings */
  const auto substitutes   = get( settings, "substitutes",   mign_substitutes_map_t() );
  const auto keep_bitmarks = get( settings, "keep_bitmarks", true );

  /* statistics */
  properties_timer t( statistics );

  mign_graph mign_new; 

  if ( keep_bitmarks && mign.bitmarks().num_layers() > 0u )
  {
    mign_new.bitmarks().init_marks( 0u, mign.bitmarks().num_layers() );
    mign_new.bitmarks().set_used( mign.bitmarks().get_used() );
  }

  /* create constant and PIs */
  auto old_to_new = init_visited_table_sub( mign, mign_new, keep_bitmarks );


  /* map nodes */
  for ( const auto& po : mign.outputs() )
  {
    mign_new.create_po( mign_rewrite_top_down_rec_sub( mign, po.first.node, mign_new, old_to_new, substitutes, keep_bitmarks ) ^ po.first.complemented, po.second );
  }

  return mign_new;
}

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
