
#include "classical/mign/mign_utils.hpp"
#include "classical/mign/mign_mffc.hpp"
#include "classical/mign/mign_rewrite.hpp"

#include <classical/mign/mign_simulate.hpp>

#include <classical/mign/minimum_ce/minimum_ce.hpp>

#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/assign/std/vector.hpp>
#include <boost/format.hpp>
#include <boost/graph/topological_sort.hpp>
#include <boost/range/algorithm.hpp>
#include <boost/range/iterator_range.hpp>
#include <boost/regex.hpp>

using namespace boost::assign;

using boost::format;

namespace cirkit {
	
	
mign_function create_mig_rec (mign_graph& mign, mign_node node, std::map<mign_node, mign_function>& old_to_new, mign_graph& mign_new)
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
    	operands.push_back(create_mig_rec( mign, x.node, old_to_new, mign_new ) ^ x.complemented);
    }
  
    const auto f = mign_new.create_maj(operands); 

    old_to_new.insert( {node, f} );
    return f;
}
	
mign_graph create_mig (std::map<mign_node, std::vector<mign_node>> mffc, mign_node node, mign_graph& mign_old,  std::map<mign_node, mign_function>& old_to_new)
{
	mign_graph mign_new;
	
    /* create constant and PIs */
    old_to_new[0] = mign_new.get_constant( false ); 
    for ( const auto& pi : mffc[node] )
    {
	  if (pi == 0) {continue; }
      const auto it = old_to_new.find( pi );
      if ( it != old_to_new.end() ){ continue; }
	  if (mign_old.is_input(pi))
		  old_to_new[pi] = mign_new.create_pi(mign_old.input_name(pi));
	  else
          old_to_new[pi] = mign_new.create_pi("x_" + pi); 
    }
	
	mign_new.create_po(create_mig_rec(mign_old, node, old_to_new, mign_new), "f");
	
	return mign_new; 
}

mign_function rewrite_map_rec (std::unordered_map<mign_node, mign_function>& rewrite_map, mign_graph& mign, std::map<mign_node, mign_function>& old_to_new, mign_node n,std::vector<mign_node> mffc_r, mign_function cc , mign_graph& mign_old)
{
	std::vector<mign_function> operands; 
	
	if (n == 0) {return mign_old.get_constant(false); }
	else if (mign_old.is_input(n)) {
		for (auto & t : mign_old.inputs())
			{
				if (t.first == n)
					return t.first; 
			}
		}
	else 
	{
		for (auto& y : mffc_r)
		{
			if (y == n)
			{
				for (auto & h : mign_old.topological_nodes())
				{
					if (h == y)
						return mign_function(h,0); 
				}
			}	
		}
	}
	
	auto c = mign_old.children(n); 
	auto children = mign.children(old_to_new[n].node); 
		
	for (auto x = 0; x < c.size(); x++)
	    {
			operands.push_back(rewrite_map_rec(rewrite_map, mign, old_to_new, c[x].node, mffc_r, children[x], mign_old) ^ cc.complemented);
			assert (old_to_new[c[x].node].node == children[x].node); 
	    }
		
	auto f = mign_old.create_maj(operands); 
	rewrite_map.insert( std::make_pair( n, f)); 
	
	return f;
}

void save_rewrite_map (std::unordered_map<mign_node, mign_function>& rewrite_map, mign_graph& mign, std::map<mign_node, mign_function>& old_to_new, mign_node root, std::vector<mign_node> mffc_r, mign_graph& mign_old)
{
	
	auto flag = 0; 
	if (root == 0) {flag = 2;}
	else if (mign_old.is_input(root)) {flag = 2;}
	
	if (flag == 0)
	{
		rewrite_map.insert( std::make_pair( root, mign_function(root,mign.outputs()[0].first.complemented))); 
		auto c = mign_old.children(root); 
		auto children = mign.children(old_to_new[root].node); 
		for (auto x = 0; x < c.size(); x++)
	    {
			assert ((old_to_new[c[x].node].node) = (children[x].node));
	    	rewrite_map_rec(rewrite_map, mign, old_to_new, c[x].node, mffc_r, children[x], mign_old); 
	    }
	}
}

/******************************************************************************
	    * Public functions                                                           *
 ******************************************************************************/

mign_graph mffrc_sat ( mign_graph& mign,
	                               const properties::ptr& settings,
	                               const properties::ptr& statistics )
{
		unsigned ce = 0u;
		mign_graph new_mign; 
		std::unordered_map<mign_node, mign_function> rewrite_map;
		
		ce = compute_ce(mign);
		std::cout << "Initial number of ce =\t" << ce << std::endl;
		
		auto mffc = mign_mffcs(mign); 
		
		for ( auto node : mign.topological_nodes() )
		{
		    if ( mffc.find( node ) == mffc.end() ) { continue; }
			else 
			{
				std::map<mign_node, mign_function> old_to_new;
				auto mign_t = create_mig (mffc, node, mign, old_to_new); 
				ce = compute_ce(mign_t);
				set( settings, "start",  ce );
				set( settings, "minimum_after", true );
				set( settings, "verbose", false );
				const mign_tt_simulator simulator{};
				auto values = simulate_mign(mign_t, simulator);
				auto spec_tt = values[mign_t.outputs()[0].first];
				mign_t = min_ce_with_sat(mign_t,spec_tt, settings, statistics); 
				save_rewrite_map(rewrite_map, mign_t, old_to_new, node, mffc[node], mign); 					
			}
		}
		
		settings->set("substitutes", mign_substitutes_map_t( rewrite_map ) );
		
		new_mign = mign_rewrite_top_down_sub(mign,settings,statistics); 
		ce = compute_ce (new_mign);
		std::cout << "Final number of ce =\t" << ce << std::endl; 

	    return new_mign; 
	
}

}