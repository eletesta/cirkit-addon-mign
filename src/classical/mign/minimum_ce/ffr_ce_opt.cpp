#include "classical/mign/minimum_ce/ffr_ce_opt.hpp" 
#include "classical/mign/minimum_ce/ce_opt_tree.hpp"
#include "classical/mign/minimum_ce/mig_utils_ce.hpp"
#include "classical/mign/minimum_ce/fanout_free_regions.hpp"

#include "classical/mig/mig.hpp"
#include "classical/mig/mig_utils.hpp"

#include "classical/mign/mign_utils.hpp"
#include "classical/mign/mig_to_mign.hpp"

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
	
bool obtain_outcompl_ffr (mig_graph mig, mig_node node)
{
	bool outcompl; 
	auto max_depth = 0u; 
	const auto& info = mig_info( mig );
	std::ostream& os = std::cout; 

	max_depth = mig_print_stats_depth(mig,os); 
	//std::cout << "siamo in outcompl e siamo al nodo " << node << std::endl; 
			for (const auto& node_bis : boost::make_iterator_range(boost::vertices(mig)))
			{
				if ( !boost::out_degree( node_bis, mig ) ) { continue; }

				    const auto children = get_children( mig, node_bis );
					for ( auto x = 0; x<3; ++x)
					{
						if (children[x].node == node)
						{ 
							outcompl = children[x].complemented; 
						}
					}

				}
		
			if (compute_level(mig,node) == max_depth)
			{ 
				for ( auto x = 0; x<info.outputs.size(); ++x)
				{
					if (node == info.outputs[x].first.node)
					{
				outcompl = info.outputs[x].first.complemented; 
				
			}
			}
			}	
			return outcompl; 
		}


std::vector<std::pair<mig_function,unsigned>> mig_rewrite_ce_ffr( const mig_graph& old_mig, mig_node node,
	                                       mig_graph& new_mig,
	                                       std::map<mig_node, std::vector<std::pair<mig_function,unsigned>>> & old_to_new, bool outcompl )
{
	 
 const auto it = old_to_new.find( node );
 if ( it != old_to_new.end() )
	    {
	      return it->second;
	    }

 	const auto c = get_children( old_mig, node );
	const auto f = change_ce( new_mig,comple(
	                           mig_rewrite_ce_ffr( old_mig, c[0].node, new_mig, old_to_new, c[0].complemented ),c[0].complemented),
	                           comple(mig_rewrite_ce_ffr( old_mig, c[1].node, new_mig,old_to_new, c[1].complemented ),c[1].complemented),
	                           comple(mig_rewrite_ce_ffr( old_mig, c[2].node, new_mig,old_to_new, c[2].complemented ) ,c[2].complemented), outcompl, old_mig, node);

	    old_to_new.insert( {node, f} );
	    return f;
 

	}
/******************************************************************************
	    * Public functions                                                           *
 ******************************************************************************/

mign_graph ffr_min_ce ( mign_graph& mign,
	                               const properties::ptr& settings,
	                               const properties::ptr& statistics )
{
		unsigned ce = 0u;
		mign_graph new_mign; 
		mig_graph new_mig; 
		ce = compute_ce(mign);
		std::cout << "Initial number of ce =\t" << ce << std::endl;
		
		auto old_mig = mign_to_mig (mign); 
	 
		std::string model_name = "minimumce";
		mig_initialize(new_mig,model_name);
		
		auto old_to_new = init_visited_table_ce( old_mig, new_mig);

		auto ffrs = fanout_free_regions(old_mig); 
		std::vector<mig_node> top( num_vertices( old_mig ) );
		boost::topological_sort( old_mig, top.begin() );
		
		for ( auto node : top )
		{
			//std::cout << " siamo al nodo = per le ffr " << node << std::endl; 
		    if ( ffrs.find( node ) == ffrs.end() ) { 
				if (e_output(old_mig,node) == false)
				continue; 
				else 
				{
  		          for ( const auto& output : mig_info( old_mig ).outputs )
  		          {
  				  if (output.first.node != node) {continue;}
  				  else {
                      //std::cout << " creiamo l'output che non e ffr" << std::endl;
  		         mig_create_po( new_mig, mig_rewrite_ce_ffr( old_mig, output.first.node, new_mig, old_to_new ,output.first.complemented)[0].first ^ output.first.complemented, output.second );
  		          }
				}
			}
				//mig_create_po( new_mig, mig_rewrite_ce_ffr( old_mig, output.first.node, new_mig, old_to_new ,output.first.complemented)[0].first ^ output.first.complemented, output.second );
				//std::cout << " nodo tot NON e nelle ffr" << std::endl; 
			}
			else 
			{
				const auto it = old_to_new.find( node );
				if ( it != old_to_new.end() ) {continue;
				//std::cout << " nodo tot e nelle ffr ma e gia statp messo" << std::endl;
				 }
				else {
					if (e_output(old_mig,node) == false)
					{
						//std::cout << " nodo tot e nelle ffr e non e un output" << std::endl; 
					const auto outcompl = obtain_outcompl_ffr(old_mig, node); 
					const auto c = get_children(old_mig,node); 
					const auto f = change_ce( new_mig,comple(
						                           mig_rewrite_ce_ffr( old_mig, c[0].node, new_mig, old_to_new, c[0].complemented ),c[0].complemented),
						                           comple(mig_rewrite_ce_ffr( old_mig, c[1].node, new_mig,old_to_new, c[1].complemented ),c[1].complemented),
						                           comple(mig_rewrite_ce_ffr( old_mig, c[2].node, new_mig,old_to_new, c[2].complemented ) ,c[2].complemented), outcompl, old_mig, node); 
				 old_to_new.insert( {node, f} );
			 }
            
			//std::cout << " questo deve farlo tutte le volte" << std::endl;
			else if (e_output(old_mig,node) == true)
			  {
                  // std::cout << " creiamo l'output che e ffr" << std::endl;
				  //std::cout << " nodo tot e nelle ffr ed e un output" << std::endl;
		          for ( const auto& output : mig_info( old_mig ).outputs )
		          {
				  if (output.first.node != node) {continue;}
				  else {
		         mig_create_po( new_mig, mig_rewrite_ce_ffr( old_mig, output.first.node, new_mig, old_to_new ,output.first.complemented)[0].first ^ output.first.complemented, output.second );
		          }
			  }
		  }
					 
                } }
        }

	    new_mig = mig_strash(new_mig); 
		
	    new_mign = mig_to_mign(new_mig);  
		
		ce = compute_ce (new_mign);
		std::cout << "Final number of ce =\t" << ce << std::endl; 

	    return new_mign; 
	
}

}