
//#include <classical/mign/minimum_ce/ce_opt_tree.hpp>

#include <classical/mign/mign_utils.hpp>
#include "ffr_ce_opt.hpp"

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
namespace cirkit
{

/*bool obtain_outcompl_ffr (mign_graph mig, mign_node node)
{
	bool outcompl; 
	auto max_depth = 0u;  
	mig.compute_levels(); 

	max_depth = evaluate_depth(mig);  
	
			for (const auto& node_bis : mig.vertices())
			{
				if ( !boost::out_degree( node_bis, mig ) ) { continue; }

				    const auto children = mig.children(node_bis);
					for ( auto x = 0; x<3; ++x)
					{
						if (children[x].node == node)
						{ 
							outcompl = children[x].complemented; 
						}
					}

				}
		
			if (level(node) == max_depth)
			{ 
				for ( auto x = 0; x<mig.outputs.size(); ++x)
				{
					if (node == mig.outputs[x].first.node)
					{
				outcompl = mig.outputs[x].first.complemented; 
				
			}
			}
			}	
			return outcompl; 
		}


std::vector<std::pair<mign_function,unsigned>> mig_rewrite_ce_ffr( const mign_graph& old_mig, mign_node node,
	                                       mign_graph& new_mig,
	                                       std::map<mign_node, std::vector<std::pair<mign_function,unsigned>>> & old_to_new, bool outcompl )
{
	 
 const auto it = old_to_new.find( node );
 if ( it != old_to_new.end() )
	    {
	      return it->second;
	    }

 	const auto c = old_mig.children(node);
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
/*	
mign_graph ffr_min_ce ( const mign_graph old_mig,
	                               const properties::ptr& settings,
	                               const properties::ptr& statistics )
{
		unsigned ce = 0u;
		mign_graph new_mig; 
	 
		ce = compute_ce(old_mig);
		std::cout << "Initial number of ce =\t" << ce << std::endl;
		
		std::string model_name = "minimumce";
		mig.set_name(model_name); 
		
		auto old_to_new = init_visited_table_ce( old_mig, new_mig);

		auto ffrs = fanout_free_regions(old_mig); 
		//std::vector<mig_node> top( num_vertices( old_mig ) );
		//boost::topological_sort( old_mig, top.begin() );
		
		for ( auto node : mig.topological_nodes() )
		{
		    if ( ffrs.find( node ) == ffrs.end() ) { 
				if (e_output(old_mig,node) == false)
				continue; 
				else 
				{
  		          for ( const auto& output : old_mig.outputs() )
  		          {
  				      if (output.first.node != node) {continue;}
  				      else {
               
  		                    new_mig.create_po(mig_rewrite_ce_ffr( old_mig, output.first.node, new_mig, old_to_new ,output.first.complemented)[0].first ^ output.first.complemented, output.second );
  		                    }
				   }
			    }
			}
			else 
			{
				const auto it = old_to_new.find( node );
				if ( it != old_to_new.end() ) {continue;
				 }
				else {
					 if (e_output(old_mig,node) == false)
					    {
					    const auto outcompl = obtain_outcompl_ffr(old_mig, node); 
					    const auto c = old_mig.children(node);  
					    const auto f = change_ce( new_mig,comple(
						                           mig_rewrite_ce_ffr( old_mig, c[0].node, new_mig, old_to_new, c[0].complemented ),c[0].complemented),
						                           comple(mig_rewrite_ce_ffr( old_mig, c[1].node, new_mig,old_to_new, c[1].complemented ),c[1].complemented),
						                           comple(mig_rewrite_ce_ffr( old_mig, c[2].node, new_mig,old_to_new, c[2].complemented ) ,c[2].complemented), outcompl, old_mig, node); 
				        old_to_new.insert( {node, f} );
			            }
            
			         else if (e_output(old_mig,node) == true)
			            {
                 
		                for ( const auto& output : old_mig.outputs() )
		                    {
				            if (output.first.node != node) {continue;}
				            else {
		                          new_mig.create_po(mig_rewrite_ce_ffr( old_mig, output.first.node, new_mig, old_to_new ,output.first.complemented)[0].first ^ output.first.complemented, output.second );
		                          }
			                }
		                 }
					 } 
              }
        }
     
	    new_mig = mig_strash(new_mig); 
		
		ce = compute_ce (new_mig);
		std::cout << "Final number of ce =\t" << ce << std::endl; 

	    return new_mig; 
	
	}*/
	
mign_graph ffr_min_ce (const mign_graph mig, const properties::ptr& settings,const properties::ptr& statistics )
	{
		mign_graph mign; 
		return mign; 
	}   
}