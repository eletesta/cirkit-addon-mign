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

#include "mign_invfree.hpp"

#include <iostream>
#include <fstream>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <ctime>

#include <classical/mign/mign_utils.hpp>

namespace cirkit
{
	

void write_path_rec( mign_graph& mign, mign_node node, std::vector<std::tuple<mign_node, std::vector<mign_node>, std::string>>& write_output, std::vector<mign_node>& nodes, mign_node output)
	{
		//std::cout << " path rec" << std::endl; 
		auto c = mign.children(node); 
		for (auto& x : c)
		{
			//std::cout << " siamo al nodo " << x.node << std::endl; 
			if (x.node == 0) {continue;}
			if (mign.is_input(x.node)) {write_output.push_back(std::make_tuple(output, nodes, mign.input_name(x.node)));}
			else 
			{
				auto prova = nodes;  
				prova.push_back(x.node); 
				write_path_rec(mign, x.node,write_output,prova,output); 
			}
			
		}
	}
	
void write_paths (mign_graph& mign)
	{ 
		mign.compute_fanout();  
		std::vector<std::tuple<mign_node, std::vector<mign_node>, std::string>> write_output; // output -> nodes -> inputs
		std::vector<mign_node> nothing(0); 
	
	
		auto count = 0; 
	    for ( auto& output : mign.outputs())
	    {
			auto c = mign.children(output.first.node); 
			for (auto& x : c)
			{
				//std::cout << " siamo al nodo " << x.node << std::endl; 
				if (x.node == 0) {continue; }
				if (mign.is_input(x.node)) {write_output.push_back(std::make_tuple(output.first.node, nothing, mign.input_name(x.node)));}
				else 
				{
					std::vector<mign_node> nodes; 
					nodes.push_back(x.node); 
					write_path_rec(mign, x.node,write_output,nodes,output.first.node); 
				}
				
			}
			//write_output[count] = std::make_tuple(output.first.node, obtain_nodes(mign,output.first.node, nodes), obtain_input(mign,output.first.node, inputs)); 	 
		}
		
		//std::cout << "nome = " << mign.name(); 
		std::ofstream Savefile (mign.name()+".txt"); 
		
		for (auto& x : write_output)
			{
				//if (mign.is_input(get<2>(x)))
				Savefile << " [i] Input " <<  get<2>(x) << ": " ; 
				
				if (get<1>(x).size() > 1)
					std::reverse(get<1>(x).begin(),get<1>(x).end()); 
				
				for ( auto& j : get<1>(x))
				{
					Savefile << " -> node w" << j - mign.inputs().size() - 1 << " ( FO = " << mign.fanout_count(j) << ")"; 
				}
				Savefile << " -> node w" << get<0>(x)-mign.inputs().size() - 1<< " = Output : " << mign.output_name( get<0>(x) ) << std::endl; 
			}  
	}

unsigned ce_count (const mign_graph& mign)
{
	unsigned ce = 0u; 
	unsigned count = 0; 
	for ( auto& node : mign.topological_nodes())
	{
		
	     if ( mign.is_input( node) ) { continue; }
		 count++; 

	      const auto c = mign.children( node);
	
	  	for (auto x :c)
	  	{
			
	  		if (( x.node != 0) && (x.complemented == 1))
			{
				++ce; 
			}
	  	}	
	}
	
	for ( const auto& output : mign.outputs())
		if ( output.first.complemented == 1)
			++ce; 
		
	std::cout << " size = " << count << std::endl; 
	std::cout << " Number of ce = " << ce << std::endl; 
	return ce; 
}

std::map<mign_node, mign_function> init_visited_table_inv( const mign_graph& mign, mign_graph& mign_new )
{
  std::map<mign_node, mign_function> old_to_new;

  old_to_new[0] = mign_new.get_constant( false ); 
  for ( const auto& pi : mign.inputs() )
  {
    old_to_new[pi.first] = mign_new.create_pi( pi.second);
  }
 
  return old_to_new;
}

mign_function mign_invfree_rec(  mign_graph& mign, mign_node node,
                                       mign_graph& mign_new,
                                       std::map<mign_node, mign_function>& old_to_new, int changed, std::vector<int>& visited, std::map<std::vector<mign_function>,mign_function>& saved_nodes)
{
	std::vector<mign_function> operands; 
	
	//std::cout << " nodo " << node << std::endl; 
	//mign_new.disable_strashing(); //TODO and check if it works ;) 
	
  /* visited are only PI here */
  const auto it = old_to_new.find( node );
  if ( it != old_to_new.end() )
  {
    return it->second;
  }
  
  if (mign_new.size() < mign.size() * 10)
 {
	 auto iter = find(visited.begin(), visited.end(),node); 
     if ( iter == visited.end() ) 
     {
		 visited.push_back(node); 
	 }
	 
	//  std::cout << " node = " << node << std::endl; 
 const auto c = mign.children(node); 
 
// std::cout << " node" << node << std::endl; 
 //
 auto p = c; 
 auto count = 0; 
 for (auto j : p) 
 {
	// std::cout << " figlio " << j.node; 

	  if (changed == 1)
	  	{ 
			p[count].complemented = !p[count].complemented; 
	    } 
		
  //std::cout << " con c " << p[count].complemented << std::endl; 
	count++; 	
 }
 
 const auto iter_tw = saved_nodes.find( p );
 if ( iter_tw != saved_nodes.end() )
 {
   return iter_tw->second;
 }
 
  for ( auto& x : c)
  {
	  if ( mign.is_input(x.node))
	  {
		 
		  if (changed == 1)
		  {
		  	operands.push_back(mign_invfree_rec( mign, x.node, mign_new, old_to_new, 0, visited, saved_nodes ) ^ !x.complemented);
		  }
		  else 
			  operands.push_back(mign_invfree_rec( mign, x.node, mign_new, old_to_new, 0 , visited, saved_nodes) ^ x.complemented);
	  }
	  else 
	  {
		  if ( changed == 1)
		 	 {
			 
				 if (x.complemented == 1)
	  			 operands.push_back(mign_invfree_rec( mign, x.node, mign_new, old_to_new, 0 ,visited, saved_nodes) ^ !x.complemented);
		 		 else 
	             operands.push_back(mign_invfree_rec( mign, x.node, mign_new, old_to_new, 1 ,visited, saved_nodes) ^ x.complemented);  
			 }
		else 
		     {
		         if (x.complemented == 0)
		         operands.push_back(mign_invfree_rec( mign, x.node, mign_new, old_to_new, 0 ,visited, saved_nodes) ^ x.complemented);
			     else 
	             operands.push_back(mign_invfree_rec( mign, x.node, mign_new, old_to_new, 1,visited , saved_nodes) ^ !x.complemented);  
	         }
	  }
	  
  }

  const auto f = mign_new.create_maj(operands); 
  saved_nodes.insert({p,f}); 
  return f;
}
else 
{
	float q = mign.size()-mign.inputs().size() - 1 ; 
	float percentage = (q - visited.size())/q * 100;  
	std::cout << " too big. The circuit has a size 10 times bigger than the original one. I have visited " << visited.size() << "/" << mign.size()-mign.inputs().size() - 1 << " = " << percentage << " %" << std::endl; 
	exit(EXIT_FAILURE); 
}

  
}

int final_check(const mign_graph& mign_new, mign_function func, unsigned& flag)
{
	
	if (mign_new.is_input(func.node)) return 0; 
	
	if (flag == 0)
		return 1; 
	
    const auto c = mign_new.children(func.node); 
    for ( auto& x : c)
    {
		if (mign_new.is_input(x.node)) continue; 
			if ( x.complemented == 1)
		{
			flag = 0; 
			return 1; 
		}
		std::cout << " final check per il nodo " << x.node << std::endl; 
		final_check(mign_new,x,flag); 
	}
	return 0; 
}

mign_graph mign_invfree ( mign_graph& mign)
{
	
	mign_graph mign_new; 
	
	std::string fo_res = "inv_free"; 
	
	if (!mign.name().empty())
	mign_new.set_name(mign.name()+fo_res); 
	else 
		mign.set_name("inv_free"); 
	
	mign.compute_fanout(); 
	auto ce = ce_count (mign); 
	
	
  /* create constant and PIs */
     auto old_to_new = init_visited_table_inv( mign, mign_new );
	 
	 std::map<std::vector<mign_function>,mign_function> save_nodes; //to map nodes and the chidlren 

	 std::vector<int> visited; 
  /* map nodes */
	 int start_s=clock();
    for ( const auto& output : mign.outputs())
    {
		if ((mign.is_input(output.first.node)) || (output.first.node == 0))
			mign_new.create_po(mign_invfree_rec(mign, output.first.node, mign_new, old_to_new, 0, visited, save_nodes) ^ output.first.complemented, output.second); 
		else{
			if ( output.first.complemented == 0)
			{
	
				mign_new.create_po(mign_invfree_rec(mign, output.first.node, mign_new, old_to_new, 0, visited, save_nodes) ^ output.first.complemented, output.second); 
			}
			else 
			{
				mign_new.create_po(mign_invfree_rec(mign, output.first.node, mign_new, old_to_new, 1, visited, save_nodes) ^ !output.first.complemented, output.second); 
			}}
		
	 }
	 int stop_s=clock();
	 std::cout << "Inv_free time = " << (stop_s-start_s)/double(CLOCKS_PER_SEC) << std::endl;
		
	 //std::cout << " Final Check: " << std::endl; 
	 
	 /*unsigned flag = 1; 
     for ( const auto& output : mign_new.outputs())
     {
		    std::cout << " output " << output.first.node << std::endl; 
 			if ( output.first.complemented == 0)
 			{
 				auto g = final_check(mign_new, output.first, flag); 
				if ( g == 1)
				{
						std::cout << " Error : output " << output.first.node << " is complemented" << std::endl; 
					   break; 
				   }
 			}
 			else 
 			{
				flag = 0; 
 				std::cout << " Error : output " << output.first.node << " is complemented" << std::endl; 
				break; 
 			}
		
 	}
	
	if (flag == 1)
	{
		std::cout << " Complemented edges:  "; 
		std::cout << " CE ok " << std::endl; 
	}*/
 	
	ce = ce_count (mign_new); 
	
	//write_paths (mign_new); 
	// non basta controllare che non ci siano nodi uguali?? Fa davvero quello che voglio?? 
	
	return mign_new; 
}

mign_function mign_fo_restricted_rec(mign_graph& mign, mign_node node, mign_graph& mign_new, std::map<mign_node, mign_function>& old_to_new)
{
	  std::vector<mign_function> operands; 
	
	  /* visited */
	  const auto it = old_to_new.find( node );
	  if ( it != old_to_new.end() )
	  {
	    return it->second;
	  }

	  std::map<mign_function, mign_function> equivalence; 
	  std::cout << " nel vecchio questo e il nodo = " << node - mign.inputs().size() - 1 << std::endl; 
	  const auto c = mign.children(node); 
	  for ( auto& x : c)
	  {
	  	operands.push_back(mign_fo_restricted_rec( mign, x.node, mign_new, old_to_new ) ^ x.complemented);
	  }
	  
	  const auto f = mign_new.create_maj_fo_restricted(operands,2); 

	  //old_to_new.insert( {node, f} );
	  return f;
}

mign_graph mign_fo_restricted ( mign_graph& mign)
{
	
	mign_graph mign_new; 
	
	std::string fo_res = "fo_rest"; 
	
	if (!mign.name().empty())
	mign_new.set_name(mign.name()+fo_res); 
	else 
		mign_new.set_name("fo_res"); 
	
	mign.compute_fanout(); 	
	
	std::cout << " theory : ";
	
	auto total = 0u; 
 	for ( auto& node : mign.topological_nodes())
 	{
		// std::cout << " siamo al nodo " << node << std::endl; 
		if(mign.is_input(node) || node == 0) {continue; }
		int a = mign.fanout_count(node) / 3; 
		int b = a - 1; 
		if (4%3 != 0)
			b++; 
		total = total + b; 
		if (mign.fanout_count(node) > 3)
			std::cout << " nodo " << node - mign.inputs().size() - 1 << " ha FO > = " << mign.fanout_count(node) << std::endl; 
	}
	
	std::cout << total << std::endl; 
	
  /* create constant and PIs */
     auto old_to_new = init_visited_table_inv( mign, mign_new );

  /* map nodes */
    for ( const auto& output : mign.outputs())
    {
		mign_new.create_po(mign_fo_restricted_rec(mign, output.first.node, mign_new, old_to_new) ^ output.first.complemented, output.second); 
	 }
	
	 std::cout << " New size = " << mign_new.size() - mign_new.inputs().size() - 1 << std::endl; 
	 
	 auto flag = 1; 
	 std::cout << " Final Check: " << std::endl; 
	 mign_new.compute_fanout(); 
	 
 	for ( auto& node : mign_new.topological_nodes())
 	{
		if(mign_new.is_input(node)) {continue; }
		if (mign_new.fanout_count(node) > 3)
			flag = 0; 
	}
	
	if (flag == 1)
	{
		//std::cout << " C:  "; 
		std::cout << " FO ok " << std::endl; 
	}
	
	write_paths (mign_new);  
	return mign_new; 
}


/************************** -- optimized method fanout *************************/





void map_critpath_node_rec (mign_graph&mign, std::map<mign_node,unsigned>& cp_node, mign_node node, unsigned depth)
{
	
      const auto it = cp_node.find( node );
      if ( it != cp_node.end() )
      {
		  if (cp_node[node] < depth)
		  {
			  cp_node.erase(it); 
			  cp_node[node] = depth; 
		  }
      }
	  else 
		  cp_node[node] = depth; 
	  
	  if (!mign.is_input(node))
	  {
	  	auto c = mign.children(node); 
		for (auto& x : c)
			map_critpath_node_rec (mign, cp_node, x.node, depth); 
	  }
}

void map_critpath_node (mign_graph& mign, std::map<mign_node,unsigned>& cp_node)
{
	for ( const auto& output : mign.outputs())
	{
		auto depth = mign.level(output.first.node); 
		cp_node[output.first.node] = depth; 
		map_critpath_node_rec (mign, cp_node, output.first.node, depth); 
	}	
}

unsigned number_of_buffer(mign_graph& mign)
{
	unsigned count = 0u; 
 	for ( auto& node : mign.topological_nodes())
 	{
		if(mign.is_input(node)) {continue; }
		auto c = mign.children(node); 
		if ((c[0].node == 0) && (c[1].node == 0))
		{
			count++; 
		}
		else if ((c[1].node == 0) && (c[2].node == 0))
		{
			count++; 
		}
		else if ((c[0].node == 0) && (c[2].node == 0))
		{
			count++; 
		}
	}
	return count; 
} 

unsigned n_parents_on_cp (mign_graph& mign, std::map<mign_node,unsigned>& cp_node, mign_node node, unsigned depth, unsigned delay)
{
	assert (depth >= delay);
	auto count = 0u; 
	auto p = mign.parents(node); 
	
	for (auto x = 0; x <p.size(); x++)
	{
		if (cp_node[p[x]] == depth - delay)
			count++; 
	}
	return count; 
}

bool other_parent_not_buffer(mign_graph& mign_new, mign_function func, mign_node node)
{
	std::cout << " entri in parents not buffer" << std::endl; 
	auto p = mign_new.parents(func.node); 
	assert (p.size() == 2); 
	
	for (auto & x : p)
	{
		if (x != node)
		{
			auto operands = mign_new.children(x); 
			if (((operands[0] == 0) && (operands[1] == 1)) || ((operands[0] == 1) && (operands[1] == 0)))
				return false; 
			else 
				return true; 
				
		}
	}
	return false; 
}
	
bool children_small_fo (mign_graph& mign_new,std::map<std::vector<mign_function>,unsigned>& count,std::vector<mign_function> children, unsigned fanout)
{
	int div = fanout/3; 
	if (fanout%3 != 0)
		div++; 
	
	for (auto& x :children)
	{
		if (mign_new.is_input(x.node) || (x.node == 0))
		{
			continue; 
		}
		auto xc = mign_new.children(x.node); 
		if (count[xc] + div <=2)
			continue; 
		else return false; 
	}
	return true; 
}
	
mign_function clean_mig_rec(mign_graph& mign_new, mign_node node, std::map<mign_node, mign_function>& old_to_new, std::map<mign_node,mign_function>& old_to_new_fo,std::map<std::vector<mign_function>,unsigned>& count, mign_graph& mign)
{
  
	  //std::cout << " node = " << node << std::endl; 
	  const auto it = old_to_new.find( node );
	  if ( it != old_to_new.end() )
	  {
	    return it->second;
	  }
	  
	   auto c = mign_new.children(node); 
	 
 
	  const auto iter = old_to_new_fo.find(node);
	  
	  if ( iter != old_to_new_fo.end() )
	  {
		  //std::cout << " esiste gia " << std::endl; 
				    auto children = mign.children(old_to_new_fo[node].node);
					if (count[children] <2)
				   {
					   count[children]++; 
					   //const auto f = mign_new.create_maj_fo_restricted(children, max_fo-1);
					   return iter->second; 
				   } 
					else 
					{
						//const auto f = mign_new.create_maj(children);
						//count[children] = 0; 
						//old_to_new[node] = f; 
					    old_to_new_fo.erase(iter);
						//return f; 
						 
					}
	  }
	  
	  //std::cout << " nuovo " << std::endl; 
	  std::vector<mign_function> operands; 
	  
	  for ( auto& x : c)
	  {
	  	operands.push_back(clean_mig_rec(mign_new, x.node, old_to_new,  old_to_new_fo, count, mign) ^ x.complemented);
	  }
	  
	  
	  if (mign_new.fanout_count(node) == 1)
	  {
		  if ((operands[0].node == 0) && (operands[1].node == 0))
		  {
		  	if ((operands[0].complemented == 0) && (operands[1].complemented == 1))
			{
				return operands[2]; 
			}
			else if ((operands[0].complemented == 1) && (operands[1].complemented == 0))
			{
				return operands[2]; 
			}
			else if ((operands[0].complemented == 0) && (operands[1].complemented == 0))
					return mign.get_constant( true);  
			else 
				    return mign.get_constant( false );  
		  }
	  }
	  
	  /*if ((operands[0].node == 0) && (operands[1].node == 0))
	  {
	      
		  auto n = (int) mign_new.fanout_count(c[2].node); 
		  int free_spaces = (3 - n) + 1; 
		  if (mign_new.fanout_count(node) <= free_spaces)
			 {
				 return operands[2]; 
			 } 
	  }*/
	  
		    const auto f = mign.create_maj_10(operands);
			count[operands] = 0; 
			old_to_new_fo[node] = f; 
		    return f;

}

mign_graph clean_mig(mign_graph& mign_new)
{
	//std::cout << " clean mig" << std::endl; 
	mign_graph mign; 
	mign.structural_hashing(1);
	mign_new.compute_fanout(); 
	mign_new.compute_parents(); 
	
	std::map<mign_node,mign_function> old_to_new_fo; 
	std::map<std::vector<mign_function>,unsigned> count;
	
	auto old_to_new = init_visited_table_inv( mign_new, mign );
	
    for ( const auto& output : mign_new.outputs())
    {
		mign.create_po(clean_mig_rec(mign_new, output.first.node, old_to_new,  old_to_new_fo, count, mign) ^ output.first.complemented, output.second); 
	}
	
	return mign; 
}

unsigned father_on_cp (mign_graph& mign,std::map<mign_node, unsigned>& cp_node, unsigned depth, mign_node node)
{
	auto count = 0u; 
	auto p = mign.parents(node); 
	for (auto& pp : p)
	{
		if (cp_node[pp] == depth)
			count++; 
	}
	return count; 
}

void update(std::map<mign_node, unsigned>& cp_node, mign_node father,mign_graph& mign)
{
	if ((mign.is_output(father)) && (mign.fanout_count(father) == 1))
		cp_node[father]++;
	else 
	{
		cp_node[father]++; 
		auto c = mign.parents(father); 
		for (auto x = 0; x < c.size(); x++)
		{
			update(cp_node,c[x],mign); 
		}
	}
	
}

mign_function mign_fo_restricted_rec_opt_depthconst(mign_graph& mign, mign_node node, mign_graph& mign_new,  std::map<mign_node, mign_function>& old_to_new, std::map<mign_node,mign_function>& old_to_new_fo, std::map<mign_node, unsigned>& cp_node, unsigned int max_fo,std::map<std::vector<mign_function>,unsigned>& count, mign_node father, std::map<mign_node, mign_function>& buffer, std::map<mign_function, unsigned>& buffer_count)
{
	 
	  //const auto depth = evaluate_depth(mign); 
	
	  std::vector<mign_function> operands;
	
	  /* visited */
	  const auto it = old_to_new.find( node );
	  if ( it != old_to_new.end() )
	  {
	    return it->second;
	  }
	   
	  auto c = mign.children(node); 
 
	  const auto iter = old_to_new_fo.find(node);
	  //const auto iterpp = old_to_new_fo.find(pp);
	  //std::cout << "nodo = " << node << std::endl; 
	  if ( iter != old_to_new_fo.end() )
	  {
				    auto children = mign_new.children(old_to_new_fo[node].node);
					if (count[children] <2)
				   {
					   count[children]++; 
					   //const auto f = mign_new.create_maj_fo_restricted(children, max_fo-1);
					   return iter->second; 
				   } 
					else 
					{
						//const auto f = mign_new.create_maj(children);
						//count[children] = 0; 
						//old_to_new[node] = f; 
					    old_to_new_fo.erase(iter);
						//return f; 
						 
					}
	  }
	  
	  
	  for ( auto& x : c)
	  {
	  	operands.push_back(mign_fo_restricted_rec_opt_depthconst( mign, x.node, mign_new, old_to_new,old_to_new_fo, cp_node,max_fo , count, node, buffer, buffer_count) ^ x.complemented);
	  }
	  
		    const auto f = mign_new.create_maj(operands);
			count[operands] = 0; 
			old_to_new_fo[node] = f; 
		    return f;
}

mign_function mign_fo_restricted_rec_opt(mign_graph& mign, mign_node node, mign_graph& mign_new,  std::map<mign_node, mign_function>& old_to_new, std::map<mign_node,mign_function>& old_to_new_fo, std::map<mign_node, unsigned>& cp_node, unsigned int max_fo,std::map<std::vector<mign_function>,unsigned>& count, mign_node father, std::map<mign_node, mign_function>& buffer, std::map<mign_function, unsigned>& buffer_count)
{
	 
	  const auto depth = evaluate_depth(mign); 
	
	  std::vector<mign_function> operands;
	
	  /* visited */
	  const auto it = old_to_new.find( node );
	  if ( it != old_to_new.end() )
	  {
	    return it->second;
	  }
	   
	  auto c = mign.children(node); 
 
	  const auto iter = old_to_new_fo.find(node);

	  if ( iter != old_to_new_fo.end() )
	  {
		  auto children = mign_new.children(old_to_new_fo[node].node);
		  if (mign.fanout_count(node) > 3)   // mign.fanout > 3
		  {
			 if (cp_node[father] >= depth)
			 {
				if (count[children] <2)
			   {
				   count[children]++; 
				   return iter->second; 
			   } 
				else 
				{
					old_to_new_fo.erase(iter);
				}
					
			 }
			else if ((father_on_cp(mign,cp_node,depth,node)%3 == 0) && ((mign.fanout_count(node)- father_on_cp(mign,cp_node,depth,node)) <=3))
			 {
 				if (count[children] <2)
 			   {
 				   count[children]++; 
 				   return iter->second; 
 			   } 
 				else 
 				{
 					old_to_new_fo.erase(iter);
 				}
			 }
			 else   // non on the critical path 
		    {
				if (children_small_fo(mign_new,count,children,mign.fanout_count(node)))
				{
					if (count[children] <2)
				   {
					   count[children]++; 
					   return iter->second; 
				   } 
					else 
					{
						old_to_new_fo.erase(iter);
					}
				}
				else   // buffer is inserted only if is not on the critical depth, if the number of children have FO large
				{
	    			  const auto a = buffer.find(node); 
	    			  if (( a != buffer.end() ) && (buffer_count[node] < 3))
	    			  {
	  				  auto g = a->second; 
	    				  buffer_count[node]++; 
	    				  update(cp_node, father, mign); 
	    				  return g; 
	    			  }
	    			  else 
	    			  {
	    				  if ( a != buffer.end())  
	    				  {
	  					  buffer.erase(a); 
	  				      }
	  					  if (count[children] < 2)
	  					  {
	  					  	  count[children]++;
	  	   					  update(cp_node, father, mign); 
	  	   	  				  std::vector<mign_function> operands_two;
	  	   	  				  operands_two.push_back(mign_new.get_constant( false ));  
	  	   	  				  operands_two.push_back(mign_new.get_constant( true ));  
	  	   	  				  operands_two.push_back(old_to_new_fo[node]); 
	  		  				  const auto g = mign_new.create_maj_10(operands_two);
	  						  buffer[node] = g; 
	  						  buffer_count[node] = 1; 
	  		  				  return g; 
	  					  }
	  					  else 
	  					  {
	    						old_to_new_fo.erase(iter);
	  					  }

	  		          }  
				  }
  			  
	 	      }
	      }
		  else   // mign.fanout <= 3 in old graph. Maybe not true anymore. 
		  {
			if (count[children] <2)
		   {
			   count[children]++; 
			   return iter->second; 
		   } 
			else 
			{
				old_to_new_fo.erase(iter);
			}
			   
		  }
	  }
	  
	  for ( auto& x : c)
	  {
	  	operands.push_back(mign_fo_restricted_rec_opt( mign, x.node, mign_new, old_to_new,old_to_new_fo, cp_node,max_fo , count, node, buffer, buffer_count) ^ x.complemented);
	  }
	 
		if (cp_node[father] >= depth)
		{
			
		    const auto f = mign_new.create_maj_10(operands);
			count[operands] = 0; 
			old_to_new_fo[node] = f; 
		    return f;
		}
		else if ((father_on_cp(mign,cp_node,depth,node)%3 == 0) && ((mign.fanout_count(node)- father_on_cp(mign,cp_node,depth,node)) <=3))
		{
		    const auto f = mign_new.create_maj_10(operands);
			count[operands] = 0; 
			old_to_new_fo[node] = f; 
		    return f;
		}
		else 
		{
			if (mign.fanout_count(node) <= 3)
			{
			    const auto f = mign_new.create_maj_10(operands);
				count[operands] = 0; 
				old_to_new_fo[node] = f; 
			    return f;
			}
			else 
			{
				if (children_small_fo(mign_new,count,operands,mign.fanout_count(node)))
				{
				    const auto f = mign_new.create_maj_10(operands);
					count[operands] = 0; 
					old_to_new_fo[node] = f; 
				    return f;
				}
				else 
				{
					//std::cout << "buffer! " << std::endl; 
  				  buffer_count[node] = 1; 
  				  count[operands]++; 
  				  const auto f = mign_new.create_maj_10(operands);
				 
  				   old_to_new_fo[node] = f; 
  				  // buffer[node] = f; 
  				  update(cp_node, father, mign); 
  				  std::vector<mign_function> operands_two;
  				  operands_two.push_back(mign_new.get_constant( false ));  
  				  operands_two.push_back(mign_new.get_constant( true ));  
  				  operands_two.push_back(f);  
  				  const auto g = mign_new.create_maj_10(operands_two);
  				  buffer[node] = g; 
  				  return g;
				}
				  
			  }
		}
		
	const auto f = mign_new.create_maj(operands); 
	return f;
}

mign_graph mign_fo_restricted_opt ( mign_graph& mign, unsigned int max_fo, bool depth_const)
{
	
	mign_graph mign_new; 
	
	std::string fo_res = "fo_rest_opt"; 
	
	if (!mign.name().empty())
	mign_new.set_name(mign.name()+fo_res); 
	else 
		mign_new.set_name("fo_rest_opt"); 
	
	mign_new.structural_hashing(1); 
	
	mign.compute_fanout(); 
	mign.compute_levels(); 
	mign.compute_parents(); 	
	
  /* create constant and PIs */
	
     auto old_to_new = init_visited_table_inv( mign, mign_new );
	 
	 std::map<mign_node,unsigned> cp_node; 
	 std::map<mign_function,unsigned> buffer_count; 
	 std::map<mign_node,mign_function> buffer;  
	 std::map<mign_node,mign_function> old_to_new_fo;
	 std::map<std::vector<mign_function>,unsigned> count;
	 
	 if (depth_const == 0)
	 {
	     int start_s=clock();
	     map_critpath_node (mign, cp_node); 
         int stop_s=clock();
	     std::cout << "Time to generate depth = " << (stop_s-start_s)/double(CLOCKS_PER_SEC)*1000 << std::endl; 
		 start_s=clock();
	     for ( const auto& output : mign.outputs())
	     {
	 		mign_new.create_po(mign_fo_restricted_rec_opt(mign, output.first.node, mign_new, old_to_new, old_to_new_fo, cp_node, max_fo, count, 0, buffer, buffer_count) ^ output.first.complemented, output.second); 
	 	 }
		 stop_s=clock();
		 std::cout << "Fo graphs = " << (stop_s-start_s)/double(CLOCKS_PER_SEC)*1000 << std::endl;
		 auto total_num_of_buffer = number_of_buffer(mign_new); 
		 std::cout << " Total number of buffers = " << total_num_of_buffer << std::endl;
		 std::cout << " New size before = " << mign_new.size() - mign_new.inputs().size() - 1 << std::endl;
		 mign_new = clean_mig(mign_new); 
		 std::cout << " New size after = " << mign_new.size() - mign_new.inputs().size() - 1 << std::endl;
		 total_num_of_buffer = number_of_buffer(mign_new); 
		 std::cout << " Total number of buffers = " << total_num_of_buffer << std::endl;
	 }
	
	  
	 else 
	 {
		 int start_s=clock();
	     for ( const auto& output : mign.outputs())
	     {
			 mign_new.create_po(mign_fo_restricted_rec_opt_depthconst(mign, output.first.node, mign_new, old_to_new, old_to_new_fo, cp_node, max_fo, count, 0, buffer, buffer_count) ^ output.first.complemented, output.second); 
	 	 }
		 int stop_s=clock();
		 std::cout << "Fo graphs = " << (stop_s-start_s)/double(CLOCKS_PER_SEC) << std::endl;
		 std::cout << " New size = " << mign_new.size() - mign_new.inputs().size() - 1 << std::endl; 
	 }
	
	 auto flag = 1; 
	 std::cout << " Final Check : " << std::endl; 
	 mign_new.compute_fanout(); 
	
	 auto node_fo = 0; 
	 
 	for ( auto& node : mign_new.topological_nodes())
 	{
		if(mign_new.is_input(node)) {continue; }
		if (mign_new.fanout_count(node) > 3)
		{
			flag = 0; 
			node_fo = node; 
			break; 
		}
	}
	
	if (flag == 1)
	{
		//std::cout << " C:  "; 
		std::cout << " FO ok " << std::endl; 
	}
	else 
		std::cout << "Error at least on FO error with node " << node_fo - mign_new.inputs().size() - 1  << std::endl;
	
	flag = 1; 
 	for ( auto& node : mign_new.topological_nodes())
 	{
		if(mign_new.is_input(node)) {continue; }
		auto c = mign_new.children(node); 
		for (auto & x : c)
		{
			if ((!mign_new.is_input(x.node)) && (x.node != 0) && (x.complemented == 1))
				flag = 0; 
		}
		if (flag == 0)
			break;
	}
	
	if (flag == 1)
	{
		//std::cout << " C:  "; 
		std::cout << " CE ok " << std::endl; 
	}
	else 
		std::cout << "Error CE"  << std::endl;
	
	//write_paths (mign_new);  
	return mign_new; 
}

}
// Local Variables:
// c-basic-offset: 2
// eval: (c-set-offset 'substatement-open 0)
// eval: (c-set-offset 'innamespace 0)
// End:
