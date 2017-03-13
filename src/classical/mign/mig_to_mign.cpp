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

#include "mig_to_mign.hpp"

#include <fstream>
#include <iostream>

#include <boost/optional/optional_io.hpp>
#include <boost/format.hpp> 
#include <classical/mig/mig_utils.hpp>

using boost::format;

namespace cirkit 
{
	
mign_graph mig_to_mign_not_strash (const mig_graph& mig)
{
	std::unordered_map<std::string, mign_function> name_to_function; 
	std::map<detail::mig_traits_t::vertex_descriptor, std::string> wires; 
	
	const auto& info = mig_info(mig); 
	mign_graph mign (3);
	
	mign.structural_hashing( true); 
	
    name_to_function.insert( {"1'b0", mign.get_constant( false )} );
    name_to_function.insert( {"1'b1", mign.get_constant( true )} );
	
	
	for (const auto& input : info.inputs)
	{
		auto name = info.node_names.at(input); 
		name_to_function.insert({name, mign.create_pi(name)}); 
	}
	auto count = 0u;
	for (const auto& node : boost::make_iterator_range( boost::vertices( mig ) ) )
	{
        if ( !boost::out_degree( node, mig ) ) { continue; }
		auto name = boost::str( format( "w%d" ) % count++);
		wires.insert({node,name}); 
	}

	
    for ( const auto& node : boost::make_iterator_range( boost::vertices( mig ) ) )
    {
      if ( !boost::out_degree( node, mig ) ) { continue; }

	  //std::cout << " stiamo scrivendo il nodo " << node << std::endl; 
      const auto children = get_children( mig, node );
	  std::vector<mign_function> operands;
	  mign_function function;  
	  
	  for ( auto& x : children)
	  {
		 
		  if (x.node == 0)
		  {
			  if (x.complemented == 1)
			  {
			  	function = name_to_function.at("1'b1");
			  }
			  else 
				function = name_to_function.at("1'b0");
		  }
		  else 
		  {
			  auto name = boost::str( format( "w%d" ) % (node-mign.inputs().size() -1));
			  if (!boost::out_degree( x.node, mig ))
			  name = info.node_names.at(x.node); 
			  else 
			  name = wires.at(x.node); 
		 	  //std::cout << " come ti chiami : " << name << std::endl; 
			  function = name_to_function.at(name); 
		  	  if (x.complemented == 1)
		  	  function = !function;  
		  }
		  operands.push_back(function);
			   
	  }
	  
	  auto name_node = wires.at(node); 
	  name_to_function.insert({name_node, mign.create_maj(operands)});
	  
  }
	  for ( const auto& output : info.outputs )
	  {
		  mign_function function;  
	  
		  auto name_out = output.second; 
		  std::string name; 
		  
		  if (!boost::out_degree( output.first.node, mig ))
		  {
			  if (output.first.node == 0)
			  {
				  if (output.first.complemented == 1)
				  {
				  	function = name_to_function.at("1'b1");
					mign.create_po(function,name_out );
				  }
				  else {
					function = name_to_function.at("1'b0"); 
				  mign.create_po(function,name_out );
			  	}
			  }
			  else {
		  name = info.node_names.at(output.first.node); 
		  if (output.first.complemented == 1)
		   mign.create_po(!name_to_function[name],name_out );
		  else 
		mign.create_po(name_to_function[name],name_out );
		  }
	  }
		  //name = info.node_names.at(output.first.node); 
		  else 
		  {
			 name = wires.at(output.first.node);
   		  if (output.first.complemented == 1)
   		   mign.create_po(!name_to_function[name],name_out );
   		  else 
   		mign.create_po(name_to_function[name],name_out );
	      }
	  }
  
	return mign; 
}
mign_graph mig_to_mign (const mig_graph& mig)
{

	std::unordered_map<std::string, mign_function> name_to_function; 
	std::map<detail::mig_traits_t::vertex_descriptor, std::string> wires; 
	
	std::string default_name = "default"; 
	
	const auto& info = mig_info(mig); 
	mign_graph mign (3);
	
	if (!info.model_name.empty())
	mign.set_name(info.model_name); 
	else 
		mign.set_name(default_name); 
		
	
    name_to_function.insert( {"1'b0", mign.get_constant( false )} );
    name_to_function.insert( {"1'b1", mign.get_constant( true )} );
	
	
	for (const auto& input : info.inputs)
	{
		auto name = info.node_names.at(input); 
		name_to_function.insert({name, mign.create_pi(name)}); 
	}
	auto count = 0u;
	for (const auto& node : boost::make_iterator_range( boost::vertices( mig ) ) )
	{
        if ( !boost::out_degree( node, mig ) ) { continue; }
		auto name = boost::str( format( "w%d" ) % count++);
		wires.insert({node,name}); 
	}

	
    for ( const auto& node : boost::make_iterator_range( boost::vertices( mig ) ) )
    {
      if ( !boost::out_degree( node, mig ) ) { continue; }

	  //std::cout << " stiamo scrivendo il nodo " << node << std::endl; 
      const auto children = get_children( mig, node );
	  std::vector<mign_function> operands;
	  mign_function function;  
	  
	  for ( auto& x : children)
	  {
		 
		  if (x.node == 0)
		  {
			  if (x.complemented == 1)
			  {
			  	function = name_to_function.at("1'b1");
			  }
			  else 
				function = name_to_function.at("1'b0");
		  }
		  else 
		  {
			  auto name = boost::str( format( "w%d" ) % (node-mign.inputs().size() -1));
			  if (!boost::out_degree( x.node, mig ))
			  name = info.node_names.at(x.node); 
			  else 
			  name = wires.at(x.node); 
		 	  //std::cout << " come ti chiami : " << name << std::endl; 
			  function = name_to_function.at(name); 
		  	  if (x.complemented == 1)
		  	  function = !function;  
		  }
		  operands.push_back(function);
			   
	  }
	  
	  auto name_node = wires.at(node); 
	  name_to_function.insert({name_node, mign.create_maj(operands)});
	  
  }
	  for ( const auto& output : info.outputs )
	  {
		  mign_function function;  
	  
		  auto name_out = output.second; 
		  std::string name; 
		  
		  if (!boost::out_degree( output.first.node, mig ))
		  {
			  if (output.first.node == 0)
			  {
				  if (output.first.complemented == 1)
				  {
				  	function = name_to_function.at("1'b1");
					mign.create_po(function,name_out );
				  }
				  else {
					function = name_to_function.at("1'b0"); 
				  mign.create_po(function,name_out );
			  	}
			  }
			  else {
		  name = info.node_names.at(output.first.node); 
		  if (output.first.complemented == 1)
		   mign.create_po(!name_to_function[name],name_out );
		  else 
		mign.create_po(name_to_function[name],name_out );
		  }
	  }
		  //name = info.node_names.at(output.first.node); 
		  else 
		  {
			 name = wires.at(output.first.node);
   		  if (output.first.complemented == 1)
   		   mign.create_po(!name_to_function[name],name_out );
   		  else 
   		mign.create_po(name_to_function[name],name_out );
	      }
	  }
  
	return mign; 
	
}

mig_graph mign_to_mig (const mign_graph& mign)
{
	mig_graph mig; 
    mig_initialize( mig, "from_mign" );

    //auto& info = mig_info( mig );
	std::unordered_map<std::string, mig_function>    node_to_function; 
	std::map<unsigned, std::string>                  name_to_function;
	std::vector<mig_function>                        functions;
	
	for (auto & input : mign.inputs())
	{
		name_to_function.insert(std::make_pair(input.first, input.second)); 
		node_to_function.insert(std::make_pair(input.second, mig_create_pi(mig, input.second))); 
	}
	
    for ( auto v : mign.topological_nodes() )
    {
      if ( mign.is_input( v ) ) { continue; }
      const auto c = mign.children( v );
	  if (c.size() != 3)
	  {
		  std::cout << " Error. Cannot be changed into mig3" << std::endl; 
	  }
	  std::vector<mig_function> children; 
	  
	  for (auto& x : c)
	  {
		  if (x.node == 0)
		  {
			  if (x.complemented == 1)
			  {
				children.push_back(mig_get_constant( mig, true )); 
			  }
			  else 
				children.push_back(mig_get_constant( mig, false )); 
		  }
		  else 
		  {
			   auto name = name_to_function[x.node];
			   if (x.complemented == 1)
			   children.push_back(!node_to_function.at(name)); 
			   else  
			   children.push_back(node_to_function.at(name)); 
		  }
	  }
	 auto wire =  boost::str( format( "w%d" ) % (v -mign.inputs().size() -1));
	name_to_function.insert(std::make_pair(v, wire)); 
	node_to_function.insert(std::make_pair(wire, mig_create_maj(mig,children[0],children[1],children[2]))); 
	
    }
	
	for (auto& output : mign.outputs() )

		{
			auto name_out = output.second;
			if (output.first.node == 0)
  			 {
  			 	 if (output.first.complemented == 1)
  			  {
  				mig_create_po (mig,mig_get_constant( mig, true ),name_out );
  			  }
  			  else 
  				mig_create_po (mig,mig_get_constant( mig, false ),name_out );
		    }
			else 
			{
				auto name = name_to_function.at(output.first.node);
				if (output.first.complemented == 1)
				mig_create_po (mig,!node_to_function[name],name_out ); 
				else 
				mig_create_po (mig,node_to_function[name],name_out ); 
			}		 
			
		}
		
		return mig; 
	}
}


// Local Variables:
// c-basic-offset: 2
// eval: (c-set-offset 'substatement-open 0)
// eval: (c-set-offset 'innamespace 0)
// End:
