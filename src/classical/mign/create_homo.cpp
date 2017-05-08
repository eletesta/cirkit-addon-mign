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

#include "create_homo.hpp"

#include <fstream>
#include <iostream>

#include <boost/optional/optional_io.hpp>
#include <boost/format.hpp>
#include <boost/graph/graphviz.hpp>

using boost::format;

namespace cirkit 
{

	// Public function 
mign_graph create_homo (const mign_graph& mign, unsigned Nin) 
{
    std::unordered_map<std::string, mign_function> name_to_function;
    std::vector<std::string> output_names;
	
  
    mign_graph mign_new(Nin); // homogeneous, each node with Nin inputs 
	
    name_to_function.insert( {"1'b0", mign_new.get_constant( false )} );
    name_to_function.insert( {"1'b1", mign_new.get_constant( true )} );

	
	 for ( const auto& input : mign.inputs() )
	 {
		 auto str = input.second; 
		 
		 name_to_function.insert ({str, mign_new.create_pi(str)});
		 //std::cout << " Dopo dei PI ;)" << std::endl; 
	 }
  
     for ( auto v : mign.topological_nodes() )
     {
       if ( mign.is_input( v ) ) { continue; }

       const auto c = mign.children( v );
	
	   // thrree cases to build homo, but they are consdiered into create_maj and create_and and or  : 1. alrady correct 2. smaller, then add some 1 and 0 to reach the right number 3. larger --> pensarci bene!! ;) 
	   // std::cout << " Prima di create_maj" << std::endl; 
	   
		name_to_function.insert( {boost::str( format( "w%d" ) % v ), mign_new.create_maj( c )} );
		//std::cout << " Dopo di create_maj" << std::endl; 
	}
	
	//std::cout << " Prima degli output" << std::endl; 
	for (const auto& output : mign.outputs())
	{
		 mign_new.create_po( output.first, output.second);
	}
	
	return mign_new; 
}

}