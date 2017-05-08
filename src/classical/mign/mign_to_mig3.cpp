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

#include "mign_to_mig3.hpp"

#include <fstream>
#include <iostream>

#include <boost/optional/optional_io.hpp>
#include <boost/format.hpp> 
#include <classical/mign/mign_io.hpp>
#include <classical/mign/arith_mign.hpp>

using boost::format;

namespace cirkit 
{
	

std::string get_name( const mign_graph& mign, const mign_function& f )
	{
	  std::string name;

	  if ( f.node == 0u )
	  {
		  if (f.complemented)
			  name = "1'b1"; 
		  else
	    name = "1'b0";
	  }
  
	  else if ( mign.is_input( f.node ) )
	  {
	    name = ( mign.input_name( f.node ) );
	  }
	  else
	  {
	    name = boost::str( format( "w%d" ) % (f.node - mign.inputs().size()-1) );
	  }

	  return name;
	}


mign_graph mign_to_mig3 (const mign_graph& mign_old)
{

	std::unordered_map<std::string, mign_function> name_to_function; 
 
	mign_graph mign_new (3);
	
	auto count = 0u; 
	
    name_to_function.insert( {"1'b0", mign_new.get_constant( false )} );
    name_to_function.insert( {"1'b1", mign_new.get_constant( true )} );
	
	//std::cout << "debugging 1" << std::endl; 
	for (const auto& input : mign_old.inputs())
	{
		auto name = input.second;  
		name_to_function.insert({name, mign_new.create_pi(name)}); 
	}
	
	std::vector<std::string> wnames;
	for (const auto& v : mign_old.vertices())
	{
		if ( mign_old.is_input( v ) ) { continue; }
		wnames.push_back( boost::str( format( "w%d" ) % (v-mign_old.inputs().size()-1) ) );
	}

	std::vector<std::string> onames; 
    for ( const auto& output : mign_old.outputs() )
    {
      onames.push_back(output.second);
    }
	//std::cout << "debugging 2" << std::endl; 
	
    for ( const auto& node : mign_old.topological_nodes() )
    {
      if ( mign_old.is_input( node ) ) { continue; }

	  const auto c = mign_old.children( node );
      if (c.size() == 3)
	  {
		   //std::cout << "mig da 3" << std::endl; 
		  std::vector<mign_function> operands;
		  for ( auto& x : c)
		  {
			  auto name = get_name(mign_old,x); 
			  if ((x.complemented == 1) && (x.node != 0))
			  {
			  	operands.push_back(!find_function (name_to_function, name));
			  }
			  else
			  operands.push_back(find_function (name_to_function, name));
		  }
		  
		  auto name_wire = boost::str( format( "w%d" ) % (node-mign_new.inputs().size()-1) ) ;
		  name_to_function.insert({name_wire,mign_new.create_maj(operands)}); 
	  }
	  
	  else if (c.size() == 5)
	  {
		  //::cout << "mig da 5" << std::endl;
		  std::vector<mign_function> operands; 
		  
		  ++count; 
		  
		  auto name_wire_one = boost::str( format( "w%d_new1" ) % (node-mign_new.inputs().size()-1) ) ; 
		  auto name_one = get_name(mign_old,c[0]);
		 // std::cout << name_one << std::endl; 
		  auto name_two = get_name(mign_old,c[1]);
		 // std::cout << name_two << std::endl; 
		  auto name_three = get_name(mign_old,c[2]);
		 // std::cout << name_three << std::endl; 
		  if ((c[0].complemented == 1) && (c[0].node != 0))
		  {
		  	operands.push_back(!find_mign(name_to_function,name_one));
		  }
		  else 
		  operands.push_back(find_mign(name_to_function,name_one));
		  
		  
		  if ((c[1].complemented == 1) && (c[1].node != 0))
		  {
		  	operands.push_back(!find_mign(name_to_function,name_two));
		  }
		  else 
		  operands.push_back(find_mign(name_to_function,name_two));
		  
		  
		  if ((c[2].complemented == 1) && (c[2].node != 0))
		  {
		  	operands.push_back(!find_mign(name_to_function,name_three));
		  }
		  else 
		  operands.push_back(find_mign(name_to_function,name_three));
		    
		  name_to_function.insert({name_wire_one,mign_new.create_maj(operands)});
		  
		  operands.erase(operands.begin(),operands.end());
		  //std::cout << "nodo 2"<< std::endl;
		  ++count; 
		  auto name_wire_two = boost::str( format( "w%d_new2" ) % (node-mign_new.inputs().size()-1) ) ;
		 // std::cout << " aggiungo il wire " << name_wire_two << std::endl; 
		  name_one = get_name(mign_old,c[1]);
		  name_two = get_name(mign_old,c[2]);
		  name_three = get_name(mign_old,c[3]);
		  if ((c[1].complemented == 1) && (c[1].node != 0))
		  {
		  	operands.push_back(!find_mign(name_to_function,name_one));
		  }
		  else 
		  operands.push_back(find_mign(name_to_function,name_one));
		  
		  if ((c[2].complemented == 1) && (c[2].node != 0))
		  {
		  	operands.push_back(!find_mign(name_to_function,name_two));
		  }
		  else 
		  operands.push_back(find_mign(name_to_function,name_two));
		  
		  
		  if ((c[3].complemented == 1) && (c[3].node != 0))
		  {
		  	operands.push_back(!find_mign(name_to_function,name_three));
		  }
		  else 
		  operands.push_back(find_mign(name_to_function,name_three));
		  
		  name_to_function.insert({name_wire_two,mign_new.create_maj(operands)});
		  
		  operands.erase(operands.begin(),operands.end());
		  
		  //second new node 
		 // std::cout << "nodo 3"<< std::endl;
		  ++count; 
		  auto name_wire_three = boost::str( format( "w%d_new3" ) % (node-mign_new.inputs().size()-1) ) ;
		 // std::cout << " aggiungo il wire " << name_wire_three << std::endl; 
		  name_one = get_name(mign_old,c[0]);
		  if ((c[0].complemented == 1) && (c[0].node != 0))
		  {
		  	operands.push_back(!find_mign(name_to_function,name_one));
		  }
		  else 
		  operands.push_back(find_mign(name_to_function,name_one));
		  
		  if ((c[3].complemented == 1) && (c[3].node != 0))
		  {
		  	operands.push_back(!find_mign(name_to_function,name_three));
		  }
		  else 
		  operands.push_back(find_mign(name_to_function,name_three));
		  
		  operands.push_back(find_mign(name_to_function,name_wire_one));
		  
		  name_to_function.insert({name_wire_three,mign_new.create_maj(operands)});
		  
		  operands.erase(operands.begin(),operands.end());
		//  std::cout << "nodo 4"<< std::endl;
		  auto name_wire = boost::str( format( "w%d" ) % (node-mign_new.inputs().size()-1) ) ;
		  //std::cout << name_wire << std::endl; 
		  name_one = get_name(mign_old,c[4]);
		  if ((c[4].complemented == 1) && (c[4].node != 0))
		  {
		  	operands.push_back(!find_mign(name_to_function,name_one));
		  }
		  else 
		  operands.push_back(find_mign(name_to_function,name_one));
		  
 		  operands.push_back(find_mign(name_to_function,name_wire_two));
 		  operands.push_back(find_mign(name_to_function,name_wire_three));
		  
		  name_to_function.insert({name_wire,mign_new.create_maj(operands)});
		  
		  
	  }
	  else if (c.size() == 7)
	  {
		   // std::cout << "mig da 7" << std::endl; 
		  std::vector<mign_function> operands; 
		  
		  // first new node 
		  ++count; 
		  auto name_wire_one = boost::str( format( "w%d_new1" ) % (node-mign_new.inputs().size()-1) ) ;
		  auto name_one = get_name(mign_old,c[0]);
		  auto name_two = get_name(mign_old,c[1]);
		  auto name_three = get_name(mign_old,c[2]);
		  if ((c[0].complemented == 1) && (c[0].node != 0))
		  {
		  	operands.push_back(!find_mign(name_to_function,name_one));
		  }
		  else 
		  operands.push_back(find_mign(name_to_function,name_one));
		  if ((c[1].complemented == 1) && (c[1].node != 0))
		  {
		  	operands.push_back(!find_mign(name_to_function,name_two));
		  }
		  else 
		  operands.push_back(find_mign(name_to_function,name_two));
		  if ((c[2].complemented == 1) && (c[2].node != 0))
		  {
		  	operands.push_back(!find_mign(name_to_function,name_three));
		  }
		  else 
		  operands.push_back(find_mign(name_to_function,name_three));
		  name_to_function.insert({name_wire_one,mign_new.create_maj(operands)});
		  
		  operands.erase(operands.begin(),operands.end());
		  
		  //second new node 
		  ++count; 
		  auto name_wire_two = boost::str( format( "w%d_new2" ) % (node-mign_new.inputs().size()-1) ) ;
		  if ((c[0].complemented == 1) && (c[0].node != 0))
		  {
		  	operands.push_back(!find_mign(name_to_function,name_one));
		  }
		  else 
		  operands.push_back(find_mign(name_to_function,name_one));
		  if ((c[1].complemented == 1) && (c[1].node != 0))
		  {
		  	operands.push_back(!find_mign(name_to_function,name_two));
		  }
		  else 
		  operands.push_back(find_mign(name_to_function,name_two));
		  
		  operands.push_back(!find_mign(name_to_function,name_wire_one));
		  
		  name_to_function.insert({name_wire_two,mign_new.create_maj(operands)});
		  
		  operands.erase(operands.begin(),operands.end());
		   
		   //third new node 
		  ++count; 
 		  auto name_wire_three = boost::str( format( "w%d_new3" ) % (node-mign_new.inputs().size()-1) );
		  name_one = get_name(mign_old,c[3]);
		  name_two = get_name(mign_old,c[4]);
		  name_three = get_name(mign_old,c[5]);
		  if ((c[3].complemented == 1) && (c[3].node != 0))
		  {
		  	operands.push_back(!find_mign(name_to_function,name_one));
		  }
		  else 
		  operands.push_back(find_mign(name_to_function,name_one));
		  if ((c[4].complemented == 1) && (c[4].node != 0))
		  {
		  	operands.push_back(!find_mign(name_to_function,name_two));
		  }
		  else 
		  operands.push_back(find_mign(name_to_function,name_two));
		  if ((c[5].complemented == 1) && (c[5].node != 0))
		  {
		  	operands.push_back(!find_mign(name_to_function,name_three));
		  }
		  else 
		  operands.push_back(find_mign(name_to_function,name_three));
		  
 		  name_to_function.insert({name_wire_three,mign_new.create_maj(operands)});
		  
		  operands.erase(operands.begin(),operands.end());
		  //forth node 
		  ++count; 
 		  auto name_wire_four = boost::str( format( "w%d_new4" ) % (node-mign_new.inputs().size()-1) );
		  if ((c[3].complemented == 1) && (c[3].node != 0))
		  {
		  	operands.push_back(!find_mign(name_to_function,name_one));
		  }
		  else 
		  operands.push_back(find_mign(name_to_function,name_one));
		  if ((c[4].complemented == 1) && (c[4].node != 0))
		  {
		  	operands.push_back(!find_mign(name_to_function,name_two));
		  }
		  else 
		  operands.push_back(find_mign(name_to_function,name_two));
 		  operands.push_back(find_mign(name_to_function,name_wire_one));
		  
 		  name_to_function.insert({name_wire_four,mign_new.create_maj(operands)});
		  
		  operands.erase(operands.begin(),operands.end());
		  
		  auto name_wire_five = boost::str( format( "w%d_new5" ) % (node-mign_new.inputs().size()-1) ) ;
		  ++count; 
		  name_one = get_name(mign_old,c[0]);
		  if ((c[0].complemented == 1) && (c[0].node != 0))
		  {
		  	operands.push_back(!find_mign(name_to_function,name_one));
		  }
		  else 
		  operands.push_back(find_mign(name_to_function,name_one));
		  
 		  operands.push_back(find_mign(name_to_function,name_wire_two));
 		  operands.push_back(find_mign(name_to_function,name_wire_three));
		  
		  name_to_function.insert({name_wire_five,mign_new.create_maj(operands)});
		  
		  operands.erase(operands.begin(),operands.end());
		  
		  auto name_wire_six = boost::str( format( "w%d_new6" ) % (node-mign_new.inputs().size()-1) ) ;
		  ++count; 
		  if ((c[5].complemented == 1) && (c[5].node != 0))
		  {
		  	operands.push_back(!find_mign(name_to_function,name_three));
		  }
		  else 
		  operands.push_back(find_mign(name_to_function,name_three));
 		  
 		  operands.push_back(find_mign(name_to_function,name_wire_one));
 		  operands.push_back(find_mign(name_to_function,name_wire_four));
		  
		  name_to_function.insert({name_wire_six,mign_new.create_maj(operands)});
		  
		  operands.erase(operands.begin(),operands.end());
		  
		  auto name_wire = boost::str( format( "w%d" ) % (node-mign_new.inputs().size()-1) ) ;
		  name_one = get_name(mign_old,c[6]);
		  if ((c[6].complemented == 1) && (c[6].node != 0))
		  {
		  	operands.push_back(!find_mign(name_to_function,name_one));
		  }
		  else 
		  operands.push_back(find_mign(name_to_function,name_one));
		  
 		  operands.push_back(find_mign(name_to_function,name_wire_five));
 		  operands.push_back(find_mign(name_to_function,name_wire_six));
		  
		  name_to_function.insert({name_wire,mign_new.create_maj(operands)});
	  }
      else if (c.size() == 9)
	  {
		std::cout << " To be done. " << std::endl; 
		// TODO 
	  }		  
      else 
	  {
		std::cout << " Error. The n for the maj is too big. It changes until 9 input maj into 3maj" << std::endl; 
	  }
  }
  
  //std::cout << "output.." << std::endl;
  for ( auto& output : mign_old.outputs())
  {
	  
	  auto name_out = output.second; 
	 // std::cout << output.second << std::endl; 
	  
	  if (output.first.node == 0)
	  {
		  if (output.first.complemented == 1)
		  {
			  auto name = "1'b1"; 
			  mign_new.create_po(find_mign(name_to_function,name),name_out);	  
		  }
		  else 
		  {
			  auto name = "1'b0"; 
			  mign_new.create_po(find_mign(name_to_function,name),name_out);	
		  }
	  }
	 else if (mign_old.is_input(output.first.node))
		{
			 auto name = mign_old.input_name(output.first.node); 
   		  //std::cout << name <<std::endl; 
   	  	  mign_new.create_po(find_mign(name_to_function,name)^output.first.complemented,name_out);
	     }
    else 
	  {
		  auto name = boost::str( format( "w%d" ) % (output.first.node-mign_new.inputs().size()-1) ) ;
		 // std::cout << name <<std::endl; 
	  	  mign_new.create_po(find_mign(name_to_function,name)^output.first.complemented,name_out);
	  }
		 
  }
  
	return mign_new; 
	
}
}

// Local Variables:
// c-basic-offset: 2
// eval: (c-set-offset 'substatement-open 0)
// eval: (c-set-offset 'innamespace 0)
// End:
