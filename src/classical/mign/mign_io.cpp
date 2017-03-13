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

#include "mign_io.hpp"

#include <chrono>
#include <fstream>
#include <unordered_map>
#include <vector>
#include <string>

#include <boost/algorithm/string/join.hpp>
#include <boost/format.hpp>

#include <classical/functions/isop.hpp>
#include <classical/utils/truth_table_utils.hpp>

#include <classical/mign/mign_simulate.hpp>
#include <classical/mign/math_utils.hpp>
#include <core/utils/string_utils.hpp>
#include <core/cube.hpp>

using boost::format;

namespace cirkit
{

/******************************************************************************
 * Types                                                                      *
 ******************************************************************************/
	tt tt_const0_child( unsigned leaf_size)
	{
		auto c = int_pow2(leaf_size); 
		return boost::dynamic_bitset<>( c, 0u );
	}

	tt tt_const1_child(unsigned leaf_size)
	{
		  return ~tt_const0_child(leaf_size);
	}
/******************************************************************************
 * Private functions                                                          *
 ******************************************************************************/

std::string escape_name_mign( const std::string& name )
{
	
	if (name[0] == ' ')
	{
		return name.substr(1u); 
	}
  //if ( name.find( "[" ) )
  //{
  //  return "\\" + name;
  //}
  //else
  //{
    return name;
  //}
}

unsigned int factorial(unsigned int n) 
{
    if (n == 0)
       return 1;
    return n * factorial(n - 1);
}


std::vector<std::string> get_input_names( const mign_graph& mign )
{
  std::vector<std::string> inames;
  for ( const auto& input : mign.inputs() )
  {
	  //std::cout << input.second << std::endl; 
    inames.push_back( escape_name_mign( input.second ) );
  }
  return inames;
}

std::vector<std::string> get_output_names( const mign_graph& mign )
{
  std::vector<std::string> onames;
  for ( const auto& output : mign.outputs() )
  {
    onames.push_back( escape_name_mign( output.second ) );
  }
  return onames;
}

std::vector<std::string> get_wire_names( const mign_graph& mign )
{
  std::vector<std::string> wnames;
  for ( auto v : mign.vertices() )
  {
    if ( mign.is_input( v ) ) { continue; }
    wnames.push_back( boost::str( format( "w%d" ) % (v-mign.inputs().size()-1) ) );
  }
  return wnames;
}


std::string get_operand_name( const mign_graph& mign, const mign_function& f )
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
    name = escape_name_mign( mign.input_name( f.node ) );
  }
  else
  {
    name = boost::str( format( "w%d" ) % (f.node - mign.inputs().size()-1) );
  }

  if (( f.complemented ) && (f.node != 0))
  {
    name = "~" + name;
  }

  return name;
}

std::string get_operand_name_blif( const mign_graph& mign, const mign_function& f )
{
  std::string name;
  
  if ( mign.is_input( f.node ) )
  {
    name = escape_name_mign( mign.input_name( f.node ) );
  }
  else
  {
    name = boost::str( format( "w%d" ) % (f.node - mign.inputs().size()-1) );
  }

  return name;
}


void write_verilog_compact( const mign_graph& mign, std::ostream& os, const properties::ptr& settings = properties::ptr())
{
	auto time   = std::chrono::system_clock::now();
	auto time_c = std::chrono::system_clock::to_time_t( time );
	const auto header_prefix      = get( settings, "header_prefix",      std::string( "written by CirKit" ) );
	os << "// " << header_prefix << " " << std::ctime( &time_c ) << std::endl;
		
  const auto default_name = get( settings, "default_name", std::string( "top" ) );
 // const auto maj_module   = get( settings, "maj_module",   false );

  /* MAJ module */
  
    os << "module majority #(parameter n=91)" << std::endl; 
	os << "(input [n-1:0] A, output maj);" << std::endl; 
	os << "integer num_ones, bit;" << std::endl; 
	os << "always @(A)" << std::endl; 
	os << "begin" << std::endl; 
	os << "num_ones = 0;" << std::endl; 
	os << "for (bit = 0; bit<n;bit=bit+1) begin" << std::endl; 
	os << "if (A[bit] == 1'b1)" << std::endl; 
	os << "num_ones=num_ones+1;" << std::endl; 
	os << "else if (A[bit]==1'b0)" << std::endl; 
	os << "num_ones=num_ones;" << std::endl;
	os << "end" << std::endl;
	os << "end" << std::endl;
	os << "assign maj=(num_ones>(n-num_ones));" << std::endl;
	os << "endmodule" << std::endl << std::endl << std::endl << std::endl;
	

  /* top module */
	auto name = mign.name().empty() ? default_name : mign.name();
	//auto name = default_name; 

  auto inames = get_input_names( mign );
  auto onames = get_output_names( mign );
  auto wnames = get_wire_names( mign );

 // wnames.insert( wnames.begin(), "zero" );

  os << format( "module %s ( %s, %s );" ) % name % boost::join( inames, ", " ) % boost::join( onames, ", " ) << std::endl;

  os << format( "input %s;" ) % boost::join( inames, ", " ) << std::endl
     << format( "output %s;" ) % boost::join( onames, ", " ) << std::endl
     << format( "wire %s;" ) % boost::join( wnames, ", " ) << std::endl;
  


  auto inst_num = 0u; 
 
  for ( auto v : mign.topological_nodes() )
  {

    if ( mign.is_input( v ) ) { continue; }

    const auto c = mign.children( v );
	
	std::vector<std::string> c_names; 
	
	for (auto x :c)
	{
		c_names.push_back(get_operand_name( mign, x )); 
	}
	
	unsigned int size_c = c.size(); 

	os << format( "majority #(%d) inst%d ({%s}, w%d);") %size_c %inst_num %boost::join( c_names, ", " ) % (v - mign.inputs().size() -1) << std::endl; 
    
	++inst_num; 
  }

  for ( const auto& output : mign.outputs() )
  {
    os << format( "assign %s = %s;" ) % escape_name_mign( output.second ) % get_operand_name( mign, output.first ) << std::endl;
	//std::cout << " siamo al nodo" << output.second << " e di output vediamo " << output.first.node<<std::endl; 
  }

  os << "endmodule" << std::endl; 
}

void write_blif( const mign_graph& mign, std::ostream& os, const properties::ptr& settings = properties::ptr())
{
	auto time   = std::chrono::system_clock::now();
	auto time_c = std::chrono::system_clock::to_time_t( time );
	const auto header_prefix      = get( settings, "header_prefix",      std::string( "written by CirKit" ) );
	os << "# " << header_prefix << " " << std::ctime( &time_c ) << std::endl;
		
  const auto default_name = get( settings, "default_name", std::string( "top" ) );
 // const auto maj_module   = get( settings, "maj_module",   false );

  /* top module */
	auto name = mign.name().empty() ? default_name : mign.name();
	//auto name = default_name; 

  auto inames = get_input_names( mign );
  auto onames = get_output_names( mign );

  os << format(".model %s") % name << std::endl; 
 
  os << format( ".inputs %s" ) % boost::join( inames, " " ) << std::endl;
   os  << format( ".outputs %s" ) % boost::join( onames, " " ) << std::endl;
 
  for ( auto v : mign.topological_nodes() )
  {

    if ( mign.is_input( v ) ) { continue; }

    const auto c = mign.children( v );
	
	std::vector<std::string> c_names; 
	
	for (auto x :c)
	{
		if (x.node == 0) {continue;}
		
		c_names.push_back(get_operand_name_blif( mign, x )); 
	}
	

	os << format( ".names %s w%d ") %boost::join( c_names, " " ) % (v - mign.inputs().size() - 1) << std::endl; 
    
    std::map<mign_node, tt> inputs;
    auto i = 0u;
	auto c_size = c.size(); 
	
    for ( auto child : c )
    {
		if (child.node == 0)
			c_size--; 
	}
	
    for ( auto child : c )
    {
		tt f; 
		if (child.node == 0)
		{
			
			if (child.complemented == 1)
			f = tt_const0(); 
			else 
				f = tt_const1(); 	
		}
		else 
	    {
			f = tt_nth_var( i++ );
			if (tt_num_vars(f) < c_size )
				tt_extend(f,c_size);
			else
				tt_shrink(f,c_size); }
		
		std::cout << f << std::endl; 
		inputs[child.node] = f; 
	
    }

    mign_tt_simulator tt_sim;
	
    mign_partial_node_assignment_simulator<tt> sim( tt_sim, inputs, tt_const0() );
  
    auto func = simulate_mign_function( mign, v, sim );
	
	tt_shrink(func,c_size); 

	std::cout << " func = " << func << std::endl; 
	
	std::vector<int> cover; 
    const auto num_vars = tt_num_vars( func );
    tt_isop( func, func, cover );
    const auto sop = cover_to_cubes( cover, num_vars );
	common_pla_print( sop );
	
  	for ( const auto& c : sop )
  	  {
  	   os << format ("%s %d") %(c.to_string()) %1 << std::endl; 
	  }

  }

  auto flag = 1; 
  for ( const auto& output : mign.outputs() )
  {
	  if (output.first.node == 0) 
	  {
		  std::string zero= "zero"; 
		  if ( flag == 1)
			  {
				  flag =2;   
			  os << format( ".names %s" ) %zero << std::endl;   
		      }
		      os << format( ".names %s %s" ) %zero % escape_name_mign( output.second )  << std::endl;
	  }
    else 
		os << format( ".names %s %s" ) % get_operand_name_blif( mign, output.first ) % escape_name_mign( output.second )  << std::endl;
	
	if (output.first.complemented == 1)
		os << format( "%d %d" ) %1 %0 << std::endl; 
	else 
		os << format( "%d %d" ) %1 %1 << std::endl;
  }

  os << ".end" << std::endl; 
}


/******************************************************************************
 * Public functions                                                           *
 ******************************************************************************/
mign_function find_function( const std::unordered_map<std::string, mign_function>& name_to_function, const std::string& name )
{
  if ( name[0] == '~' )
  {
    return !name_to_function.at( name.substr( 1u ) );
  }
  else
  {
    return name_to_function.at( name );
  }
}

mign_graph read_verilog_compact( const std::string& filename )
{
  std::unordered_map<std::string, mign_function> name_to_function;
  std::vector<std::string> output_names;
  
  //std::vector<mign_function> operands; 

  mign_graph mign;

  name_to_function.insert( {"1'b0", mign.get_constant( false )} );
  name_to_function.insert( {"1'b1", mign.get_constant( true )} );
  

  line_parser( filename, {
      {std::regex( "^input (.*);" ), [&mign, &name_to_function]( const std::smatch& m ) {
		  
		  const auto str = std::string(m[1]);
		  std::smatch match;
		  
		  
		  if ( std::regex_search( str, match, std::regex( "^\\[(.*)\\] (.*)$" ) ) )
		  {
			  auto sm0 = std::string(match[1]); 
			  auto sm1 = std::string(match[2]); 
			  
			  auto m = sm0.size() - 3;
			  auto str2 = str.substr (1,m); 
			  
			  unsigned long ul = std::stoul (str2,nullptr,0);
			  
			  for ( auto i = 0; i <= ul; ++i)
			  {
				  auto name = boost::str( boost::format( sm1+"[%d]" ) % i);
				 // std::cout << " Nome aggiunto " << name << std::endl; 
				  name_to_function.insert( {name, mign.create_pi( name )} );
			  }
		  }

		  else {
          foreach_string( m[1], ", ", [&mign, &name_to_function]( const std::string& name ) {
              name_to_function.insert( {name, mign.create_pi( name )} );
			//  std::cout << "input " << name << std::endl; 
            } 
		
		);
          }} },
		
      {std::regex( "^output (.*);" ), [&output_names]( const std::smatch& m ) {
		  
		  const auto str = std::string(m[1]);
		  std::smatch match;
		  
		 // std::cout << "output " << str << std::endl; 
		  
		
		if ( std::regex_search( str, match, std::regex( "^\\[(.*)\\] (.*)$" ) ) )
		  {
			  auto sm0 = std::string(match[1]); 
			  auto sm1 = std::string(match[2]); 
			  
			  auto ma = sm0.size() - 3;
			  auto str2 = str.substr (1,ma); 
			  
			  unsigned long ul = std::stoul (str2,nullptr,0);
			  
			  for ( auto i = 0; i <= ul; ++i)
			  {
				  auto name = boost::str( boost::format( sm1+"[%d]" ) % i);
				  output_names.push_back(name); 
			  }
		  }

		  else {
		  
          split_string( output_names, m[1], ", " );
          }}},
		{std::regex ( "majority (.*) \\(\\{(.*)\\}, (.*)\\);"), [&mign, &name_to_function, &output_names] ( const std::smatch& m) {
			const auto name_wire = std::string (m[3]); 
			//const auto value = unsigned(m[1]); 
			//std::cout << "wire " << name_wire << std::endl; 
			
			std::vector<mign_function> operands; 
			
            foreach_string( m[2], ", ", [&mign, &name_to_function, &operands, &name_wire]( const std::string name) {
				
				
				operands.push_back(find_function (name_to_function, name)); 
				  
              } );
			  name_to_function.insert( {name_wire, mign.create_maj( operands )} );
			
		}},
		
		
      {std::regex( "^assign (.*) = (.*);" ), [&mign, &name_to_function, &output_names]( const std::smatch& m ) {
          const auto name = std::string( m[1] );
		  
		  //std::cout << "output " << name << std::endl; 
          const auto expr = std::string( m[2] );
          std::smatch match;
		  
		  if ((name != "maj") && (name != "zero")) {
            assert( boost::find( output_names, name ) != output_names.end() );
            name_to_function.insert( {name, find_function( name_to_function, expr )} );
		}
          }
      }}, false );

  for ( const auto& name : output_names )
  {
    mign.create_po( name_to_function[name], name );
  }

  return mign;
}



void write_verilog_compact( const mign_graph& mign, const std::string& filename , const properties::ptr& settings)
{
	//std::ostream& os = std::cout; 
  std::ofstream os( filename.c_str(), std::ofstream::out );
  write_verilog_compact( mign, os, settings );
}

void write_blif ( const mign_graph& mign,const std::string& filename , const properties::ptr& settings)
{
  std::ofstream os( filename.c_str(), std::ofstream::out );
  write_blif( mign, os, settings );
}

}
// Local Variables:
// c-basic-offset: 2
// eval: (c-set-offset 'substatement-open 0)
// eval: (c-set-offset 'innamespace 0)
// End:
