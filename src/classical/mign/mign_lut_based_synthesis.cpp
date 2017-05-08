/* CirKit: A circuit toolkit
 * Copyright (C) 2009-2015  University of Bremen
 * Copyright (C) 2015-2016  EPFL
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

#include "mign_lut_based_synthesis.hpp"

#include <vector>

#include <core/utils/conversion_utils.hpp>
#include <core/utils/graph_utils.hpp>
#include <core/utils/range_utils.hpp>
#include <core/utils/temporary_filename.hpp>
#include <core/utils/terminal.hpp>

#include <classical/lut/lut_graph.hpp>
#include <classical/io/read_blif.hpp>
#include <classical/functions/isop.hpp>
#include <classical/utils/truth_table_utils.hpp>


#include <core/cube.hpp>

#define timer timer_class
#include <boost/progress.hpp>
#undef timer

namespace cirkit
{

/******************************************************************************
 * Manager                                                                    *
 ******************************************************************************/

class mign_lut_based_synthesis_manager
{
public:
  mign_lut_based_synthesis_manager(const lut_graph_t& lut, const properties::ptr& settings , mign_graph& mign)
    : lut( lut ),
      verbose( get( settings, "verbose", verbose ) ), 
	  mign (mign)
	  {}
  

  mign_graph run()
  {
      const auto type = get( boost::vertex_lut_type, lut );
	  const auto names = get(boost::vertex_name, lut); 
	  std::vector<lut_vertex_t> topsort( boost::num_vertices( lut ) );
	  const auto spec = get( boost::vertex_lut, lut );
	  boost::topological_sort( lut, topsort.begin() );
	  
	  
	  for ( const auto& l : topsort)
	  {
		  auto flag = 1; 
		  mign_function add; 
		  if (type[l] == lut_type_t::pi ) 
			  add = mign.create_pi (names[l]); 
		  else if (type[l] == lut_type_t::po ) 
		  {
		  	const auto fanin = boost::out_degree( l, lut );
			assert (fanin == 1) ; 
			for ( auto n : boost::make_iterator_range( adjacent_vertices( l, lut ) ) )
			{
				const auto it = computed_mign.find(n);
                mign.create_po (it->second, names[l]); 
			}
			flag = 2; 
			
		  }
		  else if (type[l] == lut_type_t::gnd ) 
			 add = mign.get_constant(false);
		  else if (type[l] == lut_type_t::vdd ) 
			  add = mign.get_constant(true);
		  else if ( boost::out_degree( l, lut ) <= 1u )
		           {
					   std::vector<mign_function> operands; 
					   operands.push_back(mign.get_constant(false)); 
					   operands.push_back(mign.get_constant(true)); 
					   std::vector<lut_edge_t> children; 
					   
		   		    for ( const auto& e : boost::make_iterator_range( boost::out_edges( l, lut ) ) )
		   				children.push_back(e); 
					
		              if ( spec[l] == "2" )
		              {
	   					  auto target = boost::target( children[0], lut ); 
	   					  const auto it = computed_mign.find(target);
	                        operands.push_back(mign_function(it->second.node,0));  
		              }
		              else if ( spec[l] == "1" )
		              {
	   					  auto target = boost::target( children[0], lut ); 
	   					  const auto it = computed_mign.find(target);
	                        operands.push_back(mign_function(it->second.node,1));  
					}
    				 add = mign.create_maj(operands);  
   				  
				}
		  else  // internal node 
		  {
			  std::cout << " l " << l << " e un internalnode " << std::endl; 
			  auto hex = spec[l]; 
			 // std::cout << " hx = " << hex << std::endl; 
			  std::vector<int> cover; 
			  
			  tt func(convert_hex2bin(hex)); 
			 // std::cout << " convert " << convert_hex2bin(hex) << std::endl; 
			  const auto num_vars = tt_num_vars( func );
			  tt_isop( func, func, cover );
			  const auto sop = cover_to_cubes( cover, num_vars );
			  
	  	      common_pla_print( sop );
			  
			  std::vector<lut_edge_t> children; 
			  std::vector<mign_function> operands_or; 
			  
  		    for ( const auto& e : boost::make_iterator_range( boost::out_edges( l, lut ) ) )
  				children.push_back(e); 
			
			//std::cout << " children size = " << children.size() << std::endl;  
			auto fanin = boost::out_degree( l, lut );
			//std::cout << " fanin = " << fanin << std::endl; 
			std::vector<std::string> str;
			
		  	for ( const auto& c : sop )
		  	  {
		  	    str.push_back(c.to_string()); 
			}
			
			//std::cout << " string.size() = " << str.size() << std::endl; 
			
			for (auto& s : str)
			{
			      std::vector<mign_function> operand_and; 	
			  
		
		  for ( auto x = 0u; x < s.size(); ++x)
			  {
				  //std::cout << " s[x]" << s[x] << std::endl; 
				  if (s[x] == '-')
					  continue; 
				  else if (s[x] == '0')
				  {
					  //std::cout << " target" << std::endl; 
					  auto target = boost::target( children[x], lut ); 
					  const auto it = computed_mign.find(target);
                      operand_and.push_back(mign_function(it->second.node,1)); 
				  }
				  else if (s[x] == '1') // chiedere se e giusto?? 
				  {
					  auto target = boost::target( children[x], lut ); 
					  const auto it = computed_mign.find(target);
                      operand_and.push_back(mign_function(it->second.node,0)); 
				  }
			  }
		  
			  
			  operands_or.push_back(mign.create_and(operand_and)); 
		  }
		  if (operands_or.size() == 1)  // chiedere se e giusto questo ?? 
			 {
				 operands_or.push_back(mign.get_constant(false)); 
				 operands_or.push_back(mign.get_constant(true)); 
			 	 add = mign.create_maj(operands_or); 	 
			 } 
		 else
		  add = mign.create_or(operands_or); 
	  }
	  if (flag != 2)
	  computed_mign.insert({l,add}); 
  }
	  	
	  
    return mign;
  }

private:
	
  mign_graph                                         mign;

  const lut_graph_t&                                 lut;

  std::unordered_map<lut_vertex_t, mign_function>    computed_mign;
  
  bool                                               verbose = false;
  
};

/******************************************************************************
 * Private functions                                                          *
 ******************************************************************************/

/******************************************************************************
 * Public functions                                                           *
 ******************************************************************************/

mign_graph mign_lut_based_synthesis( const std::string& filename, const properties::ptr& settings, const properties::ptr& statistics )
{
  /* timing */
  //properties_timer t( statistics );
  
  /* store_cubes = false if circuiti is small */
  bool store_cubes = false; 
  mign_graph mign; 
 
  const auto lut = read_blif (filename, store_cubes); // lut_graph_t& lut
  
  // dalla cover so quali sono p1 + p2 + p3 + p4... con p1 = l1 * l2 * l3 * l4... e li creo con MAJ-n

  mign_lut_based_synthesis_manager mgr( lut, settings, mign );
  return mgr.run();
}

}

// Local Variables:
// c-basic-offset: 2
// eval: (c-set-offset 'substatement-open 0)
// eval: (c-set-offset 'innamespace 0)
// End:
