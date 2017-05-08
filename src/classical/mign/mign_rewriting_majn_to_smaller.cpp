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

#include "mign_rewriting_majn_to_smaller.hpp"

#include <functional>

#include <boost/assign/std/vector.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/format.hpp>
#include <boost/range/iterator_range.hpp>

#include <core/graph/depth.hpp>
#include <core/utils/timer.hpp>
#include <core/utils/range_utils.hpp>
#include <classical/mign/mign_utils.hpp>
#include <classical/mign/mign_rewrite.hpp>

using namespace boost::assign;

namespace cirkit
{

/******************************************************************************
 * Types                                                                      *
 ******************************************************************************/

using mign_function_vec_t = std::vector<mign_function>;

class mign_rewriting_majn_to_smaller_manager
{
public:
  mign_rewriting_majn_to_smaller_manager( const mign_graph& mign, bool verbose, unsigned n);
  
  void swap_current_mign( const std::string& method );
  
  void run_changen_to_smaller();
  void schensted();
  mign_function mign_change_to_smaller( const mign_node& node );
  mign_function mign_schensted( const mign_node& node );
  
private:

inline unsigned counting (std::vector<unsigned> v)
	{
		auto count = 0; 
		for (auto x = 0; x < v.size() ; x++)
		{
			if (v[x]== 1)
				count++; 
		}
		return count; 
	}
	
inline std::vector<std::vector<unsigned>> eval_combination (unsigned h)
{  
    std::vector<unsigned> a(h+1,0u); 
    std::vector<unsigned> m(h+1,2u); 
   
	std::vector<std::vector<unsigned>> combinations;
	
    mixed_radix( a, m, [&]( const std::vector<unsigned>& a ) {
       
		auto c = a; 
		c.erase (c.begin() + 0); 
 
  	  if (counting(c) == h-1)
  	  combinations.push_back(c); 
  
          return false;
        } );
 
   return combinations; 
}
   
   
  inline mign_function mign_equaln(const mign_function_vec_t children, std::pair<unsigned,unsigned> t)
  {
	  
	  std::cout << " equal n" << std::endl; 
    std::vector<std::vector<mign_function>> operands; 
	const auto e = (children.size()-1)/2; 
	operands.resize(e); 
	
	//std::cout << operands.size() << std::endl;  // numero delle operazioni a 
	const auto k = children.size()-2;  // 5 ingressi l'uno 
	
	//std::cout << k << std::endl; 
	auto provv = children;
	
	const auto d = mign_change_to_smaller( children[t.first].node)^ children[t.first].complemented; 
	for ( auto & o : operands)
	{ 
		o.push_back(d);  // all with the common operands ;) 
	}
	
	provv.erase(provv.begin() + t.first); 
	
	provv.erase(provv.begin() + t.second - 1); 
	
	const auto size_p = provv.size(); 
	
	auto h = k - e; 
	//std::cout << "h" << h << std::endl;
	if ( h == 1) // caso da 5 a 3
	{
		operands[0].push_back(mign_change_to_smaller( provv[0].node)^ provv[0].complemented); 
	}
	else 
	{
		// ci sono h numeri da due possibilita l'uno. ovvero o c'e o non c'e 
		
		for ( auto j = 0; j < h; j++)
		{
			const auto f = mign_change_to_smaller( provv[j].node)^ provv[j].complemented;  
			operands[0].push_back(f);
		}
		
		auto comb =  eval_combination (h); 
		
		auto g = 1; 
		for (auto & c : comb)
		{
			for ( auto x = 0; x < h; ++x)
			{
				if ( c[x] == 1)
					operands[g].push_back(mign_change_to_smaller( provv[x].node)^ provv[x].complemented); 
			}
			g++;
		}
	}
	
	//std::cout << provv.size() << std::endl; 
	for ( auto j = h; j < size_p ; j++)
	{
		for ( auto o = 1; o < operands.size(); o++)
		{
			operands[o].push_back(mign_change_to_smaller( provv[j].node)^ provv[j].complemented); 
		}
	}
	
	for ( auto o = 1; o < operands.size() ; o++)
	{
		auto d = mign_current.create_maj(operands[o]); 
		std::cout << " aggiunto nodo " << d.node <<std::endl; 
		operands[0].push_back(d); 
	}
	
	//std::cout << provv.size() << std::endl;
	std::cout << " operands[0] " << operands[0].size() << std::endl; 
	return mign_current.create_maj(operands[0]); 
    
  }
  
  inline mign_function mign_equaln_inside(const mign_function_vec_t children, std::pair<unsigned,unsigned> t)
  {
	 
    std::vector<std::vector<mign_function>> operands; 
	const auto e = (children.size()-1)/2; 
	operands.resize(e); 
	
	//std::cout << operands.size() << std::endl;  // numero delle operazioni a 
	const auto k = children.size()-2;  // 5 ingressi l'uno 
	
	//std::cout << k << std::endl; 
	auto provv = children;
	
	const auto d = children[t.first]; 
	for ( auto & o : operands)
	{ 
		o.push_back(d);  // all with the common operands ;) 
	}
	
	provv.erase(provv.begin() + 0); 
	
	provv.erase(provv.begin() + 0); 
	
	const auto size_p = provv.size(); 
	
	auto h = k - e; 
	
	if ( h == 1) // caso da 5 a 3
	{
		operands[0].push_back(provv[0]); 
	}
	else 
	{
		// ci sono h numeri da due possibilita l'uno. ovvero o c'e o non c'e 
		
		for ( auto j = 0; j < h; j++)
		{
			const auto f = provv[j]; 
			operands[0].push_back(f);
		}
		
		auto comb =  eval_combination (h); 
		
		for ( auto c : comb)
		{
			for ( auto x = 0; x < c.size(); x++)
			{
				std::cout << c[x]; 
			}
			std::cout << std::endl; 
		}
		auto g = 1; 
		for (auto & c : comb)
		{
			for ( auto x = 0; x < h; ++x)
			{
				if ( c[x] == 1)
					operands[g].push_back(provv[x]); 
			}
			g++;
		}
	}
	
	//std::cout << provv.size() << std::endl; 
	for ( auto j = h; j < size_p ; j++)
	{
		for ( auto o = 1; o < operands.size(); o++)
		{
			operands[o].push_back(provv[j]); 
		}
	}

	
	for ( auto o = 1; o < operands.size() ; o++)
	{
		auto d = mign_current.create_maj(operands[o]); 
		std::cout << " aggiunto nodo " << d.node <<std::endl; 
		operands[0].push_back(d); 
	}
	
	//std::cout << provv.size() << std::endl;
	std::cout << " operands[0] " << operands[0].size() << std::endl; 
	return mign_current.create_maj(operands[0]); 
    
  }

  inline mign_function mign_notequaln(const mign_function_vec_t children)
  {
	  std::cout << " children " << std::endl; 
	  for ( auto& c : children)
	  {
		  std::cout << c.node << std::endl; 
	  }
	
	std::cout << " NOT equal n" << std::endl;
    std::vector<std::vector<mign_function>> operands; 
	operands.resize(4); // operands 0 lower maj 3, operands 1 upper maj3,operands 2 = g,  operands 3 = h h  
	//std::cout << operands.size() << std::endl; 
	const auto k = children.size()-2; 
	
	operands[0].resize(3); 
	operands[1].resize(3); 
	
	operands[0][0] = mign_change_to_smaller(children[0].node)^ children[0].complemented ; 
	auto s = mign_change_to_smaller(children[1].node)^ children[1].complemented ; 
	operands[0][1] = s; 
	
	operands[3].push_back(s); 
	operands[3].push_back(s);
	
	auto provvd =children; 
	
	std::cout << " provvd size = " << provvd.size() << std::endl;
	for ( auto x = 2; x < provvd.size(); x++)
	{
		std::cout << " provv d node " << provvd[x].node << std::endl; 
		auto f = mign_change_to_smaller(provvd[x].node)^ provvd[x].complemented ; 
		std::cout << " enrico 1" << std::endl;
		operands[2].push_back(f);
		std::cout << " enrico 2" << std::endl; 
		operands[3].push_back(f);
	}
		
	auto pair = std::make_pair(0,1); 
	
	auto g = mign_current.create_maj(operands[2]); 
	operands[0][2] = g; 
	operands[1][0] = g; 
	
	auto three = mign_current.create_maj(operands[0]); 
	operands[1][1] = three; 
	
	
	auto h = mign_equaln_inside(operands[3], pair); 
	
	operands[1][2] = h; 
	
	//std::cout << "operands 0 size" << operands[0].size() << std::endl; 
	//std::cout << "operands 1 size" << operands[1].size() << std::endl; 
	//std::cout << "operands 2 size" << operands[2].size() << std::endl; 
	//std::cout << "operands 3 size" << operands[3].size() << std::endl; 
	
	return mign_current.create_maj(operands[1]); 
    
  }
  
  inline mign_function mign_3kdecom(const mign_function_vec_t children)
  {
	  
   // std::cout << " NOT equal n" << std::endl;
    std::vector<std::vector<mign_function>> operands; 
	operands.resize(4); // operands 0 lower maj 3, operands 1 upper maj3,operands 2 = g,  operands 3 = h h  
	//std::cout << operands.size() << std::endl; 
	const auto k = children.size(); 
	
	operands[0].resize(k); 
	operands[1].resize(k); 
	operands[2].resize(k); 
	
	auto s = mign_change_to_smaller(children[0].node)^ children[0].complemented ; 
	operands[0][0] = s; 
	operands[0][1] = s;
	operands[1][2] = s;  
	
	s = mign_change_to_smaller(children[1].node)^ children[1].complemented ; 
	operands[1][0] = s; 
	operands[1][1] = s;
	operands[2][2] = s;  
	
	s = mign_change_to_smaller(children[2].node)^ children[2].complemented ; 
	operands[0][2] = s; 
	operands[2][0] = s;
	operands[2][1] = s;
	
	//std::cout << " provvd size = " << provvd.size() << std::endl;
	for ( auto x = 3; x < children.size(); x++)
	{
		//std::cout << " provv d node " << provvd[x].node << std::endl; 
		auto f = mign_change_to_smaller(children[x].node)^ children[x].complemented ; 
		operands[0][x] = f; 
		operands[1][x] = f;
		operands[2][x] = f;
	}
		
	auto pair = std::make_pair(0,1); 
	
	for ( auto x = 0; x < 3; ++x)
	{
		auto g = mign_equaln_inside(operands[x], pair);
		operands[3].push_back(g); 
	}
		
	return mign_current.create_maj(operands[3]); 
    
  }

public:
  mign_graph                        mign_old;
  mign_graph                        mign_current;
  std::map<mign_node, mign_function> old_to_new;
  bool                             verbose;
  unsigned                          n; 
  
};

/******************************************************************************
 * Private functions                                                          *
 ******************************************************************************/

std::pair<unsigned,unsigned> equal(const mign_function_vec_t& c1 )  // return the position of the two nodes that are equals
{
	auto count1 = 0, count2 = 0, flag = 0, flag_z = 0, store = 0, a = 0, b = 0; 
	
	for ( auto& node : c1)
	{
		count2 = 0; 
		flag_z = 0;
		if (node == 0)
		{
			for ( auto& n : c1)
			{
				if (( flag_z == 0) && (n == 0))
				    {
						flag_z = 1; 
						count2++;
						continue; 
					}
				else if (( flag_z == 1) && (n == 0))
				{
					store = 1; 
					a = count1; 
					b = count2; 
				}
				
				count2++; 
			}
		}
		else 
		{
		flag = 0;
		for ( auto& n : c1)
		{
			if (( flag == 0) && (n.node == node.node))
			    {
					flag = 1; 
					count2++;
					continue; 
				}
			else if (( flag == 1) && (n.node == node.node))
				return std::make_pair(count1,count2); 
			count2++; 
		}
	}
		count1++; 
	}
	if (store == 1)
		return std::make_pair(a,b); 
	
	return std::make_pair(0u,0u); 
}

mign_rewriting_majn_to_smaller_manager::mign_rewriting_majn_to_smaller_manager( const mign_graph& mign, bool verbose, unsigned n)
  : mign_current( mign ),
    verbose( verbose ), 
	n (n)
{
}

void mign_rewriting_majn_to_smaller_manager::swap_current_mign ( const std::string& method )
{
  mign_old = mign_current;
  mign_current = mign_graph();

  old_to_new.clear();
  old_to_new.insert( {0, mign_current.get_constant( false )} );
  
 for ( const auto& input : mign_old.inputs() )
 {
	 auto str = input.second; 
	 old_to_new.insert ({input.first, mign_current.create_pi(str)});
 }
 
}

void mign_rewriting_majn_to_smaller_manager::run_changen_to_smaller()
{
	swap_current_mign( "R" );
	
  for ( const auto& output : mign_old.outputs() )
  {
	  std::cout << " output " << output.first.node << std::endl; 
    mign_current.create_po(mign_change_to_smaller( output.first.node )^ output.first.complemented, output.second );
  }
}

void mign_rewriting_majn_to_smaller_manager::schensted()
{
	swap_current_mign( "S" );
	
  for ( const auto& output : mign_old.outputs() )
  {
    mign_current.create_po(mign_schensted( output.first.node )^ output.first.complemented, output.second );
  }
}

/**
 * 〈〈xyu〉〈xyv〉z〉↦〈xy〈uvz〉〉
 */
mign_function mign_rewriting_majn_to_smaller_manager::mign_change_to_smaller( const mign_node& node )
{
	
  const auto it = old_to_new.find( node );
  if ( it != old_to_new.end() ) { return it->second; }
  
  std::cout << " nodo " << node << std::endl;
  
  std::vector<mign_function> operands; 
  mign_function res; 
  const auto children = mign_old.children(node);
 // std::cout << "children size " << children.size() << std::endl; 
  
  if ( children.size() < n)
	{
		for ( auto a = 0; a < children.size() ; a++)
		{
			operands.push_back(mign_change_to_smaller( children[a].node)^ children[a].complemented ); 
		}
	    res = mign_current.create_maj(operands);
		goto cache_and_return; 
	}
  else
	  {
		  const auto t = equal(children); 
		  // std::cout << " t firs = " << t.first << std::endl; 
		  
		  if (t.first != 0) // per dire se ci sno almen due figli uguali ;) 
		   {
			   res = mign_equaln(children, t);
		       goto cache_and_return;
	       }
		   else 
		   {
			   res = mign_notequaln(children);
		       goto cache_and_return;
		   }
	  } 
	
cache_and_return:
  //std::cout << "olt to new insert "<< node << " with" << res.node << std::endl; 
  old_to_new.insert( {node, res} );
  return res;
}

mign_function mign_rewriting_majn_to_smaller_manager::mign_schensted( const mign_node& node )
{
	
  const auto it = old_to_new.find( node );
  if ( it != old_to_new.end() ) { return it->second; }
  
 // std::cout << " nodo " << node << std::endl;
  
  std::vector<mign_function> operands; 
  mign_function res; 
  const auto children = mign_old.children(node);
 // std::cout << "children size " << children.size() << std::endl; 
  
  if ( children.size() < n)
	{
		for ( auto a = 0; a < children.size() ; a++)
		{
			operands.push_back(mign_schensted( children[a].node)^ children[a].complemented ); 
		}
	    res = mign_current.create_maj(operands);
		goto cache_and_return; 
	}
  else
	  {
			   res = mign_3kdecom(children);
		       goto cache_and_return;
		   
	  } 
	
cache_and_return:
  //std::cout << "olt to new insert "<< node << " with" << res.node << std::endl; 
  old_to_new.insert( {node, res} );
  return res;
}
/******************************************************************************
 * Public functions                                                           *
 ******************************************************************************/

mign_graph mign_rewriting_majn_to_smaller( const mign_graph& mign, const unsigned n, 
                              const properties::ptr& settings,
                              const properties::ptr& statistics )
{
  /* settings */
  const auto verbose = get( settings, "verbose", true);

  /* timer */
  properties_timer t( statistics );
  
  mign_rewriting_majn_to_smaller_manager mgr( mign, verbose, n);
 
  mgr.run_changen_to_smaller(); 
  //mgr.schensted(); 
  
  return mgr.mign_current;
}


}

// Local Variables:
// c-basic-offset: 2
// eval: (c-set-offset 'substatement-open 0)
// eval: (c-set-offset 'innamespace 0)
// End:
