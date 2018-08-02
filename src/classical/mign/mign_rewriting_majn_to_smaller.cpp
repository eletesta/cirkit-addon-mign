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
  void amarel1();
  void amarel2();
  void optimal_size();
  mign_function mign_change_to_smaller( const mign_node& node );
  mign_function mign_schensted( const mign_node& node );
  mign_function mign_amarel1( const mign_node& node );
  mign_function mign_amarel2( const mign_node& node );
  mign_function mign_optimal_size( const mign_node& node );
  
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
		operands[0].push_back(d); 
	}
	
	return mign_current.create_maj(operands[0]); 
    
  }
  
  inline mign_function mign_equaln_inside(const mign_function_vec_t children, std::pair<unsigned,unsigned> t)
  {
	 
    std::vector<std::vector<mign_function>> operands; 
	const auto e = (children.size()-1)/2; 
	operands.resize(e); 
	
	const auto k = children.size()-2;  // 5 ingressi l'uno 
	
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
		operands[0].push_back(d); 
	}

	return mign_current.create_maj(operands[0]); 
    
  }

  inline mign_function mign_notequaln(const mign_function_vec_t children)
  {
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
	
	for ( auto x = 2; x < provvd.size(); x++)
	{
		auto f = mign_change_to_smaller(provvd[x].node)^ provvd[x].complemented ; 
		operands[2].push_back(f);
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
	
	return mign_current.create_maj(operands[1]); 
    
  }
  
  inline mign_function mign_amareldecomp (const mign_function_vec_t children)
  {
  	std::vector<std::vector<mign_function>> operands; 
	
	const auto k = children.size(); 
	const auto n = (k + 1)/2; 
	const auto op_size = n - 3 + 1; 
	
	operands.resize(op_size + 3); 

	auto count = 0u; 
	
	std::vector<mign_function> s; 
	
	for (auto x = 0; x < k; x++)
	{
		s.push_back(mign_amarel1(children[x].node)^ children[x].complemented); 
	}
		
	
	auto m = n; 
	while (m >= 3)
	{
		for (auto x = 1; x < children.size(); x++)
		{
			if (x == m-1)
				continue; 
			operands[count].push_back(s[x]); 
			//auto s = mign_amarel1(children[x].node)^ children[x].complemented; 
		}
		m--; 
		count++;
	}
	
	for (auto g = 1; g < n ; g++ )
	{
		operands[count].push_back(s[g]); 		
	}
	
	for (auto x = 0; x < op_size; x++)
	{
		operands[count].push_back(mign_current.create_maj(operands[x])); 
	}
	
	auto f = mign_current.create_maj(operands[count]);
	count++; 
	
	for (auto x = 2; x < k; x++)
	{
		operands[count].push_back(s[x]); 
	}
	
	auto f_1 = mign_current.create_maj(operands[count]);
	count++; 
	
	operands[count].push_back(s[0]); 
	operands[count].push_back(f_1); 
	operands[count].push_back(f); 
	
	return mign_current.create_maj(operands[count]); 
	 
  }
  
  inline mign_function mign_seven_exact (const mign_function_vec_t children)
  {
	  
  	std::vector<std::vector<mign_function>> operands; 
	
	const auto k = children.size(); 
	
	operands.resize(7); 

	std::vector<mign_function> s; 
	
	for (auto x = 0; x < k; x++)
	{
		s.push_back(mign_amarel2(children[x].node)^ children[x].complemented); // dovrebbero essere gia tutte fatte 
	}
		
	operands[0].push_back(s[0]); 
	operands[0].push_back(s[2]); 
	operands[0].push_back(s[3]); 
	
	operands[1].push_back(s[4]); 
	operands[1].push_back(s[5]); 
	operands[1].push_back(s[6]);
	
	operands[2].push_back(mign_current.create_maj(operands[0])); 
	operands[2].push_back(s[5]); 
	operands[2].push_back(s[6]);
	
	operands[3].push_back(mign_current.create_maj(operands[1])); 
	operands[3].push_back(s[0]); 
	operands[3].push_back(s[2]);
	
	operands[4].push_back(mign_current.create_maj(operands[0])); 
	operands[4].push_back(mign_current.create_maj(operands[2])); 
	operands[4].push_back(s[4]);
	
	operands[5].push_back(mign_current.create_maj(operands[1])); 
	operands[5].push_back(mign_current.create_maj(operands[3])); 
	operands[5].push_back(s[3]);
	
	operands[6].push_back(mign_current.create_maj(operands[4])); 
	operands[6].push_back(mign_current.create_maj(operands[5])); 
	operands[6].push_back(s[1]);
	
	return mign_current.create_maj(operands[6]); 
	 
  }
  
  inline mign_function mign_nine_aexact (const mign_function_vec_t children)
  {
	  
  	std::vector<std::vector<mign_function>> operands; 
	
	const auto k = children.size(); 
	
	operands.resize(15); 

	std::vector<mign_function> s; 
	
	for (auto x = 0; x < k; x++)
	{
		s.push_back(mign_amarel2(children[x].node)^ children[x].complemented); // dovrebbero essere gia tutte fatte 
	}
		
	operands[0].push_back(s[6] ^ 1); 
	operands[0].push_back(s[7]); 
	operands[0].push_back(s[8]); 
	
	operands[1].push_back(s[4]); 
	operands[1].push_back(s[6]); 
	operands[1].push_back(mign_current.create_maj(operands[0]));
	
	operands[2].push_back(s[6]); 
	operands[2].push_back(s[7]); 
	operands[2].push_back(s[8]);
	
	operands[3].push_back(s[4]); 
	operands[3].push_back(s[5]); 
	operands[3].push_back(mign_current.create_maj(operands[1]));
	
	operands[4].push_back(mign_current.create_maj(operands[1])); 
	operands[4].push_back(mign_current.create_maj(operands[2])); 
	operands[4].push_back(s[5]);
	
	operands[5].push_back(mign_current.create_maj(operands[4])); 
	operands[5].push_back(s[2]);
	operands[5].push_back(s[4]);
	
	operands[6].push_back(mign_current.create_maj(operands[3])); 
	operands[6].push_back(mign_current.create_maj(operands[5])); 
	operands[6].push_back(s[3]);
	
	operands[7].push_back(s[2]);
	operands[7].push_back(s[3]);
	operands[7].push_back(s[5] ^ 1);
	
	operands[8].push_back(s[2]);
	operands[8].push_back(s[3]);
	operands[8].push_back(s[5]);
	
	operands[9].push_back(s[6]);
	operands[9].push_back(mign_current.create_maj(operands[0]));
	operands[9].push_back(mign_current.create_maj(operands[8]));
	
	operands[10].push_back(s[5]);
	operands[10].push_back(mign_current.create_maj(operands[2]));
	operands[10].push_back(mign_current.create_maj(operands[7]));
	
	operands[11].push_back(s[4] ^ 1);
	operands[11].push_back(mign_current.create_maj(operands[9]));
	operands[11].push_back(mign_current.create_maj(operands[10]));
	
	operands[12].push_back(s[4]);
	operands[12].push_back(mign_current.create_maj(operands[9]));
	operands[12].push_back(mign_current.create_maj(operands[10]));
	
	operands[13].push_back(s[0]);
	operands[13].push_back(mign_current.create_maj(operands[6]));
	operands[13].push_back(mign_current.create_maj(operands[11]));
	
	operands[14].push_back(s[1]);
	operands[14].push_back(mign_current.create_maj(operands[12]));
	operands[14].push_back(mign_current.create_maj(operands[13]));
	
	return mign_current.create_maj(operands[14]); 
	 
  }
  
  inline mign_function mign_amarel2decomp (const mign_function_vec_t children)
  {
	  
  	std::vector<std::vector<mign_function>> operands; 
	
	const auto k = children.size(); 
	const auto n = (k + 1)/2; 
	const auto op_size = (n - 1 )* 2 + 1; 
	
	operands.resize(op_size + 2); 

	std::vector<mign_function> s; 
	
	for (auto x = 0; x < k; x++)
	{
		s.push_back(mign_amarel2(children[x].node)^ children[x].complemented); 
	}
		
	auto skip1 = s.size()-1; 
	auto skip2 = s.size()-2;
	for (auto h = 2; h < n; h++)
	{
		for (auto g = 2; g < s.size(); g++)
		{
				if ((g == skip1) || (g == skip2))
					continue; 
				operands[h-2].push_back(s[g]); 
		} 
		skip1 = skip1 - 2; 
		skip2 = skip2 - 2; ; 
	}
	for (auto g = 2; g < s.size(); g++)
	{
			if ((g == 2) || (g == 4))
				continue; 
			operands[n-2].push_back(s[g]); 
	} 
	for (auto g = 2; g < s.size(); g++)
	{
			if ((g == 2) || (g == 3))
				continue; 
			operands[n-1].push_back(s[g]); 
	} 
	
	auto h = 0; 
	for (h = 0; h < n-1; h++)
	{
		if (h == 0)
		{
			operands[n+h].push_back(s[1]); 
			operands[n+h].push_back(mign_current.create_maj(operands[h])); 
			operands[n+h].push_back(mign_current.create_maj(operands[h+1]));
		}
		else 
		{
			operands[n+h].push_back(s[1]); 
			operands[n+h].push_back(mign_current.create_maj(operands[n+h-1])); 
			operands[n+h].push_back(mign_current.create_maj(operands[h+1]));
		}
		
	}
	 
	for (auto e = 2; e < k; e++)
	{
		operands[n+h].push_back(s[e]); 
	}
	
	operands[n+h+1].push_back(s[0]); 
	operands[n+h+1].push_back(mign_current.create_maj(operands[n+h-1])); 
	operands[n+h+1].push_back(mign_current.create_maj(operands[n+h])); 
	
	return mign_current.create_maj(operands[n+h+1]); 
	 
  }
  
  inline mign_function mign_3kdecom(const mign_function_vec_t children)
  {
	
    std::vector<std::vector<mign_function>> operands; 
	operands.resize(4); 
	const auto k = children.size(); 
	
	operands[0].resize(k); 
	operands[1].resize(k); 
	operands[2].resize(k); 
	
	auto s = mign_schensted(children[0].node)^ children[0].complemented ; 
	operands[0][0] = s; 
	operands[0][1] = s;
	operands[1][2] = s;  
	
	s = mign_schensted(children[1].node)^ children[1].complemented ; 
	operands[1][0] = s; 
	operands[1][1] = s;
	operands[2][2] = s;  
	
	s = mign_schensted(children[2].node)^ children[2].complemented ; 
	operands[0][2] = s; 
	operands[2][0] = s;
	operands[2][1] = s;
	
	for ( auto x = 3; x < children.size(); x++)
	{
		auto f = mign_schensted(children[x].node)^ children[x].complemented ; 
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
  mign_graph                          mign_old;
  mign_graph                          mign_current;
  std::map<mign_node, mign_function>  old_to_new;
  bool                                verbose;
  unsigned                            n; 
  
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

void mign_rewriting_majn_to_smaller_manager::amarel1()
{
   swap_current_mign( "A1" );
	
   for ( const auto& output : mign_old.outputs() )
   {
    	mign_current.create_po(mign_amarel1( output.first.node )^ output.first.complemented, output.second );
   }
}

void mign_rewriting_majn_to_smaller_manager::amarel2()
{
   swap_current_mign( "A2" );
	
   for ( const auto& output : mign_old.outputs() )
   {
    	mign_current.create_po(mign_amarel2( output.first.node )^ output.first.complemented, output.second );
   }
}

void mign_rewriting_majn_to_smaller_manager::optimal_size()
{
   swap_current_mign( "OP" );
      
   for (auto& n : mign_old.topological_nodes())
   {
	   mign_optimal_size( n );	   		
   }
   mign_current.create_po(old_to_new[mign_old.outputs()[0].first.node] ^ mign_old.outputs()[0].first.complemented, mign_old.outputs()[0].second );
}

/**
 * 〈〈xyu〉〈xyv〉z〉↦〈xy〈uvz〉〉
 */
mign_function mign_rewriting_majn_to_smaller_manager::mign_change_to_smaller( const mign_node& node )
{
	
  const auto it = old_to_new.find( node );
  if ( it != old_to_new.end() ) { return it->second; }

  
  std::vector<mign_function> operands; 
  mign_function res; 
  const auto children = mign_old.children(node);
  
  if ( children.size() <= n)
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
  old_to_new.insert( {node, res} );
  return res;
}

mign_function mign_rewriting_majn_to_smaller_manager::mign_schensted( const mign_node& node )
{
	
  const auto it = old_to_new.find( node );
  if ( it != old_to_new.end() ) { return it->second; }
    
  std::vector<mign_function> operands; 
  mign_function res; 
  const auto children = mign_old.children(node);
 // std::cout << "children size " << children.size() << std::endl; 
  
  if ( children.size() <= n)
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

mign_function mign_rewriting_majn_to_smaller_manager::mign_amarel1( const mign_node& node )
{
	
  const auto it = old_to_new.find( node );
  if ( it != old_to_new.end() ) { return it->second; }
    
  std::vector<mign_function> operands; 
  mign_function res; 
  const auto children = mign_old.children(node);
   
  if ( children.size() <= n)
  {
		for ( auto a = 0; a < children.size() ; a++)
		{
			operands.push_back(mign_amarel1( children[a].node)^ children[a].complemented ); 
		}
	    res = mign_current.create_maj(operands);
		goto cache_and_return; 
  }
  else
  {
		res = mign_amareldecomp(children);
		goto cache_and_return;
		   
  } 
	
cache_and_return:
  //std::cout << "olt to new insert "<< node << " with" << res.node << std::endl; 
  old_to_new.insert( {node, res} );
  return res;
}

mign_function mign_rewriting_majn_to_smaller_manager::mign_amarel2( const mign_node& node )
{
	
  const auto it = old_to_new.find( node );
  if ( it != old_to_new.end() ) { return it->second; }
    
  std::vector<mign_function> operands; 
  mign_function res; 
  const auto children = mign_old.children(node);
   
  if ( children.size() <= n)
  {
		for ( auto a = 0; a < children.size() ; a++)
		{
			operands.push_back(mign_amarel2( children[a].node)^ children[a].complemented ); 
		}
	    res = mign_current.create_maj(operands);
		goto cache_and_return; 
  }
  else
  {
		res = mign_amarel2decomp(children);
		goto cache_and_return;
		   
  } 
	
cache_and_return:
  old_to_new.insert( {node, res} );
  return res;
}

mign_function mign_rewriting_majn_to_smaller_manager::mign_optimal_size( const mign_node& node )
{ 
  mign_function best_res;
  mign_graph best_mig; 
  unsigned best_size;
  
  const auto it = old_to_new.find( node );
  if ( it != old_to_new.end() ) { return it->second; }
    
  std::vector<mign_function> operands; 
  mign_function res; 
  const auto children = mign_old.children(node);
   
  if ( children.size() <= n)
  {
		for ( auto a = 0; a < children.size() ; a++)
		{
			operands.push_back(mign_optimal_size( children[a].node)^ children[a].complemented ); 
		}
	    res = mign_current.create_maj(operands);
  }
  else if (children.size() == 7)
  {
  	  res = mign_seven_exact(children);
  }
  /*else if (children.size() == 9)
  {
  	  res = mign_nine_aexact(children);
  }*/
  else
  {
	  std::vector<mign_graph> migs(3); 
	  std::vector<mign_function> ress(3); 
	  auto mign_p = mign_current;
	  
	  ress[0] = mign_amarel2decomp(children);
	  migs[0] = mign_current; 
	  
	  migs[1] = mign_p;
	  mign_current = migs[1];
	  
	  ress[1] = mign_amareldecomp(children);
	  migs[1] = mign_current; 
	  
	  migs[2] = mign_p;
	  mign_current = migs[2];
	  
	  ress[2] = mign_3kdecom(children);
	  migs[2] = mign_current; 
	  
	  best_res = ress[0]; 
	  best_mig = migs[0]; 
	  best_size = migs[0].size();
	  
	  for (auto f = 1; f < migs.size(); f++)
	  {
		  if (migs[f].size() < best_size)
		  {
			  best_size = migs[f].size();
			  best_res = ress[f]; 
			  best_mig = migs[f]; 
		  }
	  }
	  
	  res = best_res; 
	  mign_current = best_mig; 	  
  } 
	
  old_to_new.insert( {node, res} );
  return res;
}

/******************************************************************************
 * Public functions                                                           *
 ******************************************************************************/

mign_graph mign_rewriting_majn_to_smaller( mign_graph& mign, const unsigned n, 
                              const properties::ptr& settings,
                              const properties::ptr& statistics )
{
  /* settings */
  const auto verbose = get( settings, "verbose", true);
  const auto option  = get (settings, "option", 0u);

  /* timer */
  properties_timer t( statistics );
  
  mign_rewriting_majn_to_smaller_manager mgr ( mign, verbose, n);
 
  if (option == 0)
  {
	  mgr.run_changen_to_smaller();  
  }
	  	 
  else if (option == 1)
  {
	  mgr.schensted(); 
  }
  
  else if (option == 2)
  {
  	 mgr.amarel1();  
  }
 
  else if (option == 3)
  {
  	 mgr.amarel2();  
  }
  
  else if (option == 4)
  {
	  if (mign.outputs().size() > 1)
		  {
			  std::cout << " Error! Works for single-output" << std::endl; 
			  return mign;
		  }
  	 mgr.optimal_size();  
  }
 	  
  else 
	  {
	  	std::cout << " not valid option" << std::endl; 
		return mign;
	  } 
  
  return mgr.mign_current;

}


}

// Local Variables:
// c-basic-offset: 2
// eval: (c-set-offset 'substatement-open 0)
// eval: (c-set-offset 'innamespace 0)
// End:
