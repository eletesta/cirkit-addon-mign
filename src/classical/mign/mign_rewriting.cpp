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

#include "mign_rewriting.hpp"

#include <functional>

#include <boost/assign/std/vector.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/format.hpp>
#include <boost/range/iterator_range.hpp>
#include <boost/range/algorithm_ext/iota.hpp>

#include <core/graph/depth.hpp>
#include <core/utils/timer.hpp>
#include <core/utils/combinations.hpp>

#include <classical/mign/mign_utils.hpp>
#include <classical/mign/mign_rewrite.hpp>
#include <classical/mign/math_utils.hpp>

using namespace boost::assign;

namespace cirkit
{

/******************************************************************************
 * Types                                                                      *
 ******************************************************************************/

using mign_function_vec_t = std::vector<mign_function>;


mign_function make_mign_function (mign_function f, bool complemented)
{
	if (complemented == 1) 
	return !f; 
	else 
		return f; 
}

class mign_rewriting_manager
{
public:
  mign_rewriting_manager( const mign_graph& mign, bool verbose );

  void swap_current( const std::string& method );
  //inline unsigned depth() const { return max_depth; }

  void run_mign_distributivity();
  void run_mign_associativity();

  mign_function mign_distributivity( const mign_node& f );
  mign_function mign_associativity( const mign_node& f );
  bool initial_check( const std::vector<mign_function> children );

private:
  inline mign_function mign_distributivity_apply( const std::vector<std::vector<unsigned>> uguali,
                                                const std::vector<mign_function_vec_t>& children,
                                                const std::vector<mign_function>& other )
  {
	std::vector<mign_function> operands, operands_two; 
	
	auto count = 0; 
	for ( auto &c:children)
	{
		for (auto x = 0; x <c.size(); ++x)
		{
			const auto it = std::find(uguali[count].begin(), uguali[count].end(), x);  
			if ( it == uguali[count].end() )
			{
				operands_two.push_back(make_mign_function( mign_distributivity( c[x].node ), c[x].complemented ));
				//std::cout <<" aggiungiamo il nodo " << x << " figlio di " << count << std::endl; 
			}
		}
		++count; 
	}
	
	for (auto x = 0; x < other.size(); ++x)
	{
		operands_two.push_back(make_mign_function( mign_distributivity( other[x].node ), other[x].complemented ) ); 
    }
	operands.push_back(mign_current.create_maj(operands_two)); 
	
	for (auto x = 0; x <children[0].size(); ++x)
	{
		const auto it = std::find(uguali[0].begin(), uguali[0].end(), x);  
		if ( it != uguali[0].end() )
			operands.push_back(make_mign_function( mign_distributivity( children[0][x].node ), children[0][x].complemented ));
	}
	
	return mign_current.create_maj(operands); 
    
  }

  inline mign_function mign_associativity_apply( const mign_function_vec_t& grand_children,
                                           const std::vector<mign_function>& common,
                                           const mign_function& extra )
  {
	  std::vector<mign_function> operands; 
	  std::vector<mign_function> operands_two; 
	  for (auto x = 0; x <common.size(); ++x)
	  {
	  	const auto common_f = make_mign_function( mign_associativity( common[x].node ), common[x].complemented );
		operands.push_back(common_f); 
		operands_two.push_back(common_f); 
	  }
   
	operands.push_back(make_mign_function( mign_associativity( grand_children[0u].node ), grand_children[0u].complemented )); 
	operands_two.push_back(make_mign_function( mign_associativity( grand_children[1u].node ), grand_children[1u].complemented )); 
	operands_two.push_back(make_mign_function( mign_associativity( extra.node ), extra.complemented)); 
	operands.push_back(mign_current.create_maj(operands_two)); 
	
	return mign_current.create_maj(operands); 
   
  }
 


public:
  mign_graph                        mign_old;
  mign_graph                        mign_current;
  std::map<mign_node, mign_function> old_to_new;
  bool                             verbose;
  std::vector<unsigned>            indegree;
  std::vector<unsigned>            depths;
  unsigned                         max_depth;

  bool                             use_distributivity       = true;
  bool                             use_associativity        = true;
  /* statistics */
  unsigned                         distributivity_count       = 0u;
  unsigned                         associativity_count        = 0u;
};

/******************************************************************************
 * Private functions                                                          *
 ******************************************************************************/

inline bool is_terminal( const mign_graph& mign, const mign_function& f )
{
  return mign.fanin_count(f.node) == 0u;
}

inline bool is_regular_nonterminal( const mign_graph& mign, const mign_function& f )
{
	if ((f.complemented == 0) && (mign.fanin_count(f.node) != 0))
		return 1; 
	else return 0; 
}

std::vector<boost::dynamic_bitset<>> get_pair_pattern( const std::vector<mign_function_vec_t>& c1)
{
	//std::cout << " get pattern " << std::endl; 
  const auto n = c1[0u].size();  
  const auto enne = c1.size(); 
  auto f = fact(enne-1); 
  std::vector<boost::dynamic_bitset<>> equals(f);
  
  auto count_c = 1;
  auto count_go = 0u;  
  for (auto & c : c1)
  {
	  //std::cout << " siamo al nodo " << count_c - 1 << std::endl; 
	  for (auto x = count_c; x <= (enne - 1); ++x)
	  {
		  
		  for (auto z = 0; z < n; ++z)
		  {
			  //std::cout <<  " figlio " << z << std::endl; 
			  for ( auto y = 0; y<n; ++y)
			  {
				  
			  	equals[count_go].push_back(c[z] == c1[x][y]);
				
			  }
			 
		  }
		  ++count_go;
	  }
	  ++count_c; 
  }
  
  return equals;
}

std::vector<std::vector<unsigned>> get_children_groups( std::vector<mign_function_vec_t>& c1)
{
  
	const auto enne = c1.size(); 
  const auto n = c1[0u].size(); 
  auto f = fact(enne-1); 
  
  std::vector<std::vector<unsigned>> uguali(enne);
  //std::cout << "calcolo patern" << std::endl; 
  const auto pattern = get_pair_pattern(c1);
  auto pattern_provv = pattern; 
  auto uguali_n = 0u; 
  // tutti uguali ;)
  for ( auto & line : pattern)
  {
  	auto pos_a = line.find_first();
	auto count = 0u; 
	 while ( pos_a != boost::dynamic_bitset<>::npos )
	 {
	 	pos_a = line.find_next(pos_a);
		++count; 
	 }
	 //std::cout << " cunt" << count << std::endl; 
	 if (count == (n-1))
	 {
		 ++uguali_n; 
	 }
  }
  
 // std::cout << " f " << f << std::endl; 
 // std::cout << " uguali " << uguali_n << std::endl; 
  
  if (uguali_n == f)
  {
	  //std::cout << " uguale!!!" << std::endl; 
	  for (auto x = 1; x <=(enne - 1); ++x)
	  {
		  if ( x == 1) // devo salvare anche quelli dello 0
		  {
			 // std::cout << " caso x = 1" << std::endl; 
			  auto pos_a = pattern[x-1].find_first(); 
		      while ( pos_a != boost::dynamic_bitset<>::npos )
		      {
				  uguali[0].push_back((unsigned)pos_a / n); 
				  uguali[x].push_back((unsigned)pos_a % n); 
				  pos_a = pattern[x-1].find_next( pos_a );
		      }
		  }
		  else 
		  {
			  auto pos_a = pattern[x-1].find_first(); 
		      while ( pos_a != boost::dynamic_bitset<>::npos )
		      {
				  //uguali[0].push_back((unsigned)pos_a / n); 
				  uguali[x].push_back((unsigned)pos_a % n); 
				  pos_a = pattern[x-1].find_next( pos_a );
		      }
		  }
	  }
  }
  
  else 
  {
	  uguali.erase(uguali.begin()); 
  }
  
  //std::cout << " uguali size " << std::endl; 
  //std::cout << " uguali " << std::endl; 
  return uguali;
}


mign_rewriting_manager::mign_rewriting_manager( const mign_graph& mign, bool verbose )
  : mign_current( mign ),
    verbose( verbose )
{
}

void mign_rewriting_manager::swap_current( const std::string& method )
{
  mign_old = mign_current;
  mign_current = mign_graph();

  old_to_new.clear();
  old_to_new.insert( {0, mign_current.get_constant( false )} );
  //old_to_new.insert( {"1'b1", mign_current.get_constant( true )} );


 for ( const auto& input : mign_old.inputs() )
 {
	 auto str = input.second; 
	 
	 old_to_new.insert ({input.first, mign_current.create_pi(str)});
 }
 
  /* depth */
 
  auto max_depth = evaluate_depth( mign_old);

  if ( verbose )
  {
	   std::cout << "[i] current depth: " << max_depth << ", run " << method << std::endl;
  }
  
}

void mign_rewriting_manager::run_mign_distributivity()
{
	swap_current( "D" );

  for ( const auto& output : mign_old.outputs() )
  {
    mign_current.create_po( make_mign_function( mign_distributivity( output.first.node ), output.first.complemented ), output.second );
  }
}

void mign_rewriting_manager::run_mign_associativity()
{
	swap_current( "A" );

  for ( const auto& output : mign_old.outputs() )
  {
    mign_current.create_po( make_mign_function( mign_associativity( output.first.node ), output.first.complemented ), output.second );
  }
}

bool mign_rewriting_manager::initial_check( const std::vector<mign_function> children )
{
	auto s = children.size(); 
	for (auto & c:children)
	{
		auto x = mign_old.children(c.node); 
		if (x.size() == s)
			continue; 
		else return false; 
	}
	return true; 
}
/**
 * 〈〈xyu〉〈xyv〉z〉↦〈xy〈uvz〉〉
 */
mign_function mign_rewriting_manager::mign_distributivity( const mign_node& node )
{
	
    const auto it = old_to_new.find( node );
    if ( it != old_to_new.end() ) { return it->second; }

    std::vector<mign_function> operands; 
    const auto children = mign_old.children(node);
	 mign_function res;
	 
    if (initial_check(children)) // node and its children of the same size!! 
	{
   
	int limit = children.size()/2; 
	//std::cout << " limit = " << limit << std::endl; 

    mign_old.compute_fanout(); 
	
	signed value = children.size()-2; 
	while(value >= limit)
	{
	
	for ( auto j = 0; j <children.size(); ++j)
	{
	//std::cout << " siamo al figlio " << j << " del nodo " << node << std::endl; 
	auto children_provv = mign_old.children(node);
	const signed n = children_provv.size() - 1; 
	children_provv.erase(children_provv.begin() + j); 
  
	
    if ( is_regular_nonterminal( mign_old, children[j] ) && mign_old.fanout_count(children[j].node) == 1u )
    {
		
	    std::vector<unsigned> numbers( children.size());
	    boost::iota( numbers, 0u );
		//std::cout << " nodo j " << j << " is reguNT" << std::endl; 
		const auto it = find(numbers.begin(),numbers.end(),j); 
		numbers.erase(it); // elimino il nodo j da quelli su cui fare le combinazioni 
		
		
		// coefficiente binomiale n! / k! (n-k)!
		auto fact_value = fact(value); 
		//std::cout << fact_value << std::endl; 
		auto fact_size = fact(n); 
		//std::cout << fact_size << std::endl; 
		auto fact_nk = fact(n-value); 
		auto number_c = fact_size/(fact_value*fact_nk); 
		//std::cout << " num" << number_c << std::endl; 
		std::vector<std::vector<unsigned>> combinations(number_c); 
		std::vector<std::vector<unsigned>> excluded_c(n-value); 
		auto count = 0u, counting = 0u; 
		
		
	    boost::unofficial::for_each_combination( numbers.begin(), numbers.begin() + value, numbers.end(),
	                                             [&combinations,&count, &counting, &value]( decltype( numbers )::const_iterator first,
	                                                 decltype( numbers )::const_iterator last ) {
	                                                
												//count = 0u; 
													
	                                               while ( first != last )
	                                               {
													   if (( counting != 0) && (counting%value == 0))
													   {
														   ++count; 
														   counting = 0; 
													   }
													   	combinations[count].push_back(*first);
														++counting; 
													   *first++; 
	                                               }
	                                
	                                               return false;
	                                             } );
												
		for ( auto & element : combinations)										 
		{
			std::vector<mign_function_vec_t> children_abc;
			children_abc.push_back(mign_old.children(children[j].node));
			auto excluded = children_provv; 
			std::vector<mign_function> combination; 
			for (auto x = 0; x < element.size(); ++x)
				{
					//std::cout << " eleemnt " << element[x] << std::endl; 
					combination.push_back(children[element[x]]); 
					const auto it = find(excluded.begin(), excluded.end(),children[element[x]]); 
					//assert (it != excluded.end()); 
						excluded.erase(it); 
				}	
				
				
					for (auto& boh:combination)
					{
				        if ( is_regular_nonterminal( mign_old, boh ) && mign_old.fanout_count(boh.node) == 1u )
				        {
				           children_abc.push_back(mign_old.children(boh.node));
				  	    }
	  
					}	
					
					//std::cout << "children abc size " << children_abc.size() << std::endl; 
					//std::cout << " value + 1" << value + 1 << std::endl; 
		  	  	  if (children_abc.size() == value + 1 ) // significa che sono tutti nodini che vanno bene 
		  			{
		  				//std::cout << " uttti i nodi vanno bene e cacolo groups " << std::endl; 
		  			const auto groups = get_children_groups( children_abc);
	
		  	        if ( groups.size() == value + 1 )
		  	        {
		  				//std::cout << " e empty? " << std::endl; 
		  	          res = mign_distributivity_apply( groups, children_abc, excluded );
		  	          goto cache_and_return;
		  	        }
		  		}								 		
		}
		 
  }
  
}
value = value - 1;
}
   }
   // recur 
	 for ( auto a = 0; a < children.size(); ++a)
	 {
	 	operands.push_back(make_mign_function( mign_distributivity( children[a].node ), children[a].complemented )); 
	 }
    res = mign_current.create_maj(operands); 
	
	
cache_and_return:
  old_to_new.insert( {node, res} );
  return res;
}

/**
 * 〈xu〈yuz〉〉↦〈zu〈yux〉〉
 */
mign_function mign_rewriting_manager::mign_associativity( const mign_node& node )
{
    const auto it = old_to_new.find( node );
    if ( it != old_to_new.end() ) { return it->second; }

    const auto s = mign_old.children(node).size();
	const auto children_or = mign_old.children(node);
    mign_function res;
    std::vector<mign_function> operands; 

    mign_old.compute_fanout(); 
	
	if (initial_check(children_or))
	{
  
    for ( auto i = 0u; i < s; ++i )
    {
	  auto children = mign_old.children(node);
	  auto grand_children = mign_old.children(children[i].node);
	  std::vector<mign_function> store_children; 
	  auto count = 0u; 
	 
	  
      if ( is_regular_nonterminal( mign_old, children[i] ) && mign_old.fanout_count(children[i].node) == 1u )
      {
        for ( auto j = 0u; j < children.size(); ++j )
        {
          if ( i == j ) { continue; 	  
		  }
		 
          const auto it = boost::find( grand_children, children[j] );
          if ( it != grand_children.end() )
          {
            grand_children.erase( it );
			store_children.push_back(children[j]);
			++count;
	 	  }
	      
          }
	
		  if (count == s - 2)
		  {
			  const auto iter = boost::find( children, children[i] );
			  children.erase(iter); 
			  for (auto x = 0; x < store_children.size(); ++x)
			  {
				  const auto itera = boost::find( children, store_children[x]);
				  children.erase(itera); 
			  }
			  assert (children.size() == 1);
              res = mign_associativity_apply( grand_children, store_children, children[0] );
              goto cache_and_return;
		  }
        }
      }
  }  
   
    /* recur */
	for (auto x = 0; x < s; ++x)
	{
		operands.push_back(make_mign_function( mign_associativity( children_or[x].node ), children_or[x].complemented )); 
	}
    res = mign_current.create_maj(operands); 
	
cache_and_return:
  old_to_new.insert( {node, res} );
  return res;
}

/******************************************************************************
 * Public functions                                                           *
 ******************************************************************************/


mign_graph mign_rewriting( const mign_graph& mign,
                              const properties::ptr& settings,
                              const properties::ptr& statistics )
{
  /* settings */
  const auto effort  = get( settings, "effort",  3u );
  const auto verbose = get( settings, "verbose", true);

  /* timer */
  properties_timer t( statistics );
  
  mign_rewriting_manager mgr( mign, verbose );

  for ( auto k = 0u; k < effort; ++k )
  {
	  mgr.run_mign_distributivity();
	  mgr.run_mign_associativity();
  }
  
  //set( statistics, "distributivity_mign_count",       mgr.distributivity_count );
  //set( statistics, "associativity_mign_count",        mgr.associativity_count );
  
  mign_graph mign_new; // to ensure that we are removing all useless nodes ;) 
  mign_new = mign_rewrite_top_down(mgr.mign_current,settings,statistics); 
  
  return mign_new;
}


}

// Local Variables:
// c-basic-offset: 2
// eval: (c-set-offset 'substatement-open 0)
// eval: (c-set-offset 'innamespace 0)
// End:
