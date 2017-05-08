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

#include "mign_rewriting_maj3.hpp"

#include <functional>

#include <boost/assign/std/vector.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/format.hpp>
#include <boost/range/iterator_range.hpp>

#include <core/graph/depth.hpp>
#include <core/utils/timer.hpp>
#include <classical/mign/mign_utils.hpp>
#include <classical/mign/mign_rewrite.hpp>

using namespace boost::assign;

namespace cirkit
{

/******************************************************************************
 * Types                                                                      *
 ******************************************************************************/

using mign_function_vec_t = std::vector<mign_function>;

struct children_mign_pair_t
{
  unsigned left_a;
  unsigned left_b;
  unsigned right_a;
  unsigned right_b;
};

inline std::pair<unsigned, unsigned> three_without( unsigned x )
{
  return std::make_pair( x == 0u ? 1u : 0u, x == 2u ? 1u : 2u );
}

mign_function make_migt_function (mign_function f, bool complemented)
{
	if (complemented == 1) 
	return !f; 
	else 
		return f; 
}

class mign_rewriting_maj3_manager
{
public:
  mign_rewriting_maj3_manager( const mign_graph& mign, bool verbose );

  void swap_current( const std::string& method );
  //inline unsigned depth() const { return max_depth; }

  void run_mign_distributivity();
  void run_mign_associativity();

  mign_function mign_distributivity( const mign_node& f );
  mign_function mign_associativity( const mign_node& f );

private:
  inline mign_function mign_distributivity_apply( const children_mign_pair_t& pair,
                                                const mign_function_vec_t& children_a,
                                                const mign_function_vec_t& children_b,
                                                const mign_function& other )
  {
    assert( children_a[pair.left_a] == children_b[pair.right_a] );
    assert( children_a[pair.left_b] == children_b[pair.right_b] );

    const auto x = pair.left_a;
    const auto y = pair.left_b;
    const auto u = 3u - pair.left_a - pair.left_b;
    const auto v = 3u - pair.right_a - pair.right_b;

    if ( verbose )
    {
      std::cout << boost::format( "[i] x: %d, y: %d, u: %d, v: %d") % x % y % u % v << std::endl;
    }
	  //std::cout << " mign _ apply " << std::endl; 
	std::vector<mign_function> operands, operands_two; 
	operands.push_back(make_migt_function( mign_distributivity( children_a[x].node ), children_a[x].complemented ));
	operands.push_back(make_migt_function( mign_distributivity( children_a[y].node ), children_a[y].complemented ));
	operands_two.push_back(make_migt_function( mign_distributivity( children_a[u].node ), children_a[u].complemented ));
	operands_two.push_back(make_migt_function( mign_distributivity( children_b[v].node ), children_b[v].complemented ));
	operands_two.push_back(make_migt_function( mign_distributivity( other.node ), other.complemented ) ); 
	operands.push_back(mign_current.create_maj(operands_two)); 
	return mign_current.create_maj(operands); 
    
  }

  inline mign_function mign_associativity_apply( const mign_function_vec_t& grand_children,
                                           const mign_function& common,
                                           const mign_function& extra )
  {
	  std::vector<mign_function> operands; 
	  std::vector<mign_function> operands_two; 
	  
	  	const auto common_f = make_migt_function( mign_associativity( common.node ), common.complemented );
		operands.push_back(common_f); 
		operands_two.push_back(common_f); 
	  
   
	operands.push_back(make_migt_function( mign_associativity( grand_children[0u].node ), grand_children[0u].complemented )); 
	operands_two.push_back(make_migt_function( mign_associativity( grand_children[1u].node ), grand_children[1u].complemented )); 
	operands_two.push_back(make_migt_function( mign_associativity( extra.node ), extra.complemented)); 
	operands.push_back(mign_current.create_maj(operands_two)); 
	
	return mign_current.create_maj(operands); 
   
  }
  
  /*8
   * @return (a,b) where z = node[a][b]
   */
  //boost::optional<std::pair<unsigned, unsigned>> find_depth_distributivity_candidate( const mign_node& node ) const;

  /**
   * @return (a,b,c) where x = node[a] and z = node[b][c]
   */
  //boost::optional<std::tuple<unsigned, unsigned, unsigned>> find_depth_associativity_candidate( const mign_node& node ) const;


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

bool is_regular_nonterminal( const mign_graph& mign, const mign_function& f )
{
	if ((f.complemented == 0) && (mign.fanin_count(f.node) != 0))
		return 1; 
	else return 0; 
}

boost::dynamic_bitset<> get_pair_pattern( const mign_function_vec_t& c1,
                                          const mign_function_vec_t& c2 )
{
  boost::dynamic_bitset<> equals( 9u );

  equals.set( 0u, c1[0u] == c2[0u] );
  equals.set( 1u, c1[0u] == c2[1u] );
  equals.set( 2u, c1[0u] == c2[2u] );
  equals.set( 3u, c1[1u] == c2[0u] );
  equals.set( 4u, c1[1u] == c2[1u] );
  equals.set( 5u, c1[1u] == c2[2u] );
  equals.set( 6u, c1[2u] == c2[0u] );
  equals.set( 7u, c1[2u] == c2[1u] );
  equals.set( 8u, c1[2u] == c2[2u] );

  return equals;
}

std::vector<children_mign_pair_t> get_children_pairs( const mign_function_vec_t& c1,
                                                 const mign_function_vec_t& c2 )
{
  std::vector<children_mign_pair_t> pairs;

  const auto pattern = get_pair_pattern( c1, c2 );

  auto pos_a = pattern.find_first(); // find the first 1 

  while ( pos_a != boost::dynamic_bitset<>::npos )
  {
    auto pos_b = pattern.find_next( pos_a );
    while ( pos_b != boost::dynamic_bitset<>::npos )
    {
      pairs.push_back({(unsigned)pos_a / 3u,
                       (unsigned)pos_b / 3u,
                       (unsigned)pos_a % 3u,
                       (unsigned)pos_b % 3u} );
      pos_b = pattern.find_next( pos_b );
    }
    pos_a = pattern.find_next( pos_a );
  }

  return pairs;
}


mign_rewriting_maj3_manager::mign_rewriting_maj3_manager( const mign_graph& mign, bool verbose )
  : mign_current( mign ),
    verbose( verbose )
{
}

void mign_rewriting_maj3_manager::swap_current( const std::string& method )
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

void mign_rewriting_maj3_manager::run_mign_distributivity()
{
	swap_current( "D" );

  for ( const auto& output : mign_old.outputs() )
  {
    mign_current.create_po( make_migt_function( mign_distributivity( output.first.node ), output.first.complemented ), output.second );
  }
}

void mign_rewriting_maj3_manager::run_mign_associativity()
{
	swap_current( "A" );

  for ( const auto& output : mign_old.outputs() )
  {
    mign_current.create_po( make_migt_function( mign_associativity( output.first.node ), output.first.complemented ), output.second );
  }
}


/**
 * 〈〈xyu〉〈xyv〉z〉↦〈xy〈uvz〉〉
 */
mign_function mign_rewriting_maj3_manager::mign_distributivity( const mign_node& node )
{
  //node is terminal 
  const auto it = old_to_new.find( node );
  if ( it != old_to_new.end() ) { return it->second; }

  std::vector<mign_function> operands; 
  const auto children = mign_old.children(node);
  mign_function_vec_t children_a, children_b, children_c;
  mign_function res;

  mign_old.compute_fanout(); 

  if ( is_regular_nonterminal( mign_old, children[0u] ) && mign_old.fanout_count(children[0u].node) == 1u )
  {
    // first child is not terminal 
    children_a = mign_old.children(children[0u].node);

    // check (0,1) 

    if ( is_regular_nonterminal( mign_old, children[1u] ) && mign_old.fanout_count(children[1u].node) == 1u )
    {
      children_b = mign_old.children(children[1u].node);

      const auto pairs = get_children_pairs( children_a, children_b );

	  
      if ( !pairs.empty() )
      {
        res = mign_distributivity_apply( pairs.front(), children_a, children_b, children[2u] );
        goto cache_and_return;
      }
    }

    // check (0,2) 
    if ( is_regular_nonterminal( mign_old, children[2u] ) && mign_old.fanout_count(children[2u].node) == 1u )
    {
      children_c = mign_old.children(children[2u].node);

      const auto pairs = get_children_pairs( children_a, children_c );

      if ( !pairs.empty() )
      {
        res = mign_distributivity_apply( pairs.front(), children_a, children_c, children[1u] );
        goto cache_and_return;
      }
    }
  }
  //std::cout << " qui devi entrare " << std::endl;
  // check (1,2) 
  if ( is_regular_nonterminal( mign_old, children[1u] ) && is_regular_nonterminal( mign_old, children[2u] ) &&
      (mign_old.fanout_count(children[1u].node) == 1u) && (mign_old.fanout_count(children[2u].node) == 1u))
  {
	 // std::cout << " qui devi entrare " << std::endl; 
    if ( children_b.empty() ) { children_b = mign_old.children(children[1u].node); }
    if ( children_c.empty() ) { children_c = mign_old.children(children[2u].node); }

    const auto pairs = get_children_pairs( children_b, children_c );

    if ( !pairs.empty() )
    {
      res = mign_distributivity_apply( pairs.front(), children_b, children_c, children[0u] );
      goto cache_and_return;
    }
  }

	 for ( auto a = 0; a < children.size(); ++a)
	 {
	 	operands.push_back(make_migt_function( mign_distributivity( children[a].node ), children[a].complemented )); 
	 }
    res = mign_current.create_maj(operands); 
	
	
cache_and_return:
  old_to_new.insert( {node, res} );
  return res;
}

/**
 * 〈xu〈yuz〉〉↦〈zu〈yux〉〉
 */
mign_function mign_rewriting_maj3_manager::mign_associativity( const mign_node& node )
{
  //node is terminal 
  const auto it = old_to_new.find( node );
  if ( it != old_to_new.end() ) { return it->second; }

  const auto children = mign_old.children(node);
  mign_function res;
  std::vector<mign_function> operands; 
  auto s = children.size(); 

  mign_old.compute_fanout(); 
  
  for ( auto i = 0u; i < 3u; ++i )
  {
    if ( is_regular_nonterminal( mign_old, children[i] ) && mign_old.fanout_count(children[i].node) == 1u )
    {
      for ( auto j = 0u; j < 3u; ++j )
      {
        if ( i == j ) { continue; }

        auto grand_children = mign_old.children(children[i].node);

        const auto it = boost::find( grand_children, children[j] );
        if ( it != grand_children.end() )
        {
          grand_children.erase( it );
          res = mign_associativity_apply( grand_children, children[j], children[3u - j - i] );
          goto cache_and_return;
        }
      }
    }
  }
    /* recur */
 
	for (auto x = 0; x < s; ++x)
	{
		operands.push_back(make_migt_function( mign_associativity( children[x].node ), children[x].complemented )); 
	}
    res = mign_current.create_maj(operands); 
	
cache_and_return:
  old_to_new.insert( {node, res} );
  return res;
}

/******************************************************************************
 * Public functions                                                           *
 ******************************************************************************/

mign_graph mign_rewriting_maj3( const mign_graph& mign,
                              const properties::ptr& settings,
                              const properties::ptr& statistics )
{
  /* settings */
  const auto effort  = get( settings, "effort",  1u );
  const auto verbose = get( settings, "verbose", true);

  /* timer */
  properties_timer t( statistics );
  
  mign_rewriting_maj3_manager mgr( mign, verbose );

  for ( auto k = 0u; k < effort; ++k )
  {
	mgr.run_mign_distributivity();
	mgr.run_mign_associativity();
  }
  
  //set( statistics, "distributivity_mign_count",       mgr.distributivity_count );
  //set( statistics, "associativity_mign_count",        mgr.associativity_count );
  
  return mgr.mign_current;
}


}

// Local Variables:
// c-basic-offset: 2
// eval: (c-set-offset 'substatement-open 0)
// eval: (c-set-offset 'innamespace 0)
// End:
