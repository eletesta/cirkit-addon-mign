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

#include "mign_cover.hpp"

#include <fstream>
#include <iostream>

#include <boost/dynamic_bitset.hpp>
#include <boost/format.hpp>

#include <classical/mign/mign_utils.hpp>

namespace cirkit
{


/******************************************************************************
 * mign_cover                                                                  *
 ******************************************************************************/
mign_cover::mign_cover()
{
}

mign_cover::mign_cover( unsigned cut_size, const mign_graph& mign )
  : _cut_size( cut_size ),
    offset( mign.size(), 0u ),
	weights (mign.size()), 
	threshold (mign.size(), 0u), 
	neg_un (mign.size()),
	almost_node (mign.size()),
	almost_node_is_maj (mign.size()),
    leafs( 1u )
{
}

void mign_cover::add_cut( mign_node n, const mign_cuts_paged::cut& cut, const unsigned T, const std::vector<unsigned> w , const std::vector<bool> N_un, const int A_node, const int is_maj)
{
	    assert( offset[n] == 0u );
	    assert( threshold[n] == 0u );
        assert( almost_node[n] == 0 );

	    neg_un[n] = N_un; 
	    threshold[n]= T; // = T; 
	    weights[n] = w; // = w; 
		almost_node[n] = A_node; 
		almost_node_is_maj[n] = is_maj;
	    offset[n] = leafs.size();
  
	    leafs.push_back( cut.size() );
         
	    for ( auto l : cut )
	    {
		
	      leafs.push_back( l );
	    }

	    ++count;
}

bool mign_cover::has_cut( mign_node n) const
{
	
  return offset[n] != 0u;
}

unsigned mign_cover::has_threshold (mign_node n) const
{
	return threshold[n]; 
}

std::vector<unsigned> mign_cover::has_weights (mign_node n) const 
{
	
	return weights[n]; 
}

std::vector<bool> mign_cover::neg_una (mign_node n) const
{

	return neg_un[n]; 
}

int mign_cover::has_almost (mign_node n) const
{
	return almost_node[n]; 
}

int mign_cover::has_is_maj (mign_node n) const
{
	return almost_node_is_maj[n]; 
}

mign_cover::index_range mign_cover::cut( mign_node n) const
{
  return boost::make_iterator_range( leafs.begin() + offset[n] + 1u,
                                     leafs.begin() + offset[n] + 1u + leafs[offset[n]] );
}

mign_function mign_to_threshold_add_cut( const mign_graph& mign_old, mign_graph& mign_new, mign_node n,
                         std::map<mign_node, mign_function>& old_to_new )
{
  std::vector<mign_function> operands; 
  /* visited */
  const auto it = old_to_new.find( n );
  if ( it != old_to_new.end() )
  {
    return it->second;
  }
  
   mign_function f; 
   
   if ((mign_old.cover().has_threshold(n) != 0) || (mign_old.cover().has_almost(n) < 0)) 
   {
	  
   assert( mign_old.cover().has_cut( n ) );
	  
   auto neg = mign_old.cover().neg_una(n); 
  
   auto leaf_c = 0u; 
   for ( auto l : mign_old.cover().cut( n ) )
   {
    operands.push_back(mign_to_threshold_add_cut( mign_old, mign_new, l,old_to_new ) ^ neg[leaf_c]);
	leaf_c++; 
   }
    f = mign_new.create_threshold(operands, mign_old.cover().has_threshold(n), 0, mign_old.cover().has_weights(n)); 
  }

  else if (mign_old.cover().has_threshold(n) == 0)
  {
    const auto c = mign_old.children(n); 
    for ( auto& x : c)
    {
   	   operands.push_back(mign_to_threshold_add_cut( mign_old, mign_new, x.node, old_to_new) ^ x.complemented);
	}

     f  = mign_new.create_maj(operands); 
   }
	
  else if ((mign_old.cover().has_threshold(n) != 0) && (mign_old.cover().has_almost(n) >= 0))
  {
      assert( mign_old.cover().has_cut( n ) );
	  
      auto neg = mign_old.cover().neg_una(n); 
  
      auto leaf_c = 0u; 
      for ( auto l : mign_old.cover().cut( n ) )
      {
       operands.push_back(mign_to_threshold_add_cut( mign_old, mign_new, l,old_to_new ) ^ neg[leaf_c]);
   	   leaf_c++; 
      }
       auto f1 = mign_new.create_threshold(operands, mign_old.cover().has_threshold(n), 0, mign_old.cover().has_weights(n)); 
	   std::vector<mign_function> operands_t; 
	   operands_t.push_back(f1); 
	   mign_node node_a; 
	   node_a = mign_old.cover().has_almost(n); 
	   operands_t.push_back({node_a,0});
	    
	   if (mign_old.cover().has_is_maj(n) == 0)
	   {
		   f = mign_new.create_or(operands_t); 
	   }
	   else if (mign_old.cover().has_is_maj(n) == 1)
	   {
		   f = mign_new.create_and(operands_t); 
	   }
  }

old_to_new.insert( {n, f} );
return f;
}

mign_function mign_to_threshold_add_cut_multi( mign_graph& mign_old, mign_graph& mign_new, mign_node n,
                         std::map<mign_node, mign_function>& old_to_new , unsigned a)
{
	//std::cout << "NODO = " << n << std::endl; 
  std::vector<mign_function> operands; 
  /* visited */
  const auto it = old_to_new.find( n );
  if ( it != old_to_new.end() )
  {
    return it->second;
  }
  
   mign_function f; 
   // can we assure that all circuits will have threshold diversa da 0??
 // std::vector<mign_node> leafs;
  if (mign_old.take_one_cover(a).has_threshold(n) != 0) // ha una threshold diverrsa da 0 IN TEORIA TUTTI OVREBBERO AVERLA divERSA da 0 
  {
	  
  assert( mign_old.take_one_cover(a).has_cut( n ) );
	  
  auto neg = mign_old.take_one_cover(a).neg_una(n); 
  
  auto leaf_c = 0u; 
  for ( auto l : mign_old.take_one_cover(a).cut( n ) )
  {
    operands.push_back(mign_to_threshold_add_cut_multi( mign_old, mign_new, l,old_to_new,a ) ^ neg[leaf_c]);
    //leafs.push_back( l );
	leaf_c++; 
  }
  
	  //std::cout << " cover nodo " << n << " ha th diversa da 0 = " << mign_old.cover().has_threshold(n) << std::endl; 

	 f = mign_new.create_threshold(operands, mign_old.take_one_cover(a).has_threshold(n), 0, mign_old.take_one_cover(a).has_weights(n)); 
  }

else 
{
	//std::cout << " cover nodo " << n << " ha th = 0 = " << mign_old.cover().has_threshold(n) << std::endl; 
    const auto c = mign_old.children(n); 
    for ( auto& x : c)
    {
    	operands.push_back(mign_to_threshold_add_cut_multi( mign_old, mign_new, x.node, old_to_new,a  ) ^ x.complemented);
			
    }

 f  = mign_new.create_maj(operands); 
}

old_to_new.insert( {n, f} );
return f;
}

mign_graph mign_cover_write (const mign_graph& mign_old)
{
	
	mign_graph mign_new;
	
    std::map<mign_node, mign_function> old_to_new;

    old_to_new[0] = mign_new.get_constant( false ); 
    for ( const auto& pi : mign_old.inputs() )
    {
      old_to_new[pi.first] = mign_new.create_pi( pi.second);
	  //node_to_node[pi.first] = name_to_function[name]; 
    }
	
    for ( const auto& output : mign_old.outputs())
    {	 
      mign_new.create_po(mign_to_threshold_add_cut( mign_old, mign_new, output.first.node, old_to_new ) ^ output.first.complemented, output.second );
    }
	   	
	return mign_new; 
}

std::vector<mign_graph> mign_cover_write_multi (mign_graph& mign_old)
{
	
	std::vector<mign_graph> mign_results; 
	
	assert (mign_old.has_multi_cover()); 
	
	for ( auto x = 0; x < mign_old.size_multi_cover(); ++x)
	{
		
	mign_graph mign_new;
	
    std::map<mign_node, mign_function> old_to_new;

    old_to_new[0] = mign_new.get_constant( false ); 
    for ( const auto& pi : mign_old.inputs() )
    {
      old_to_new[pi.first] = mign_new.create_pi( pi.second);
	  //node_to_node[pi.first] = name_to_function[name]; 
    }
	
    for ( const auto& output : mign_old.outputs())
    {
	 
      mign_new.create_po(mign_to_threshold_add_cut_multi( mign_old, mign_new, output.first.node, old_to_new, x) ^ output.first.complemented, output.second );
    }
	
	mign_results.push_back(mign_new); 
   }
	   	
	return mign_results; 
}
/******************************************************************************
 * private functions                                                          *
 ******************************************************************************/

}

// Local Variables:
// c-basic-offset: 2
// eval: (c-set-offset 'substatement-open 0)
// eval: (c-set-offset 'innamespace 0)
// End:
