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
#include <classical/mign/mign_from_string.hpp>

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
	tt_is_maj(mign.size()),
	tt_is_almostmaj_or_maj(mign.size()),
	tt_is_reminder(mign.size()),
	tt_is_exact(mign.size()),
    leafs( 1u )
{
}

void mign_cover::add_cut( mign_node n, const mign_cuts_paged::cut& cut, const unsigned T, const std::vector<unsigned> w , const std::vector<bool> N_un)
{
	    assert( offset[n] == 0u );
	    assert( threshold[n] == 0u );


	    neg_un[n] = N_un; 
	    threshold[n]= T; // = T; 
	    weights[n] = w; // = w; 
	    offset[n] = leafs.size();
  
	    leafs.push_back( cut.size() );
         
	    for ( auto l : cut )
	    {
		
	      leafs.push_back( l );
	    }

	    ++count;
}

void mign_cover::add_cut_maj( mign_node n, const mign_cuts_paged::cut& cut, const boost::dynamic_bitset<> tt_maj, const unsigned tt_almostmaj_or_maj, tt reminder, std::string exact)
{

	    offset[n] = leafs.size();
		tt_is_maj[n] = tt_maj; 
		tt_is_almostmaj_or_maj[n] = tt_almostmaj_or_maj; 
		tt_is_reminder[n] = reminder;
		tt_is_exact[n] = exact;  
  
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

boost::dynamic_bitset<> mign_cover::has_tt_maj (mign_node n) const
{

	return tt_is_maj[n]; 
}

unsigned mign_cover::has_tt_almostmaj_or_maj (mign_node n) const
{

	return tt_is_almostmaj_or_maj[n]; 
}
tt mign_cover::has_tt_reminder (mign_node n) const
{

	return tt_is_reminder[n]; 
}

std::string mign_cover::has_tt_exact (mign_node n) const
{
	return tt_is_exact[n]; 
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
   
   if ((mign_old.cover().has_threshold(n) != 0)) 
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

old_to_new.insert( {n, f} );
return f;
}

mign_function mign_to_majn_add_cut( const mign_graph& mign_old, mign_graph& mign_new, mign_node n,
                         std::map<mign_node, mign_function>& old_to_new )
{
	
  std::vector<mign_function> operands; 
  std::vector<mign_function> operands_almost;
  /* visited */
  const auto it = old_to_new.find( n );
  if ( it != old_to_new.end() )
  {
    return it->second;
  }
  
   mign_function f; 
   mign_function f1;
   if ((mign_old.cover().has_tt_almostmaj_or_maj(n) == 0) || (mign_old.cover().has_tt_almostmaj_or_maj(n) == 1))
   {
	   auto leaf_c = 0u; 
	   for ( auto l : mign_old.cover().cut( n ) )
	    {
	     operands.push_back(mign_to_majn_add_cut( mign_old, mign_new, l,old_to_new ) ^ mign_old.cover().has_tt_maj(n)[leaf_c]);
		 leaf_c++; 
	    }
	   if (mign_old.cover().has_tt_almostmaj_or_maj(n) == 1)
	   {
		   assert (leaf_c %2 != 0); 
	   	   operands.push_back(mign_new.get_constant(mign_old.cover().has_tt_maj(n)[leaf_c]));
	   }
	
	   f = mign_new.create_maj(operands); 
   }
   else 
   {
	   auto leaf_c = 0u; 
	   for ( auto l : mign_old.cover().cut( n ) )
	    {
	     operands.push_back(mign_to_majn_add_cut( mign_old, mign_new, l,old_to_new ) ^ mign_old.cover().has_tt_maj(n)[leaf_c]);
		 leaf_c++; 
	    }
	   if (mign_old.cover().has_tt_almostmaj_or_maj(n) %2 != 0)
		    operands.push_back(mign_new.get_constant(mign_old.cover().has_tt_maj(n)[leaf_c]));
	   f1 = mign_new.create_maj(operands); 

		operands_almost.push_back(f1);
		
		/* Exact synthesis */ 
		
		operands_almost.push_back(f1); // qui ci va quello ESATTO - ottenuto con il reminder 
		std::map<char, mign_function> inputs_map; 
		
		
		for (auto l = 0; l < leaf_c; l++)
		{
			char c = 'a' + l; 
			inputs_map.insert(std::make_pair(c, operands[l])); 
		}
	    auto settings = std::make_shared<properties>();
		auto statistics = std::make_shared<properties>();
	    settings->set( "variable_map", inputs_map );
		f1 = mign_from_string( mign_new, mign_old.cover().has_tt_exact(n),settings,statistics ); 
		
		if ((mign_old.cover().has_tt_almostmaj_or_maj(n) == 2) || (mign_old.cover().has_tt_almostmaj_or_maj(n) == 3))
			f = mign_new.create_and(operands_almost); 
		else 
			f = mign_new.create_or(operands_almost); 
   }
   
old_to_new.insert( {n, f} );
return f;
}

mign_function mign_to_threshold_add_cut_multi( mign_graph& mign_old, mign_graph& mign_new, mign_node n,
                         std::map<mign_node, mign_function>& old_to_new , unsigned a)
{

  std::vector<mign_function> operands; 
  /* visited */
  const auto it = old_to_new.find( n );
  if ( it != old_to_new.end() )
  {
    return it->second;
  }
  
   mign_function f; 

  if (mign_old.take_one_cover(a).has_threshold(n) != 0) 
  {
	  
  assert( mign_old.take_one_cover(a).has_cut( n ) );
	  
  auto neg = mign_old.take_one_cover(a).neg_una(n); 
  
  auto leaf_c = 0u; 
  for ( auto l : mign_old.take_one_cover(a).cut( n ) )
  {
    operands.push_back(mign_to_threshold_add_cut_multi( mign_old, mign_new, l,old_to_new,a ) ^ neg[leaf_c]);

	leaf_c++; 
  }
  

	 f = mign_new.create_threshold(operands, mign_old.take_one_cover(a).has_threshold(n), 0, mign_old.take_one_cover(a).has_weights(n)); 
  }

else 
{
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
    }
	
    for ( const auto& output : mign_old.outputs())
    {	 
      mign_new.create_po(mign_to_threshold_add_cut( mign_old, mign_new, output.first.node, old_to_new ) ^ output.first.complemented, output.second );
    }
	   	
	return mign_new; 
}

mign_graph mign_cover_write_maj (const mign_graph& mign_old)  // for almost majority 
{
	
	mign_graph mign_new;
	
    std::map<mign_node, mign_function> old_to_new;

    old_to_new[0] = mign_new.get_constant( false ); 
    for ( const auto& pi : mign_old.inputs() )
    {
      old_to_new[pi.first] = mign_new.create_pi( pi.second);
    }
	
    for ( const auto& output : mign_old.outputs())
    {	 
      mign_new.create_po(mign_to_majn_add_cut( mign_old, mign_new, output.first.node, old_to_new ) ^ output.first.complemented, output.second );
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
