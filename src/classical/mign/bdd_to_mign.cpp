/* CirKit: A circuit toolkit
 * Copyright (C) 2009-2015  University of Bremen
 * Copyright (C) 2015-2017  EPFL
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

#include "bdd_to_mign.hpp"

#include <classical/mign/mign_utils.hpp>

#include <cudd.h>
#include <cuddInt.h>

#include <boost/format.hpp>


namespace cirkit
{
	
mign_function mign_from_add_rec( mign_graph& mign, DdManager* dd, DdNode* node, std::unordered_map<DdNode*,mign_function>& visited )
{
	    if ( visited.find(node) != visited.end() ) { return visited[node]; }
	
		mign_function f; 

	    if ( Cudd_IsConstant( node ) )
	    {
			if (node == DD_ONE(dd))
				f = mign.get_constant(true);
			else 
				f = mign.get_constant(false);
		  return f; 
	    }
	    else
	    {
		  auto index = Cudd_NodeReadIndex( node);
		  std::vector<mign_function> operands; 

	      operands.push_back(mign_from_add_rec( mign, dd, Cudd_E(node), visited));
	      operands.push_back(mign_from_add_rec( mign, dd, Cudd_T(node), visited));
		  operands.push_back(create_mign_function(mign.inputs()[index].first, 0)); 

	      f = mign.create_maj(operands);
	    }
	
		visited.insert( std::make_pair(node,f) );
	    return f; 
	}
	
mign_function mign_from_bdd_rec( mign_graph& mign, DdManager* dd, DdNode* node, std::unordered_map<DdNode*,mign_function>& visited )
{
	auto * r = Cudd_Regular (node); 
	auto is_complement = Cudd_IsComplement(node);
	
	if ( visited.find(node) != visited.end() ) { 
		std::cout << " already visited " << std::endl; 
		return visited[node]; }
	
	mign_function f; 

    if (Cudd_IsConstant( node ) )
    {
		if (!is_complement)
		   f = mign.get_constant(true); 
		else 
	       f = mign.get_constant(false); 
		return f; 
    }
    else
    {
      
      auto index = Cudd_NodeReadIndex (node);
	  std::cout << "node with index = " << index << std::endl; 
	  std::vector<mign_function> operands; 

      operands.push_back(mign_from_bdd_rec( mign, dd, Cudd_E( node )  , visited));
      operands.push_back(mign_from_bdd_rec( mign, dd, Cudd_T( node )  , visited));
	  operands.push_back(create_mign_function(mign.inputs()[index].first, 0)); 

      f = mign.create_maj(operands);
    }
	
	
	visited.insert( std::make_pair(node,f) );
	
    return is_complement ? !f : f;
}

//only single output bdds into MIG-3

mign_graph bdd_to_mign (const bdd_function_t& bdd, const bool ce_on, unsigned order_option)
{
   mign_graph mign; 
   std::unordered_map<DdNode*,mign_function> visited; 
   
   const unsigned n = bdd.first.ReadSize();
   const unsigned m = bdd.second.size();
   
   for (int i = 0; i < m; i++)
   {
	   if (is_monotone( bdd.first, bdd.second[i]) == false)
		   { 
			   continue; 
		   }
	   else 
		   {
			   std::cout << " Function " << i << " not monotone" << std::endl; 
	           return mign; 
	       }
   	
   }
   
   switch ( order_option )
   {
       case 0: 
         break;

       case 4: //kfdd_synthesis_reordering_exact_dtl_friedman:
         Cudd_ReduceHeap( bdd.first.getManager(),CUDD_REORDER_SIFT , 0 );
         break;
		 
	   case 21:
          Cudd_ReduceHeap( bdd.first.getManager(), CUDD_REORDER_EXACT, 0 );
       break;
   }
     // change variables order in the BDD
   
   for ( auto i = 0; static_cast<int>( i ) < n; ++i )
   {
     mign.create_pi( boost::str( boost::format( "x%d" ) % i ) );
   }
   
   if (ce_on == true)
   {
	   std::cout << " ce on true" << std::endl; 
	   for (auto i = 0; static_cast<int>(i) < m; ++i)
	   {
	   	 auto name = boost::str( boost::format( "f%d" ) % i);
		 auto manager = bdd.first.getManager(); 
		 auto f = bdd.second[i].getNode();
	     mign.create_po(mign_from_bdd_rec( mign, manager, f, visited),name);	 
	   }
   }
   else 
   {
	   std::cout << " ce on false" << std::endl; 
	   for (auto i = 0; static_cast<int>(i) < m; ++i)
	   {
	   	 auto name = boost::str( boost::format( "f%d" ) % i);
		 auto manager = bdd.first.getManager(); 
		 const auto add = bdd.second[i].Add();
		 auto f = add.getNode();
	     mign.create_po(mign_from_add_rec( mign, manager, f, visited),name);	 
	   }
   }
   return mign; 	  
}

}