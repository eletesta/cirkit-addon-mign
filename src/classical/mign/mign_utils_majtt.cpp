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

#include "mign_utils_majtt.hpp"

#include <fstream>
#include <unordered_map>
#include <vector>

#include <boost/format.hpp>
#include <classical/mign/mign_simulate.hpp>
#include <classical/mign/math_utils.hpp>

using boost::format;

namespace cirkit
{

/******************************************************************************
 * Types                                                                      *
 ******************************************************************************/


/******************************************************************************
 * Private functions                                                          *
 ******************************************************************************/
	

/******************************************************************************
 * Public functions                                                           *
 ******************************************************************************/

bool all_larger(tt func, tt tt_p)
{		

	for (auto x = 0; x < func.size(); x++)
	{
		if (func[x] == tt_p[x])
			continue; 
		else if ((func[x] == 1) && ( tt_p[x] == 0))
				continue; 
		else 
			return false; 
	}
	return true; 
}

bool all_smaller(tt func, tt tt_p)
{
	tt_align(func, tt_p); 
	for (auto x = 0; x < func.size(); x++)
	{
		if (func[x] == tt_p[x])
			continue; 
		else if ((func[x] == 0) && ( tt_p[x] == 1))
			continue; 
		else 
			return false; 
	}
	return true; 
}

tt calculate_tt_reminder_smaller (tt func, tt tt_p)
{
	tt reminder (func.size()); 
	for (auto x = 0; x < func.size(); x++)
	{
		if ( func[x] == tt_p[x])
			reminder[x] = 0; 
		else 
			reminder[x] = 1; 
	}
	return reminder; 
}

tt calculate_tt_reminder_larger (tt func, tt tt_p)
{
	tt reminder (func.size()); 
	for (auto x = 0; x < func.size(); x++)
	{
		if ( func[x] == tt_p[x])	
			reminder[x] = 1; 
		else 
			reminder[x] = 0; 
	}
	return reminder; 
}

/*  Value = 0 MAJ , value = 1 MAJ con 0 , value = 2 Almost MAJ larger , value = 3 Almost MAJ con 0 larger, value = 4 Almost MAJ smaller , value = 5 Almost MAJ con 0 smaller */
std::tuple<unsigned,boost::dynamic_bitset<>,tt> is_almost_maj (tt func)
{
	auto num_vars = tt_num_vars( func );
    unsigned value = 2; 
	
	if (num_vars %2 == 0)
	{
		tt func_p; 
		func_p = func; 
		tt_resize(func_p, num_vars + 1); 
		
		for (auto x = 0; x < func.size(); x++)
		{
			func_p[x] = func[x]; 
			func_p[x + func.size()] = func[x];		
		}
		tt_resize(func, num_vars + 1); 
		func = func_p; 
		value = 3; 
		num_vars = num_vars + 1;
	}
	
    tt m( 1u << num_vars );
	boost::dynamic_bitset<> it( num_vars, 0 );
	boost::dynamic_bitset<> it_out( num_vars, 0 );
	tt reminder (1u << num_vars); 
	if (value == 3)
	{
		unsigned position = 0u; 
		boost::dynamic_bitset<> it_p( num_vars -1 , 0 );
	    do
	    {
	       m[position] = it.count() > ( num_vars >> 1u );
	       inc( it );
		   inc(it_p);
		   position++; 
	    } while ( it_p.any() ); // create MAJ function 
	    for (auto t= 0; t < position; t++)
		{
			m [t + position] = m [t];
		}
	}
	else 
	{
	    do
	    {
	       m[it.to_ulong()] = it.count() > ( num_vars >> 1u );
	       inc( it );
	    } while ( it.any() ); // create MAJ function 
	}
		do
		{
		tt tt_p( 1u << num_vars );
		tt_p = m; 
		for (auto x = 0; x < it_out.size(); x++)
	      {
			if (it_out[x] == 1)
				tt_p = tt_flip (tt_p,x);
		  }
		  if (all_larger(func,tt_p))
		  {
			  reminder = calculate_tt_reminder_larger (func,tt_p); 
			  return std::make_tuple(value,it_out,reminder); 
		  }
		  else if (all_smaller(func,tt_p))
		  {
			  reminder = calculate_tt_reminder_smaller (func,tt_p); 
			  return std::make_tuple(value+2,it_out,reminder); 
		  }
		  inc(it_out); 
	    } while ( it_out.any() ); 
    
return std::make_tuple(0,func, reminder); 
}

std::pair<unsigned,boost::dynamic_bitset<>> is_maj (tt func)
{
	auto num_vars = tt_num_vars( func );
	unsigned value = 1; 
	
	if (num_vars% 2 == 0)
	{
		tt func_p; 
		func_p = func; 
		tt_extend(func_p, num_vars+ 1); 

		for (auto x = 0; x < func.size(); x++)
		{
			func_p[x] = func[x]; 
			func_p[x + func.size()] = func[x];		
		}
		tt_extend(func, num_vars + 1); 
		func = func_p; 
		value = 2; 
		num_vars = num_vars + 1;	
	}
	tt m( 1u << num_vars );
	boost::dynamic_bitset<> it( num_vars, 0 );
	boost::dynamic_bitset<> it_out( num_vars, 0 );
	    /*do
	    {
	       m[it.to_ulong()] = it.count() > 1;
	       inc( it );
	    } while ( it.any() ); // create AND function
		boost::dynamic_bitset<> it_out_and( num_vars, 0 );
		do
		{ 
		  auto tt_p = m; 
		  for (auto x = 0; x < it_out_and.size(); x++)
	      {
			if (it_out_and[x] == 1)
				tt_p = tt_flip (tt_p,x);
		  }
		  if (func == tt_p)
			  return std::make_pair(1,it_out_and); 
		  
		  inc(it_out_and); 
	    } while ( it_out_and.any() );
	    do
	    {
	       m[it.to_ulong()] = it.count() > 0;
	       inc( it );
	    } while ( it.any() ); // create OR function 
		boost::dynamic_bitset<> it_out_or( num_vars, 0 );
		do
		{ 
		  auto tt_p = m;
		  for (auto x = 0; x < it_out_or.size(); x++)
	      {
			if (it_out_or[x] == 1)
				tt_p = tt_flip (tt_p,x);
		  }
		  if (func == tt_p)
		  {
			  return std::make_pair(2,it_out_or); 
		  }
		  inc(it_out_or); 
	    } while ( it_out_or.any() );*/
	if (value == 2)
	{
		unsigned position = 0u; 
		boost::dynamic_bitset<> it_p( num_vars -1 , 0 );
	    do
	    {
	       m[position] = it.count() > ( num_vars >> 1u );
	       inc( it );
		   inc(it_p);
		   position++; 
	    } while ( it_p.any() ); // create MAJ function 
	    for (auto t= 0; t < position; t++)
		{
			m [t + position] = m [t];
		}
	}
	else 
	{
	    do
	    {
	       m[it.to_ulong()] = it.count() > ( num_vars >> 1u );
	       inc( it );
	    } while ( it.any() ); // create MAJ function 
	}
	   
	do
	{
	auto tt_p = m; 
	for (auto x = 0; x < it_out.size(); x++)
    {
		if (it_out[x] == 1)
			tt_p = tt_flip (tt_p,x);
	}
		if (func == tt_p)
	   {
			return std::make_pair(value,it_out); 
	    }
		inc(it_out); 
    } while ( it_out.any() ); 
	
return std::make_pair(0,func); 
}

}
// Local Variables:
// c-basic-offset: 2
// eval: (c-set-offset 'substatement-open 0)
// eval: (c-set-offset 'innamespace 0)
// End:
