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

#include "mig_plasmon_area.hpp"

#include <fstream>
#include <iostream>

#include <boost/optional/optional_io.hpp>
#include <boost/format.hpp> 
#include <classical/mig/mig_utils.hpp>

using boost::format;

namespace cirkit 
{
	

unsigned mig_plasmon_area ( mign_graph& mign)
{
	//std::cout << " area" << std::endl;
	//std::pair<unsigned, unsigned> costs;  
	unsigned area = 0u; 
	//unsigned delay = 0u; 
	mign.compute_levels(); 
	mign.compute_fanout(); 
	
	std::cout << " [i] Area"<< std::endl; 
	for ( auto v : mign.topological_nodes() )
    {
	   if ( mign.is_input( v ) ) { continue; }
	   auto c = mign.children(v); 
	   
	  if (!mign.is_output(v))
	  {
	   if (mign.level(v) == 1)
		  {
			  std::cout << " level 1" << mign.fanout_count(v) << std::endl; 
			  area = area + 0.636*mign.fanout_count(v); 
			  for (auto x = 0; x <c.size(); ++x)
			  {
				  if (c[x].complemented == 1)
					   area = area + (0.636/2)*mign.fanout_count(v); 
			  }
		  } 
   	   if (mign.level(v) == 2)
   		  {
			  std::cout << " level 2 con fanout " << mign.fanout_count(v)<< std::endl; 
   			  area = area + 4.66*mign.fanout_count(v); 
			  for (auto x = 0; x <c.size(); ++x)
			  {
				  if (c[x].complemented == 1)
					   area = area + (4.66/2)*mign.fanout_count(v);  
			  }
   		  } 
		  
		  if (mign.level(v) == 3)
		  {
				  std::cout << " level 3" << mign.fanout_count(v) << std::endl; 
      			  area = area + 38.24*mign.fanout_count(v); 
				  for (auto x = 0; x <c.size(); ++x)
				  {
					  if (c[x].complemented == 1)
						   area = area + (38.24/2)*mign.fanout_count(v);  
				  }
      	   }
		  }
	 else 
	 {
  	   if (mign.level(v) == 1)
  		  {
  			  std::cout << " level 1" << mign.fanout_count(v) << std::endl; 
  			  area = area + 0.636; 
  			  for (auto x = 0; x <c.size(); ++x)
  			  {
  				  if (c[x].complemented == 1)
  					   area = area + (0.636/2); 
  			  }
  		  } 
     	   if (mign.level(v) == 2)
     		  {
  			  std::cout << " level 2 con fanout " << mign.fanout_count(v)<< std::endl; 
     			  area = area + 4.66; 
  			  for (auto x = 0; x <c.size(); ++x)
  			  {
  				  if (c[x].complemented == 1)
  					   area = area + (4.66/2);  
  			  }
     		  } 
		  
  		  if (mign.level(v) == 3)
  		  {
  				  std::cout << " level 3" << mign.fanout_count(v) << std::endl; 
        			  area = area + 38.24; 
  				  for (auto x = 0; x <c.size(); ++x)
  				  {
  					  if (c[x].complemented == 1)
  						   area = area + (38.24/2);  
  				  }
        	   }
  		  }
	 }
	 
	 std::cout << "[i] Delay " << std::endl; 
	 
			  
		  
     return area; 
}
}


// Local Variables:
// c-basic-offset: 2
// eval: (c-set-offset 'substatement-open 0)
// eval: (c-set-offset 'innamespace 0)
// End:
