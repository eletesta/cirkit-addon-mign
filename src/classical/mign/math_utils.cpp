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

#include "math_utils.hpp"

#include <fstream>
#include <unordered_map>
#include <vector>



namespace cirkit
{
	// calculate n! 
double fact (signed n)
{ 
		auto fact = 1u; 
		if (n <= 1)
			return fact; 
		while ( n > 1)
		{
			fact *= n; 
			n = n-1; 
		}
		return fact; 
}

unsigned int int_pow2 (unsigned int x)
{
	auto ans = 1u; 
	
	while (x > 0){
		ans = ans*2 ;
		x = x -1; 
	}
	return ans; 
	
}
}
// Local Variables:
// c-basic-offset: 2
// eval: (c-set-offset 'substatement-open 0)
// eval: (c-set-offset 'innamespace 0)
// End:
