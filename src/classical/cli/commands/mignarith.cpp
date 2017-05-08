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

#include "mignarith.hpp"

#include <boost/program_options.hpp>

#include <classical/cli/stores_mign.hpp>
#include <classical/mign/arith_mign.hpp>
#include <classical/mign/mign_rewrite.hpp>

using namespace boost::program_options;

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

mignarith_command::mignarith_command( const environment::ptr& env )
  : cirkit_command( env, "Create MIG-n for arithmetic components" )
{
  opts.add_options()
    ( "bitcount,b", value( &bits ),          "create bitcount circuit" )
	( "exact_counter,e", value(&numbers),    "create exact counter with count m of n input n:m")
    ( "add2,a", value( &bits),               "create add of two number of n bits, constant depth" )
	( "addm,am", value (&numbers),           "create add of m number of n bits m:n")
    ( "mult2,m", value (&bits),              "create mult of two numbers of n bits, constant depth")
    ( "threshold,t", value(&numbers),        "create threshold t gate of n number with polarity p n:t:p, weights can be inserted (1,1,1)")
    ;
}

bool mignarith_command::execute()
{
  auto& migns = env->store<mign_graph>();
  auto statistics = std::make_shared<properties>();
  auto settings = std::make_shared<properties>();

  if ( is_set( "bitcount" ) )
  {
    migns.extend();
    migns.current() = bitcount( bits);
  }
  if ( is_set( "exact_counter" ) )
  {
	  std::string str, str1;
	  
	 if (numbers[1] != ':')
	 {
	  str[0] = numbers[0]; 
	  str[1] = numbers[1]; 
	 }
	 else 
		 str[0] = numbers[0]; 
	 
	auto bits_exact = std::stoi(str);
	  
	  if (numbers[numbers.size()-2] != ':')
	  {
		  str1[0] = numbers[numbers.size() - 2];
		  str1[1] = numbers[numbers.size() - 1];
	  }
	  else 
		str1[0] = numbers[numbers.size() - 1];  
	  
	auto count = std::stoi(str1);
    migns.extend();
    migns.current() = exact_count( bits_exact ,count);
  }
  
  if ( is_set( "addm" ) )
  {
	  std::string str, str1;
	  
	 if (numbers[1] != ':')
	 {
	  str[0] = numbers[0]; 
	  str[1] = numbers[1]; 
	 }
	 else 
		 str[0] = numbers[0]; 
	 
	auto m = std::stoi(str);
	  
	  if (numbers[numbers.size()-2] != ':')
	  {
		  str1[0] = numbers[numbers.size() - 2];
		  str1[1] = numbers[numbers.size() - 1];
	  }
	  else 
		str1[0] = numbers[numbers.size() - 1];  
	  
	auto n = std::stoi(str1);
	auto mignadd = addm(m,n); 
    migns.extend();
    migns.current() = mign_rewrite_top_down (mignadd,settings,statistics);
  }
    
  if (is_set ("add2"))
  {
      migns.extend();
      migns.current() = add2_luca( bits);
  }
  
  if (is_set ("mult2"))
  {
	  auto mignmult= mult2( bits);
      migns.extend();
      migns.current() = mign_rewrite_top_down (mignmult,settings,statistics);
  }
     
  if ( is_set("threshold"))
  {
      std::string str, str1, str2, str3;
	  std::vector<unsigned> weights; 
	
	  auto flag = 1u; 
      
      if (numbers[1] != ':')
      {
          str[0] = numbers[0];
          str[1] = numbers[1];
		  
	      if (numbers[4] != ':')
	      {
	          str1[0] = numbers[3];
	          str1[1] = numbers[4];
			  str2[0] = numbers[6];
			  if (numbers.size() > 7)
			  {
				  flag = 0; 
				  for ( auto x = 0; x<(numbers.size()-8); ++x)
				  {
					  if (numbers[8+x] != ',')
					  {
						  str3[0] = numbers[8+x];
					  weights.push_back(std::stoi(str3)); 
				      }
					  ++x;
				  }
			  }
	      }
		  str1[0] = numbers[3];
		  str2[0] = numbers[5];
		  
		  if (numbers.size() > 6)
		  {
			  flag = 0; 
			  for ( auto x = 0; x<(numbers.size()-7); ++x)
			  {
				  if (numbers[7+x] != ',')
				  {
					  str3[0] = numbers[7+x];
				  weights.push_back(std::stoi(str3)); 
			      }
				  ++x;
			  }
		  }
      }
      else
	  {
          str[0] = numbers[0];
	      if (numbers[3] != ':')
	      {
	          str1[0] = numbers[2];
	          str1[1] = numbers[3];
			  str2[0] = numbers[5];
			  if (numbers.size() > 6)
			  {
				  flag = 0; 
				  for ( auto x = 0; x<(numbers.size()-7); ++x)
				  {
					  if (numbers[7+x] != ',')
					  {
						  str3[0] = numbers[7+x];
					  weights.push_back(std::stoi(str3)); 
				      }
					  ++x;
				  }
			  }
	      }
		  str1[0] = numbers[2];
		  str2[0] = numbers[4];
		  if (numbers.size() > 5)
		  {
			  flag = 0; 
			  for ( auto x = 0; x<(numbers.size()-6); ++x)
			  {
				  if (numbers[6+x] != ',')
				  {
					  str3[0] = numbers[6+x];
				  weights.push_back(std::stoi(str3)); 
			      }
				  ++x;
			  }
		  }
	  }
      
      auto input = std::stoi(str);
      auto thre = std::stoi(str1);
      auto polarity = std::stoi(str2);
      
      migns.extend();
	  migns.current() = threshold( input ,thre, polarity, weights);
  }

  return true;
}

}

// Local Variables:
// c-basic-offset: 2
// eval: (c-set-offset 'substatement-open 0)
// eval: (c-set-offset 'innamespace 0)
// End:
