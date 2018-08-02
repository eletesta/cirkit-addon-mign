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

#include "arithmetic_io.hpp"

#include <chrono>
#include <fstream>
#include <unordered_map>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <gmp.h>
#include <gmpxx.h>

#include <boost/algorithm/string/join.hpp>
#include <boost/format.hpp>

#include </usr/local/include/primesieve.hpp>

#include <classical/mign/math_utils.hpp>
#include <core/utils/string_utils.hpp>

using boost::format;

namespace cirkit
{

/******************************************************************************
 * Types                                                                      *
 ******************************************************************************/
mpz_class power (mpz_class g, unsigned j, mpz_class mod)
{
	mpz_class pow;
	pow.set_str(std::to_string(1),10); 
	for (auto i = 1; i <= j; i++)
	{
		mpz_class p_p, p_prov; 
		p_p.set_str(std::to_string(pow.get_ui() * g.get_ui()),10); 
		p_prov = p_p % mod; 
		pow.set_str(std::to_string(p_prov.get_ui()),10); 
	}
	return pow; 
}
/******************************************************************************
 * Private functions                                                          *
 ******************************************************************************/

void write_iterated ( std::ostream& os, const properties::ptr& settings = properties::ptr())
{
	auto time   = std::chrono::system_clock::now();
	auto time_c = std::chrono::system_clock::to_time_t( time );
	const auto header_prefix      = get( settings, "header_prefix",      std::string( "written by CirKit -- mign " ) );
	os << "// " << header_prefix << " " << std::ctime( &time_c ) << std::endl;
		
   const auto default_name  = get( settings, "default_name", std::string( "top" ) );
   const auto number_inputs = get( settings, "inputs",   4);
   const auto number_bits   = get( settings, "bits",  2);
   
   std::vector<unsigned> coefficients; 
   coefficients.push_back(pow(2,2) - 1); 
   unsigned i = 4;
   auto flag = 0; 
   auto prod_M = coefficients[0]; 
     
     auto max_number = pow(pow(2,number_bits) - 1, number_inputs); 
     std::cout << max_number <<std::endl; 
       
     std::vector<int> primes;

     primesieve::generate_primes(10000, &primes);
     auto flag_2 = 0; 
       
     while (flag == 0) 
     {	  
    	  flag_2 = 0; 
   	  bool valid = true; 
   	  unsigned j = 0; 
   	  //unsigned prov = (pow(2, i) - 1); 
      unsigned prov = (i); 
   	  auto it = find (primes.begin(), primes.end(), i);
   	  if (it == primes.end())
   		  flag_2 = 1;
   	  
   	  if (flag_2 == 0)
   	  {
   		/*for (auto & c : coefficients)
   		   	  {
   		   		  if (std::__gcd(c,prov) == 1)
   		   		  {
   		   			  j = j +1; 
   		   		  }	  
   		   		  else 
   		   		  {
   		   			  break; 
   		   		  }
   		   	  }*/
   		   	  //if (j == coefficients.size())
   		   		 // {
   		   		  	  coefficients.push_back(prov); 
   		   		  	  prod_M = prod_M * prov; 
   		   		//  }
   		   	  if ( prod_M > max_number) 
   		   		  flag = 1; 
   		   	  
   	  }
    	i = i + 1; 
     } 
     
     
     for (auto & c : coefficients)
     {
   	  std::cout << " m " << c << std::endl;
     }
     
   // auto max_bits = floor( log2(pow(coefficients[coefficients.size()-1] -1 ,number_inputs)) + 1 ); 
    auto max_bits = floor( log2(max_number) + 1 ); 
    std::cout << max_bits <<std::endl; 
     
     assert (prod_M > max_number); 
     
     std::vector<unsigned> gen; 
      for (auto i= 1; i <= coefficients.size(); i++)
      {
    	  flag = 1;
    	  //os << format ("m%d = m%d_fix;" ) % i % i << std::endl; 
    	 
    	  // Calcolo dei generators 
    	  for (auto h = 2; h < coefficients[i-1]; h++ )
    	  {
    		  if (flag == 0)
    		  		  break;
    		  unsigned int k = 0; 
    		  double gen_p = 0; 
    		  for (auto k = 1; k < coefficients[i-1] - 1; k++)
    		  {
    			  gen_p = (unsigned)pow(h,k) % coefficients[i-1]; 
    			  if (gen_p == 1)
    				 break;  
    		  }
    		  if (gen_p != 1)
    		  {
    			  gen.push_back(h); 
    			  flag = 0; 
    		  }
    	  }
    	 // os << format ("g%d = %d;" ) %i  %gen[i-1] << std::endl; 
      }
  
  /* sum module etc etc */
  
    os << format("module sum%d #(parameter n=1001, m=1001)") %coefficients.size() << std::endl; 
    os << "(" ; 
    for (auto y = 1; y <= coefficients.size(); y++)
    			os << format ("input [n-1:0]a_%d,") %y; 
    
	os << format("output [(n-1)*%d:0]res );") %coefficients.size() << std::endl; 
	os << "assign res = ";
	unsigned yy = 0; 
	for (yy = 1; yy <= coefficients.size() -1; yy++)
	    			os << format ("a_%d + ") %yy; 
	os << format ("a_%d ; ") %yy << std::endl; 
	os << "endmodule" << std::endl << std::endl;
	
	//os << "module mult #(parameter n=1001, m=1001)" << std::endl; 
	//os << "(input [n-1:0] a, input [m-1:0] b, output [m+n-1:0] mult);" << std::endl;
	//os << "assign mult = a * b; " << std::endl;
	//os << "endmodule" << std::endl << std::endl;
	
	//os << "module pow #(parameter n=33, out_pow = 128)" << std::endl;
	//os << "(a, b, pow);" << std::endl;
	//os << "input  [n-1:0] a;" << std::endl;
	//os << "input [31:0] b;" << std::endl;
	//os << "output reg [out_pow-1:0]pow;" << std::endl;
	//os << "integer i;" << std::endl; 
	//os << "always @(a,b) begin " << std::endl; 
	//os << "(for i = 0; i < b; i=i+1) begin" << std::endl; 
	//os << "pow = a**b; " << std::endl;
	//os << "endmodule " << std::endl << std::endl;
	
	// CONVERSION 
	
	os << "module conversionBIT #(parameter n=31, m = 31)" << std::endl;
		os << "(input [n-1:0] a, input [m-1:0] m1, output reg [m-1:0] crr1_a);" << std::endl;
		os << "wire [n-1:0]part_m1_a[n-1:0];" << std::endl;
		os << "wire  [n-1:0]prov; " << std::endl;
		os << "reg [n*2:0]crr1_a_p; " << std::endl;
		os << "reg [n*2:0]crr_1_sum[n-1:0]; " << std::endl;
		os << "integer j; " << std::endl;
		os << "generate" <<  std::endl;
		os << "genvar i;" <<  std::endl;
		os << "for (i=0; i < n;i=i+1) begin" << std::endl;
		os << "assign prov = 1; " << std::endl;
		os << "assign part_m1_a[i] = ((prov << (i)) * a[i]) % m1;" << std::endl;
		os << "end" <<  std::endl;
		os << "endgenerate" <<  std::endl;
		os << "always @(crr1_a_p,j,part_m1_a[0],a) begin" << std::endl;
		os << "crr1_a_p = 0; " << std::endl;
		os << "for (j=0; j < n;j=j+1) begin" << std::endl;
		os << "crr_1_sum[j] = crr1_a_p; " << std::endl;
		os << "crr1_a_p = part_m1_a[j] + crr_1_sum[j]; " << std::endl;
	    os << "crr1_a = crr1_a_p % m1; " << std::endl;
		os << "if (crr1_a_p < m1)" << std::endl;
		os << "crr1_a = crr1_a_p;" << std::endl;
		os << "end" <<  std::endl;
		os << "end" << std::endl;
		os << "endmodule" << std::endl << std::endl;
		
		
	for (auto t = 1; t <= coefficients.size(); t++)
	{
		os << format("module mux%d #( parameter n=33, m1 = 33)") %coefficients[t-1] << std::endl;
		os << "(" ;
		for (auto y = 1; y < coefficients[t-1]; y++)
			os << format ("input [n-1:0]lx_%d,") %y; 
		os << "input [n-1:0] a, output reg [n-1:0]result);" << std::endl;
		os << " integer k;" << std::endl; 
					
		os << "always @(a, "; 
				unsigned yt;
		for (yt = 1; yt < coefficients[t-1]-1; yt++)
		{
			os << format("lx_%d, ") %yt << std::endl;
		}
		os << format("lx_%d ) begin") %yt << std::endl;
		os << "result=0; " << std::endl;
		for (auto y = 1; y < coefficients[t-1]; y++)
		{
			if (y == 1)
			os << format("if ( a == lx_%d) begin") % y << std::endl;
			else
				os << format("else if ( a == lx_%d) begin") %y << std::endl;
			os << format("result = %d; ") % (y-1) << std::endl;
			os << "end" << std::endl;
		}
		os << "else " << std::endl;
		os << "result=0; " << std::endl;
		os << "end" << std::endl; 
		os << "endmodule" << std::endl << std::endl;
	}
		
	for (auto i = 0; i < coefficients.size(); i++)
	{
		os << format("module crr_mult%d_%d #(parameter m =33, i = 30)") %number_inputs % coefficients[i] << std::endl;
			os << "("; 
			for (auto y = 0; y < number_inputs; y++)
			os << format("input [m:0]result_1_%d, " ) % (y+1) ; 
			os <<	"input [m:0]m1,input [m:0]g1, output reg [m:0]abc_1);" << std::endl;
			os << format("wire [m*%d:0] sum_abc_1; ") %number_inputs << std::endl; 
			os << format("reg [m:0] gen_1[%d -1 :0]; " ) % ((coefficients[i]-1) * number_inputs ) << std::endl;
			//os << " wire [out_pow-1:0] gen_pow_1; " << std::endl;
			os << " integer sum_b; " << std::endl;

			os << format("sum%d #(m+1) inst57 ( ") %number_inputs; 
			for (auto y = 0; y < number_inputs; y++)
				os << format("result_1_%d, " ) % (y+1) ; 
			os <<		"sum_abc_1); " << std::endl;
			
			os << format("always@(sum_abc_1, sum_b, ") << std::endl; 
			for (auto k = 0; k < ((coefficients[i]-1) * number_inputs) -1; k++)
				os << format("gen_1[%d],") %k ; 
			os << format("gen_1[%d] ) begin ") % (((coefficients[i]-1) * number_inputs) -1) << std::endl; 
			os << "sum_b = sum_abc_1 ; " << std::endl;
			os << "abc_1 = 0; " << std::endl;
			for (auto k = 0; k < (coefficients[i]-1) * number_inputs; k++)
			{
				 mpz_class z; 
				 mpz_class g, genere; 
				 mpz_class mod; 
				 g.set_str(std::to_string(gen[i]),10); 
		         mod.set_str(std::to_string(coefficients[i]),10); 
				 auto p = power (g,k,mod); 
			     genere = p % mod; 
				 auto gj = genere.get_ui(); 
				 os << format ("gen_1[%d] = %d;") %(k) % gj << std::endl  ; 
			}
			for (auto y = 0; y < (coefficients[i]-1) * number_inputs; y++)
			{
				if (y == 0)
				os << format("if ( sum_b == %d) begin") % y << std::endl;
				else
				os << format("else if ( sum_b == %d) begin") %y << std::endl;
				os << format("abc_1 = gen_1[%d]; ") % y << std::endl;
						os << "end" << std::endl;
			}
			//os << "pow #(m+1) inst58 (g1, sum_b, gen_pow_1); " << std::endl;

			
			os << "end" << std::endl;
			os << "endmodule" << std::endl;
	}
	

   /* MAIN module top */
	
	os << "module top ( " ; 
  for (auto i = 0; i < number_inputs ; i++)
  {
	  os << format("input [%d-1:0] x_%d, " ) % number_bits % i ; 
  }
  os << format("output reg [%d-1:0] result_final);" ) % max_bits << std::endl; // (number_bits*number_inputs) << std::endl; 
  
  // We have number of bits as input. The highest number we can rpresent 
  // is 2^b -1. We have number_of_inputs inputs thus the maximum multiplication is 2^b - 1 ^ number_of_inputs
  
  // Scelta dei coefficienti per la base m1 m2 .. mk
  
  os << format ("parameter BIT = %d;" ) % max_bits << std::endl << std::endl; 
  std::string always; 
  
  std::vector<std::string> cij; 
  std::vector<unsigned> cij_v; 
    
  for (auto i= 1; i <= coefficients.size(); i++)
  {
	  os << format ("parameter m%d_fix = %d; " ) % i % coefficients[i-1]<< std::endl;
	  os << format ("reg  [BIT:0]m%d; " ) % i << std::endl; 
	  for (auto j = 1; j <= number_inputs; j++)
	  {
		  auto str = format ("crr_%d_%d, " ) % i %j;
		  always += str.str(); 
		  auto str1 = format ("result_%d_%d, " ) % i %j;
		  always += str1.str(); 
		  os << format ("wire [BIT:0]crr_%d_%d; " ) % i % j<< std::endl; // first number is the i which means the coefficinet, second is which input
		  os << format ("wire [BIT:0]result_%d_%d; ") % i % j << std::endl;  
	  }
	  os << format ("reg  [BIT:0]g%d; " ) % i << std::endl; 
	  os << format ("reg  [BIT:0]v%d; " ) % i << std::endl;
	  auto g = coefficients[i-1]-1;
	  os << format ("reg  [BIT:0]ge_i_mod_%d[%d:0]; " ) % i % g << std::endl;
	  if ( i != 1)
	  os << format ("reg  [BIT:0]prov_%d[%d:0]; " ) % i % (i) << std::endl;
	  
	  auto str = format ("v%d, " ) % i ;
	  always+=str.str(); 
	  os << format ("wire [BIT:0]abc_%d;  reg [BIT:0]abc_%d_mux; reg [BIT:0]x_base7_%d; ") % i % i % i << std::endl; 
	 // os << format ("reg [BIT:0]result_%d[%d-1:0];") % i %coefficients.size() << std::endl; 
	  auto str1 = format ("abc_%d, abc_%d_mux, x_base7_%d, " ) % i % i % i;
	  always+=str1.str(); 
	  
	  for (auto k = i+1; k <= coefficients.size(); k ++)
	  {
		  std::string str = "reg  [BIT:0]"; 
		  flag = 1; 
		  auto str1 = format (" c%d%d") % i % k; 
		  str += str1.str();
		  str += ";";
		  os << str << std::endl; 
		  cij.push_back( str1.str()); 
		  for (auto j = 1; j <= coefficients[k-1]; j++)
		  {
			  if (flag == 0) break; 
			  if (((j*coefficients[i-1]) % coefficients[k-1]) == 1)
			  {
				  flag = 0; 
				  cij_v.push_back(j); 
			  }
		  }
	  }
	  
	  os << std::endl; 
  }
  
  std::vector<unsigned> emmeprod(coefficients.size() -2); 
  std::vector<std::string> ms(coefficients.size() -2); 
  
  for (auto y = coefficients.size() -1; y > 1 ; y--)
  {
	  std::string str; 
	  emmeprod[coefficients.size() -1 - y] = 1; 
	  os << format ("reg  [BIT:0]" ); 
	  for (auto i = 1; i <=y; i++)
	  {
		  auto str1= format ("m%d" ) % i ;
		  str+=str1.str();
		  emmeprod[coefficients.size() -1 - y] = emmeprod[coefficients.size() -1 -y] * coefficients[i-1];
	  }
	  ms[coefficients.size() -1 -y]= str;
	  os << str << ";" << std::endl; 
  }
  
  unsigned int h=0; 
  for (h = 0; h < number_inputs-1; h++)
  {
	  auto str = format ("x_%d, " ) % h; 
	  always += str.str(); 
  }
  auto str = format ("x_%d" ) % h ; 
  always += str.str();  
  os << format ("always@ (%s) begin" ) % always << std::endl; 
  
 
  //std::vector<unsigned> gen; 
  for (auto i= 1; i <= coefficients.size(); i++)
  {
	  flag = 1;
	  os << format ("m%d = m%d_fix;" ) % i % i << std::endl; 
	 
	  // Calcolo dei generators 
	  /*for (auto h = 2; h < coefficients[i-1]; h++ )
	  {
		  if (flag == 0)
		  		  break;
		  unsigned int k = 0; 
		  double gen_p = 0; 
		  for (auto k = 1; k < coefficients[i-1] - 1; k++)
		  {
			  gen_p = (unsigned)pow(h,k) % coefficients[i-1]; 
			  if (gen_p == 1)
				 break;  
		  }
		  if (gen_p != 1)
		  {
			  gen.push_back(h); 
			  flag = 0; 
		  }
	  }*/
	  os << format ("g%d = %d;" ) %i  %gen[i-1] << std::endl; 
  }
  
  for (auto i= 0; i < coefficients.size()-2; i++)
  {
	  os << format ("%s = %d;" ) % ms[i] % emmeprod[i] << std::endl; 
  }
  for (auto i= 0; i < cij.size(); i++)
   {
 	  os << format ("%s = %d;" ) % cij[i] % cij_v[i] << std::endl; 
   }
  
  std::vector<std::vector<unsigned>> g_i_mod(coefficients.size()); 
  
  for (auto i = 0; i < coefficients.size(); i++)
  {
	  for (auto j = 0; j < coefficients[i] -1; j++ )
	  {  
		  mpz_class z; 
		  mpz_class g, genere; 
		  mpz_class mod; 
		 // auto g = mpz_class(unsigned);
		 // auto mod = mpz_class(unsigned);
		  g.set_str(std::to_string(gen[i]),10); 
		  //std::cout << g.get_ui() << std::endl; 
		  mod.set_str(std::to_string(coefficients[i]),10); 
		  auto p = power (g,j,mod); 
		  genere = p % mod; 
		  //mpz_set_ui(mod,coefficients[i]); 
		  //mpz_set_ui(z,pow(g,j)); 
		  //g_i_mod[i].push_back((unsigned)% coefficients[i]); 
		  //mpz_powm_ui (z, g, j, mod); 
		  //auto value = mpz_get_ui (z); 
		  g_i_mod[i].push_back(genere.get_ui()); 
		  os << format ("ge_i_mod_%d[%d] = %d;") %(i+1) % (j) % g_i_mod[i][j] << std::endl  ; 
	  }
  }
  
  os << std::endl; 
  for (auto i= 1; i <= coefficients.size(); i++)
  {
	  os << format ("if (" ); 
	  unsigned j = 1; 
	  
	  for (j = 1; j < number_inputs; j++)
	 {
		 os << format ("(crr_%d_%d == 0) || ") % (i) % (j); 	  
	  }
	  os << format ("(crr_%d_%d == 0))" )  % (i) % (j) << std::endl; 
	  os << format ("abc_%d_mux = 0;" ) % (i)  << std::endl; 
	  os << "else " << std::endl; 
	  os << format ("abc_%d_mux = abc_%d;" ) % (i)  % (i)<< std::endl; 
	  os << format ("x_base7_%d = abc_%d_mux;") % i % i<< std::endl; 
	  
	  for (j = 1; j <= number_inputs; j++)
	  {
	 // os << format ("result_%d[%d] = result_%d_%d;") % i % (j-1) %i %j<< std::endl;
	  }
  }
  
  os << std::endl; 
  
  os << "v1 = x_base7_1 % m1;" << std::endl; 
  
  for (auto i = 2; i <= coefficients.size(); i++)
  {
	  for (auto h = 0; h < i; h++)
	  {
		  if (h == 0)
			  os << format ("prov_%d[%d] = x_base7_%d; ") % i % h %i << std::endl; 
		  else if (h == 1)
		  { 
			  os << format ("if (prov_%d[%d] >= v%d)") % i % (h-1) % (h) << std::endl; 
			  os << format ("prov_%d[%d] = prov_%d[%d] - v%d;") % i % h % i % (h-1) % h << std::endl; 
			  os << " else" << std::endl; 
			  os << format ("prov_%d[%d] = prov_%d[%d] + m%d - v%d;") % i % h % i % (h-1) % i% h << std::endl; 
		  }
		  else 
		  {
			  os << format ("if (prov_%d[%d] * c%d%d >= v%d)") % i % (h-1) % (h-1) % i % (h) << std::endl; 
			  os << format ("prov_%d[%d] = ((prov_%d[%d] * c%d%d) - v%d)") % i % h % i % (h-1) % (h-1) % i % h << "%" << format( "m%d;") % i << std::endl; 
			  os << " else" << std::endl; 
			  os << format ("prov_%d[%d] = (prov_%d[%d] * c%d%d + m%d - v%d)") % i % h % i % (h-1) % (h-1) % i % i %h  << "%" << format( "m%d;") %  i<< std::endl; 
		  }
	  }
	  os << std::endl; 
	  os << format ("v%d = (prov_%d[%d] * c%d%d) " ) %i % i % (i-1) % (i-1) % i << " % " << format ("m%d;") %i << std::endl; 
  }
  
 os << "result_final = v1 + m1*v2 + " ; 
 
 unsigned tt = 0; 
 unsigned jj = 3; 
 for (jj = 3; jj < coefficients.size(); jj++ )
 {
	 os << format ("v%d*%s + ") % jj % ms[coefficients.size()-3-tt]; 
	 tt++;
 }
 os << format ("v%d*%s;") % jj % ms[coefficients.size()-3-tt] << std::endl << std::endl; 
 os << "end" << std::endl << std::endl; 
 
 auto counter = 0; 
 for (auto i = 0; i < coefficients.size(); i++)
 {
	 for (auto j = 0; j < number_inputs; j++)
	 {
		 os << format("conversionBIT #(%d,BIT+1) inst%d (x_%d,m%d,crr_%d_%d); ") % number_bits %counter % j % (i+1) % (i+1) % (j+1) << std::endl;
		 counter++; 
		 os << format("mux%d #(BIT+1,m%d_fix) inst%d  (") %coefficients[i] %(i+1) %counter ; 
		 for (auto y = 1; y <= coefficients[i] -1; y++)
			 os << format("ge_i_mod_%d[%d], " ) %(i+1) %(y-1); 
			os << format("crr_%d_%d, result_%d_%d); ")%(i+1) % (j+1) % (i+1) % (j+1) << std::endl;
		 counter++;
	 }
 }
 
 for (auto i = 1; i <= coefficients.size(); i++)
 {
	 os << format("crr_mult%d_%d #(BIT, %d) inst%d (") % number_inputs % coefficients[i-1] % number_inputs %counter; 
	 for (auto y = 1; y <= number_inputs; y++)
		os << format(	 "result_%d_%d, ") %i % y; 
			os << format("m%d, g%d, abc_%d); ")  %i %i %i<< std::endl;
	 counter++;
 }
 
 os << "endmodule" << std::endl; 
}

/******************************************************************************
 * Public functions                                                           *
 ******************************************************************************/

void write_iterated_mult_verilog ( const std::string& filename , const properties::ptr& settings)
{
  std::ofstream os( filename.c_str(), std::ofstream::out );
  write_iterated(  os, settings );
}

}
// Local Variables:
// c-basic-offset: 2
// eval: (c-set-offset 'substatement-open 0)
// eval: (c-set-offset 'innamespace 0)
// End:
