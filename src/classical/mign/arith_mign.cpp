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

#include "arith_mign.hpp"

#include <fstream>
#include <iostream>

#include <boost/optional/optional_io.hpp>
#include <boost/format.hpp>
#include <boost/graph/graphviz.hpp>

#include <core/utils/string_utils.hpp>
#include <classical/mig/mig.hpp> 
#include <classical/mig/mig_utils.hpp>
#include <classical/mign/math_utils.hpp>

namespace cirkit 
{
class store 
	{
	public:
		mign_graph mign; 
		std::unordered_map<std::string, mign_function> name_to_function; 
		
	};
	
class store_final
		{
		public:
			mign_graph mign; 
			std::unordered_map<std::string, mign_function> name_to_function; 
			std::vector<std::vector<std::string>> outnames; 
		
		};
	// function to compute the log2 of an unsigned number 
unsigned int int_log2 (unsigned int x)
{
	auto ans = 0u; 
	
	while (x > 0){
		x = x/2; 
		ans++; 
	}
	return ans; 
	
}

//compute the power of 2



store write_bit_count (unsigned Nin, std::vector<mign_function> operands, std::vector<std::string> outnames, unsigned int id, mign_graph mign,std::unordered_map<std::string, mign_function> name_to_function, unsigned w)
{
	
	store results; 
	//std::cout << Nin << std::endl; 
	
	auto start = 0u; 
	signed int a, b/*, count*/; 
	const auto ilog = int_log2(Nin); 
	////std::cout << ilog << std::endl; 
	auto i = 0; 
	
	std::vector<std::string> wire_names; 
	
	unsigned mat[Nin + 1][ilog]; 
	unsigned int vect[ilog]; 
	unsigned int vg[Nin+1]; 
	unsigned int vl[Nin+1]; 
	
	for ( i = 0; i < Nin + 1; ++ i)
	{
		vg[i] = vl[i]= 0; 
		for ( auto j = 0; j <ilog; ++j)
		{
			mat[i][j] = 0; 
		}
		auto k = i; 
		auto q = 0; 
		while (k > 0)
		{
			mat[i][q] = k%2; 
			k = k/2; 
			++q; 
		}
		
	}
	
	vg[Nin] = 0; 
	vl[Nin] = 0; 
		
	// copied from MAC by Luca Amaru 
	
	for(auto j=0;j<ilog;++j) {
		
		start=0;
		vect[j]=0;
		//auto i = 0; 
		for(i=0;i<Nin+1;++i) {
			if(mat[i][j]==1) {
				if(start==0) {
					a=i;
					start=1;
				}
			}
			else {
				b=i-1;
				if(start==1) {
					vg[a]=1;
					vl[b]=1;
					//auto name = boost::str( boost::format( "ge%dle%d" ) % a % b ); 
					////std::cout << " aggiunto " << name << std::endl; 
					//wire_names.push_back(name); 						
				}
				start=0;
				
			}
		}
		
		if(start==1) {
			b=i-1;
			//auto name = boost::str( boost::format( "ge%dle%d" ) % a %b ); 
			////std::cout << " aggiunto " << name << std::endl; 
			//wire_names.push_back(name); 	
			vg[a]=1;
			vl[b]=1;					
		}
	}
			
	
		for ( i = 1; i <Nin + 1; ++i)
		{
			if (vg[i] == 1)
			{
				auto name = boost::str( boost::format( "ge%d_%d_%d" ) % i %id %w ); 
				//std::cout << " aggiunto " << name << std::endl; 
				wire_names.push_back(name);
				std::vector<unsigned> weights; 
				for (auto x = 0; x<operands.size(); ++x)
				{
					weights.push_back(1u);
				}
				name_to_function.insert ({name, mign.create_threshold(operands, i,0,weights)});
			}
		}
		
		for ( i = 1; i <Nin + 1; ++i)
		{
			if (vl[i] == 1)
			{
				auto name = boost::str( boost::format( "le%d_%d_%d" ) % i %id %w); 
				//std::cout << " aggiunto " << name << std::endl; 
				wire_names.push_back(name);
				std::vector<unsigned> weights; 
				for (auto x = 0; x<operands.size(); ++x)
				{
					weights.push_back(1u);
				}
				name_to_function.insert ({name, mign.create_threshold(operands, i,1,weights)});
			}
		}
		
		for ( auto j = 0; j <ilog; ++j)
		{
			//std::cout << " qui entri oppure no? " << std::endl; 
			start = 0; 
			for ( i =0; i<(Nin + 1); ++i)
			{
				if (mat[i][j] == 1)
				{
					if (start == 0)
					{
						a = i; 
						start = 1; 
					}
				}
				else {
					b = i-1; 
					if (start == 1)
					{
						auto name = boost::str( boost::format( "ge%dle%d_%d_%d" ) % a %b %id %w);
						//auto name = outnames[j];  
						//std::cout << " aggiunto " << name << std::endl; 
						auto sm0  = boost::str( boost::format( "ge%d_%d_%d" ) % a %id %w );
						auto sm1  =  boost::str( boost::format( "le%d_%d_%d" ) %b %id %w);
						
						//std::cout << " con figli " << sm0 << " " << sm1 << std::endl; 
						std::vector<mign_function> operands_small; 
						operands_small.push_back (find_mign(name_to_function,sm0)); 
						operands_small.push_back (find_mign(name_to_function,sm1)); 
						operands_small.push_back(mign.get_constant( false )); 
						name_to_function.insert({name, mign.create_maj(operands_small)}); 
						////std::cout << " inserito = " << name << std::endl; 
						++vect[j]; 
					}
					start = 0; 
					 
				}
			}
			if (start == 1)
			{
				b = i -1; 
				auto name = boost::str( boost::format( "ge%dle%d_%d_%d" ) % a %b %id %w );
				//auto name = outnames[j];  
				//std::cout << " aggiunto " << name << std::endl; 
				auto sm0  = boost::str( boost::format( "ge%d_%d_%d" ) % a %id %w );
				auto sm1  =  boost::str( boost::format( "le%d_%d_%d" ) %b %id %w );
				//std::cout << " con figli " << sm0 << " " << sm1 << std::endl; 
				
				std::vector<mign_function> operands_small; 
				operands_small.push_back (find_mign(name_to_function,sm0)); 
				operands_small.push_back (find_mign(name_to_function,sm1)); 
				operands_small.push_back(mign.get_constant( false )); 
				name_to_function.insert({name, mign.create_maj(operands_small)}); 
				////std::cout << " inserito = " << name << std::endl; 
				++vect[j];
			}
		}
		
		for ( auto j = 0; j <ilog; ++j)
		{
			start = 0; 
			if (vect[j] != 1)
			{
				std::vector<mign_function> operands_small;  
				for ( i = 0; i <Nin + 1; ++i)
				{
					if (mat[i][j] == 1)
					{
						if (start == 0)
						{
							a = i; 
							start = 1; 
						}
					}
					else {
						b = i -1; 
						if (start == 1)
						{
							auto name = boost::str( boost::format( "ge%dle%d_%d_%d" ) % a %b % id %w);
							//auto name = outnames[j]; 
							//std::cout << " aggiunto " << name << std::endl; 
							//std::cout << " cerco = " << name << std::endl; 
							operands_small.push_back (find_mign(name_to_function,name)); 
		
						}
						start = 0; 
					}
				}
				////std::cout <<"flag 6" << std::endl; 
				if (start == 1)
				{
					b = i-1; 
					auto name = boost::str( boost::format( "ge%dle%d_%d_%d" ) % a %b % id %w);
					//auto name = outnames[j]; 
					//std::cout << " aggiunto " << name << std::endl; 
					//std::cout << " cerco = " << name << std::endl; 
					operands_small.push_back (find_mign(name_to_function,name)); 	
				}
				
				for (i =0; i<vect[j]-2; ++i)
				{
					operands_small.push_back(mign.get_constant(true)); 
				}
				auto name = boost::str( boost::format( "fnl_or_%d_%d_%d" ) % j % id %w);
				//auto name = outnames[j]; 
				//std::cout << " aggiunto " << name << std::endl; 
				//++j; 
				wire_names.push_back(name);
				////std::cout << " Penso che il problerma sia qui" <<std::endl; 
				name_to_function.insert({name, mign.create_or(operands_small)}); 
				////std::cout << " Penso che il problerma sia qui" <<std::endl; 
			}
			
		}
		
		
		for ( auto j = 0; j<ilog; ++j)
		{
			start = 0; 
			if (vect[j] == 1) {
				for ( i = 0; i < Nin + 1; ++i)
				{
					if (mat[i][j] == 1)
					{
						if (start == 0)
						{
							a = i; 
							start = 1; 
						}
					}
					else {
						b = i-1; 
						if (start ==1)
						{
							auto name = boost::str( boost::format( "ge%dle%d_%d_%d" ) % a %b %id %w);
							auto name_out = outnames[j]; 
							name_to_function.insert({name_out, name_to_function[name]}); 
							//std::cout << "name out = " << name_out <<std::endl;
							//std::cout << " name" << name << std::endl;
							//mign.create_po(name_to_function[name],name_out ); 
							
						}
						start = 0;
					}
				}
				if (start == 1)
				{
					b = i -1; 
					auto name = boost::str( boost::format( "ge%dle%d_%d_%d" ) % a %b % id %w);
					auto name_out = outnames[j]; 
					//mign.create_po(name_to_function[name],name_out );
					name_to_function.insert({name_out, name_to_function[name]}); 
					//std::cout << "name out = " << name_out <<std::endl;
					//std::cout << " name" << name << std::endl;
					
				}
			}
			else 
			{
				auto name = boost::str( boost::format( "fnl_or_%d_%d_%d" ) % j %id %w);
				auto name_out = outnames[j]; 
				name_to_function.insert({name_out, name_to_function[name]}); 
				//std::cout << "name out = " << name_out <<std::endl;
				//std::cout << " name" << name << std::endl;
				//mign.create_po(name_to_function[name],name_out );
			}
		}
		
		results.name_to_function = name_to_function; 
		results.mign = mign; 
		return results; 
	
}

store write_add2 (mign_graph mign,std::unordered_map<std::string, mign_function> name_to_function, std::vector<std::vector<std::string>> names, unsigned opt)
{
	if (opt == 1)
	{
	store results; 
	auto input = names[0].size(); 
	//std::cout << input << std::endl; 
	std::vector<std::string> output_names;
	//std::vector<std::string> wire_names; 
	std::vector<mign_function> operands;
	
    name_to_function.insert( {"1'b0", mign.get_constant( false )} );
    name_to_function.insert( {"1'b1", mign.get_constant( true )} );
		
	for (auto i = 0; i <= input; ++i)	
	{
		auto name = boost::str( boost::format( "sum_%d" ) % i );
		output_names.push_back(name); 
	}
	
	////std::cout << " siamo in adder 2" << std::endl; 
	auto name = boost::str( boost::format( "MAJ0" ));
		
	operands.push_back (find_mign(name_to_function,names[0][0])); 
	operands.push_back (find_mign(name_to_function,names[1][0])); 
	operands.push_back (mign.get_constant( false )); // Cin
		
	for ( auto i = 1; i < input; ++i)
	{
		auto nand = boost::str( boost::format( "f_AND_%u" ) % i );
		auto nor = boost::str( boost::format( "f_OR_%u" ) % i );
		auto sm0  = names[0][i];  
		////std::cout << sm0 << std::endl; 
		auto sm1  = names[1][i];  
		////std::cout << sm1 << std::endl; 
		
		std::vector<mign_function> operand; 
		operand.push_back (find_mign(name_to_function,sm0)); 
		operand.push_back (find_mign(name_to_function,sm1)); 
		
		name_to_function.insert({nand, mign.create_and(operand)});
		name_to_function.insert({nor, mign.create_or(operand)});
	}
	
	name_to_function.insert({name, mign.create_maj(operands)});
	
	for (auto i = 1; i < input; ++i)	
	{
		std::vector<mign_function> operand_two; 
		for (auto j = 0; j < i; ++j)
		{
			std::vector<mign_function> operand; 
			for (auto k = j; k<= i; ++k)
			{
				if (k == j)
				{
					if (k == 0)
					{
						operand.push_back(find_mign(name_to_function,name)); 
					}
					else 
					{
						auto names = boost::str( boost::format( "f_AND_%u" ) % k );
						////std::cout << " qui entri " << std::endl; 
						operand.push_back(find_mign(name_to_function,names)); 
					}
				}
				else 
				{
					auto namet = boost::str( boost::format( "f_OR_%u" ) % k );
					operand.push_back(find_mign(name_to_function,namet)); 
				}
			}
			for (auto k = (i -j + 1); k <= (2*(i-j+1)-2); ++k)
			{
				operand.push_back(mign.get_constant(false));
			}
			auto namef = boost::str( boost::format( "f_%u_%u" ) % i % j );
			name_to_function.insert({namef, mign.create_maj(operand)});
			operand_two.push_back(find_mign(name_to_function,namef));
		}
		
		auto namee = boost::str( boost::format( "f_AND_%u" ) % i);
		operand_two.push_back(find_mign(name_to_function,namee));
		for ( auto k = i + 1; k <= (2*(i + 1)-2); ++k)
		{
			operand_two.push_back(mign.get_constant(true));
		}
		auto cout = boost::str( boost::format( "c_%u" ) % i);
		name_to_function.insert({cout, mign.create_maj(operand_two)});
	}
	
	for ( auto i = 0; i < input; ++i)
	{
		if ( i == 0) {
			operands.push_back (!find_mign(name_to_function,name)); 
			operands.push_back (!find_mign(name_to_function,name)); 
			name_to_function.insert({"out_0", mign.create_maj(operands)});
			auto name_out = boost::str( boost::format( "sum_%d" ) % i ); 
			mign.create_po(name_to_function["out_0"],name_out );
		}
	//	else if ( i == input)
		else {
			
			std::vector<mign_function> operand; 
			auto name_two = boost::str( boost::format( "aux_%d" ) % 0 );
			if ( i == 1)
			{
				auto sm0  = names[0][i];
				operand.push_back (find_mign(name_to_function,sm0));
				auto sm1  = names[1][i];
				operand.push_back (find_mign(name_to_function,sm1));
				operand.push_back (find_mign(name_to_function,name));
				name_two = boost::str( boost::format( "aux_%d" ) % i );
				name_to_function.insert({name_two, mign.create_maj(operand)});
				operand.push_back (!find_mign(name_to_function,name_two));
				operand.push_back (!find_mign(name_to_function,name_two));
				auto out = boost::str( boost::format( "out_%d" ) % i );
				name_to_function.insert({out, mign.create_maj(operand)});
				auto name_out = boost::str( boost::format( "sum_%d" ) % i ); 
				mign.create_po(name_to_function[out],name_out );
			}
			else 
			{
				auto sm0  = names[0][i];
				operand.push_back (find_mign(name_to_function,sm0));
				auto sm1  = names[1][i];
				operand.push_back (find_mign(name_to_function,sm1));
				auto sm2  = boost::str( boost::format( "c_%f" ) % (i-1) );
				operand.push_back (find_mign(name_to_function,sm2)); 
				name_two = boost::str( boost::format( "aux_%d" ) % i );
				name_to_function.insert({name_two, mign.create_maj(operand)});
				operand.push_back (!find_mign(name_to_function,name_two));
				operand.push_back (!find_mign(name_to_function,name_two));
				auto out = boost::str( boost::format( "out_%d" ) % i );
				name_to_function.insert({out, mign.create_maj(operand)});
				auto name_out = boost::str( boost::format( "sum_%d" ) % i ); 
				mign.create_po(name_to_function[out],name_out );
			}
		}
	}
	
	auto out  = boost::str( boost::format( "c_%d" ) % (input-1) );
	auto name_out = boost::str( boost::format( "sum_%d" ) % input ); 
	mign.create_po(name_to_function[out],name_out );
	
	results.name_to_function= name_to_function; 
	results.mign = mign; 
	return results; 
}
else 
{
	store results; 
	auto input = names[0].size(); 
	const auto num_output = 2*opt;  // opt mi dice anche quanti bit di output avro ;) 

	std::vector<std::string> output_names;
	
	std::vector<mign_function> operands;
	
    name_to_function.insert( {"1'b0", mign.get_constant( false )} );
    name_to_function.insert( {"1'b1", mign.get_constant( true )} );
		
	for (auto i = 0; i < num_output; ++i)	
	{
		auto name = boost::str( boost::format( "sum_%d" ) % i );
		output_names.push_back(name); 
	}
	
	////std::cout << " siamo in adder 2" << std::endl; 
	auto name = boost::str( boost::format( "MAJ0" ));
		
	operands.push_back (find_mign(name_to_function,names[0][0])); 
	operands.push_back (find_mign(name_to_function,names[1][0])); 
	operands.push_back (mign.get_constant( false )); // Cin
		
	for ( auto i = 1; i < input; ++i)
	{
		auto nand = boost::str( boost::format( "f_AND_%u" ) % i );
		auto nor = boost::str( boost::format( "f_OR_%u" ) % i );
		auto sm0  = names[0][i];  
		////std::cout << sm0 << std::endl; 
		auto sm1  = names[1][i];  
		////std::cout << sm1 << std::endl; 
		
		std::vector<mign_function> operand; 
		operand.push_back (find_mign(name_to_function,sm0)); 
		operand.push_back (find_mign(name_to_function,sm1)); 
		
		name_to_function.insert({nand, mign.create_and(operand)});
		name_to_function.insert({nor, mign.create_or(operand)});
	}
	
	name_to_function.insert({name, mign.create_maj(operands)});
	
	for (auto i = 1; i < input; ++i)	
	{
		std::vector<mign_function> operand_two; 
		for (auto j = 0; j < i; ++j)
		{
			std::vector<mign_function> operand; 
			for (auto k = j; k<= i; ++k)
			{
				if (k == j)
				{
					if (k == 0)
					{
						operand.push_back(find_mign(name_to_function,name)); 
					}
					else 
					{
						auto names = boost::str( boost::format( "f_AND_%u" ) % k );
						////std::cout << " qui entri " << std::endl; 
						operand.push_back(find_mign(name_to_function,names)); 
					}
				}
				else 
				{
					auto namet = boost::str( boost::format( "f_OR_%u" ) % k );
					operand.push_back(find_mign(name_to_function,namet)); 
				}
			}
			for (auto k = (i -j + 1); k <= (2*(i-j+1)-2); ++k)
			{
				operand.push_back(mign.get_constant(false));
			}
			auto namef = boost::str( boost::format( "f_%u_%u" ) % i % j );
			name_to_function.insert({namef, mign.create_maj(operand)});
			operand_two.push_back(find_mign(name_to_function,namef));
		}
		
		auto namee = boost::str( boost::format( "f_AND_%u" ) % i);
		operand_two.push_back(find_mign(name_to_function,namee));
		for ( auto k = i + 1; k <= (2*(i + 1)-2); ++k)
		{
			operand_two.push_back(mign.get_constant(true));
		}
		auto cout = boost::str( boost::format( "c_%u" ) % i);
		name_to_function.insert({cout, mign.create_maj(operand_two)});
	}
	
	for ( auto i = 0; i < num_output; ++i)
	{
		if ( i == 0) {
			operands.push_back (!find_mign(name_to_function,name)); 
			operands.push_back (!find_mign(name_to_function,name)); 
			name_to_function.insert({"out_0", mign.create_maj(operands)});
			auto name_out = boost::str( boost::format( "sum_%d" ) % i ); 
			if (name_to_function["out_0"] != 0)
			mign.create_po(name_to_function["out_0"],name_out );
		}
	//	else if ( i == input)
		else {
			
			std::vector<mign_function> operand; 
			auto name_two = boost::str( boost::format( "aux_%d" ) % 0 );
			if ( i == 1)
			{
				auto sm0  = names[0][i];
				operand.push_back (find_mign(name_to_function,sm0));
				auto sm1  = names[1][i];
				operand.push_back (find_mign(name_to_function,sm1));
				operand.push_back (find_mign(name_to_function,name));
				name_two = boost::str( boost::format( "aux_%d" ) % i );
				name_to_function.insert({name_two, mign.create_maj(operand)});
				operand.push_back (!find_mign(name_to_function,name_two));
				operand.push_back (!find_mign(name_to_function,name_two));
				auto out = boost::str( boost::format( "out_%d" ) % i );
				name_to_function.insert({out, mign.create_maj(operand)});
				auto name_out = boost::str( boost::format( "sum_%d" ) % i ); 
				if (name_to_function[out] != 0)
				mign.create_po(name_to_function[out],name_out );
			}
			else 
			{
				auto sm0  = names[0][i];
				operand.push_back (find_mign(name_to_function,sm0));
				auto sm1  = names[1][i];
				operand.push_back (find_mign(name_to_function,sm1));
				auto sm2  = boost::str( boost::format( "c_%f" ) % (i-1) );
				operand.push_back (find_mign(name_to_function,sm2)); 
				name_two = boost::str( boost::format( "aux_%d" ) % i );
				name_to_function.insert({name_two, mign.create_maj(operand)});
				operand.push_back (!find_mign(name_to_function,name_two));
				operand.push_back (!find_mign(name_to_function,name_two));
				auto out = boost::str( boost::format( "out_%d" ) % i );
				name_to_function.insert({out, mign.create_maj(operand)});
				auto name_out = boost::str( boost::format( "sum_%d" ) % i ); 
				if (name_to_function[out] != 0)
				mign.create_po(name_to_function[out],name_out );
			}
		}
	}
	

	results.name_to_function= name_to_function; 
	results.mign = mign; 
	return results; 
}
}

store_final write_addm (mign_graph mign,std::unordered_map<std::string, mign_function>& name_to_function, std::vector<std::vector<std::string>> names, unsigned w)
{
	//std::cout << " *********** " << std::endl; 
	auto m = names.size(); 
	//std::cout << m << std::endl; 
	auto n = names[0].size();
	//std::cout << n << std::endl;  
	auto mlog = int_log2(m); 
	//std::cout << mlog << std::endl; 
	auto new_w = mlog; 
	
	std::vector<mign_function> operands; 
	std::vector<std::vector<std::string>> new_names; 
	new_names.resize(new_w); 
	
	store_final final_results; 
	
	for (auto j = 0; j < n; ++j)
	{
		
		std::vector<mign_function> operands; 
	 for ( auto i = 0; i < m; ++i)
	 {
		 //std::cout << names[i][j] << std::endl; 
		operands.push_back(name_to_function[names[i][j]]);
	 }
	
	
	std::vector<std::string> outnames;
	
	for ( auto g=0; g <(mlog); ++g)
	{
		auto name = boost::str( boost::format( "bitc_%u_%u_%u" ) % j %g %w);
		//std::cout << "output " << name << std::endl; 
		outnames.push_back(name); 
	}
	
	////std::cout << " ciao1" << std::endl; 
	auto count = 0u; 
	auto results = write_bit_count (m,operands,outnames,j,mign, name_to_function,w); 
	mign = results.mign; 
	name_to_function = results.name_to_function; 
 
 	////std::cout << "prima dei nometti" << std::endl; 
	if (j == 0)
	{
		for ( auto k = 0; k < mlog; ++k)
		{
			if ( k == 0){
				new_names[k].push_back(outnames[count]); 
				++count; 
			}
			else if (k > 0)
			{
					for ( auto g = 0; g < k; ++g)
						{
									new_names[k].push_back("1'b0"); 
										}
								new_names[k].push_back(outnames[count]);
									++count; 
					}
		}
	}
	
 	else 
	{
		for ( auto k = 0; k < mlog; ++k)
		{
				new_names[k].push_back(outnames[count]); 
				++count; 
		}
	}
	

	
    }
	
	auto count = 0u; 
	for (auto & e :new_names)
	{
		if (e.size() < (mlog + n - 1))
		{
			for (auto g = 0; g< (mlog + n - e.size()); ++g)
			{
				new_names[count].push_back("1'b0"); 
			}
		}
		++count; 
	}
	

	//std::cout << " proprio alla fine" << std::endl; 
		final_results.outnames = new_names; 
		final_results.mign = mign; 
		final_results.name_to_function = name_to_function; 
		return final_results; 
}
// Public Functions 

mign_function find_mign(const std::unordered_map<std::string, mign_function>& name_to_function, const std::string& name )
{
	 
	return name_to_function.at( name );
}
    
mign_graph exact_count (unsigned Nin, unsigned count)
{
	mign_graph mign; 
	mign.set_name("top"); 
	//mign.set_name("top"); ; 
	std::unordered_map<std::string, mign_function> name_to_function;
	std::vector<std::string> output_names;
	//std::vector<std::string> wire_names; 
	std::vector<mign_function> operands;
	
    name_to_function.insert( {"1'b0", mign.get_constant( false )} );
    name_to_function.insert( {"1'b1", mign.get_constant( true )} );
	
	for ( auto i = 0; i < Nin; ++ i)
	{
		auto name = boost::str( boost::format( "pi_%d" ) % i ); 
		name_to_function.insert ({name, mign.create_pi(name)}); 
	}
	
		auto out_name = boost::str( boost::format( "count_e" )); 
		output_names.push_back(out_name); 
		
		auto nameg = boost::str( boost::format( "Gthan%d" ) % count); 
		//wire_names.push_back(nameg);
		auto namel = boost::str( boost::format( "Lthan%d" ) % count); 
		//wire_names.push_back(namel);
		
		for (auto& input : mign.inputs())
		{
			operands.push_back(input.first); 
		}
		
		std::vector<unsigned> weights; 
		for (auto x = 0; x<operands.size(); ++x)
		{
			weights.push_back(1u);
		}
		name_to_function.insert ({nameg, mign.create_threshold(operands, count,0,weights)}); // grater equal count - 1
		name_to_function.insert ({namel, mign.create_threshold(operands, count,1,weights)}); // less eqaul count + 1
		
		std::vector<mign_function> operands_or;
		operands_or.push_back(find_mign(name_to_function,nameg)); 
		operands_or.push_back(find_mign(name_to_function,namel)); 
		
		const auto last = "last"; 
		name_to_function.insert ({last, mign.create_or(operands_or)}); 
		
		mign.create_po(name_to_function[last],out_name );
		
		return mign; 
}

mign_graph bitcount (unsigned Nin)
{
	mign_graph mign; 
	auto start = 0u; 
	mign.set_name("top"); ; 
	signed int a, b/*, count*/; 
	const auto ilog = int_log2(Nin); 
	auto i = 0; 
	
	std::unordered_map<std::string, mign_function> name_to_function;
    std::vector<std::string> output_names;
	std::vector<std::string> wire_names; 
	
	unsigned mat[Nin + 1][ilog]; 
	unsigned int vect[ilog]; 
	unsigned int vg[Nin+1]; 
	unsigned int vl[Nin+1]; 
	
	for ( i = 0; i < Nin + 1; ++ i)
	{
		vg[i] = vl[i]= 0; 
		for ( auto j = 0; j <ilog; ++j)
		{
			mat[i][j] = 0; 
		}
		auto k = i; 
		auto q = 0; 
		while (k > 0)
		{
			mat[i][q] = k%2; 
			k = k/2; 
			++q; 
		}
		
		////std::cout << " i = " << i << std::endl; 
		//for ( auto j = 0; j <ilog; ++j)
		//{
		//	//std::cout <<"mat[i][j] =" << mat[i][j]<< std::endl; 
		//}
	}
			
    name_to_function.insert( {"1'b0", mign.get_constant( false )} );
    name_to_function.insert( {"1'b1", mign.get_constant( true )} );
	
	for ( i = 0; i < Nin; ++ i)
	{
		auto name = boost::str( boost::format( "pi_%d" ) % i ); 
		name_to_function.insert ({name, mign.create_pi(name)}); 
	}
	
	for ( i = 0; i < ilog; ++i )
	{
		auto name = boost::str( boost::format( "bitcount_%d" ) % i ); 
		output_names.push_back(name); 
	}
		
	
	// copied from MAC by Luca Amaru 
	
	for(auto j=0;j<ilog;++j) {
		
		start=0;
		vect[j]=0;
		//auto i = 0; 
		for(i=0;i<Nin+1;++i) {
			if(mat[i][j]==1) {
				if(start==0) {
					a=i;
					start=1;
				}
			}
			else {
				b=i-1;
				if(start==1) {
					vg[a]=1;
					vl[b]=1;
					auto name = boost::str( boost::format( "ge%dle%d" ) % a % b ); 
					wire_names.push_back(name); 						
				}
				start=0;
				
			}
		}
		
		if(start==1) {
			b=i-1;
			auto name = boost::str( boost::format( "ge%dle%d" ) % a %b ); 
			wire_names.push_back(name); 	
			vg[a]=1;
			vl[b]=1;					
		}
	}
	
	
		std::vector<mign_function> operands;
		for ( const auto& input_node : mign.inputs())
		{
		operands.push_back(input_node.first);  
		} 
		

		
	for( i=0;i<Nin+1;++i) {
		if(vg[i]==1) {
			auto name = boost::str( boost::format( "ge%d" ) % i ); 
			wire_names.push_back(name); 
			std::vector<unsigned> weights; 
			for (auto x = 0; x<operands.size(); ++x)
			{
				weights.push_back(1u);
			}
			name_to_function.insert ({name, mign.create_threshold(operands, i,0,weights)});	
			//std::cout << " inserito = " << name << std::endl; 
		}
	}
	
	
	start=0;
	
	for(i=0;i<Nin+1;++i) {
		if(vl[i]==1) {
			auto name = boost::str( boost::format( "le%d" ) % i ); 
			wire_names.push_back(name);
			std::vector<unsigned> weights; 
			for (auto x = 0; x<operands.size(); ++x)
			{
				weights.push_back(1u);
			}
			name_to_function.insert ({name, mign.create_threshold(operands, i,1,weights)});
			//std::cout << " inserito = " << name << std::endl; 	
		}
			
		}
		
		for ( auto j = 0; j <ilog; ++j)
		{
			start = 0; 
			//if (vect[j] == 1){
			for ( i =0; i<Nin + 1; ++i)
			{
				if (mat[i][j] == 1)
				{
					if (start == 0)
					{
						a = i; 
						start = 1; 
					}
				}
				else {
					b = i-1; 
					if (start == 1)
					{
						auto name = boost::str( boost::format( "ge%dle%d" ) % a %b ); 
						auto sm0  = boost::str( boost::format( "ge%d" ) % a );
						auto sm1  =  boost::str( boost::format( "le%d" ) %b );
						
						std::vector<mign_function> operands; 
						operands.push_back (find_mign(name_to_function,sm0)); 
						operands.push_back (find_mign(name_to_function,sm1)); 
						operands.push_back(mign.get_constant( false )); 
						name_to_function.insert({name, mign.create_maj(operands)}); 
						//std::cout << " inserito = " << name << std::endl; 
						++vect[j]; 
					}
					start = 0; 
					 
				}
			}
			if (start == 1)
			{
				b = i -1; 
				auto name = boost::str( boost::format( "ge%dle%d" ) % a %b ); 
				auto sm0  = boost::str( boost::format( "ge%d" ) % a );
				auto sm1  =  boost::str( boost::format( "le%d" ) %b );
				
				std::vector<mign_function> operands; 
				operands.push_back (find_mign(name_to_function,sm0)); 
				operands.push_back (find_mign(name_to_function,sm1)); 
				operands.push_back(mign.get_constant( false )); 
				name_to_function.insert({name, mign.create_maj(operands)}); 
				//std::cout << " inserito = " << name << std::endl; 
				++vect[j];
			}
			//}
		}
		
		for ( auto j = 0; j <ilog; ++j)
		{
			start = 0; 
			if (vect[j] != 1)
			{
				std::vector<mign_function> operands;  
				for ( i = 0; i <Nin + 1; ++i)
				{
					if (mat[i][j] == 1)
					{
						if (start == 0)
						{
							a = i; 
							start = 1; 
						}
					}
					else {
						b = i -1; 
						if (start == 1)
						{
							auto name = boost::str( boost::format( "ge%dle%d" ) % a %b );
							//std::cout << " inserito = " << name << std::endl; 
							//std::cout << " cerco = " << name << std::endl; 
							operands.push_back (find_mign(name_to_function,name)); 
		
						}
						start = 0; 
					}
				}
				////std::cout <<"flag 6" << std::endl; 
				if (start == 1)
				{
					b = i-1; 
					auto name = boost::str( boost::format( "ge%dle%d" ) % a %b );
					//std::cout << " inserito = " << name << std::endl; 
					//std::cout << " cerco = " << name << std::endl; 
					operands.push_back (find_mign(name_to_function,name)); 	
				}
				
				for (i =0; i<vect[j]-2; ++i)
				{
					operands.push_back(mign.get_constant(true)); 
				}
				auto name = boost::str( boost::format( "fnl_or_%d" ) % j);
				//std::cout << " inserito = " << name << std::endl; 
				//++j; 
				wire_names.push_back(name);
				////std::cout << " Penso che il problerma sia qui" <<std::endl; 
				name_to_function.insert({name, mign.create_or(operands)}); 
				////std::cout << " Penso che il problerma sia qui" <<std::endl; 
			}
			
		}
		
		
		for ( auto j = 0; j<ilog; ++j)
		{
			start = 0; 
			if (vect[j] == 1) {
				for ( i = 0; i < Nin + 1; ++i)
				{
					if (mat[i][j] == 1)
					{
						if (start == 0)
						{
							a = i; 
							start = 1; 
						}
					}
					else {
						b = i-1; 
						if (start ==1)
						{
							auto name = boost::str( boost::format( "ge%dle%d" ) % a %b);
							auto name_out = boost::str( boost::format( "bitcount_%d" ) % j ); 
							//std::cout << "name out = " << name_out <<std::endl; 
							//std::cout << " name" << name << std::endl; 
							mign.create_po(name_to_function[name],name_out ); 
							
						}
						start = 0;
					}
				}
				if (start == 1)
				{
					b = i -1; 
					auto name = boost::str( boost::format( "ge%dle%d" ) % a %b);
					auto name_out = boost::str( boost::format( "bitcount_%d" ) % j ); 
					mign.create_po(name_to_function[name],name_out );
					//std::cout << "name out = " << name_out <<std::endl; 
					//std::cout << " name" << name << std::endl; 
					
				}
			}
			else 
			{
				auto name = boost::str( boost::format( "fnl_or_%d" ) % j);
				auto name_out = boost::str( boost::format( "bitcount_%d" ) % j ); 
				mign.create_po(name_to_function[name],name_out );
			}
		}
	
		////std::cout << " size = " << mign.size() << std::endl; 
		//for (auto output : mign.outputs())
		//{
			////std::cout << " siamo al nodo" << output.second << " e di output vediamo " << output.first.node<<std::endl;  
			//}
	    
	return mign; 
}

mign_graph add2_luca (unsigned input) // probabilmente la differenza sta nella depth 
{
	mign_graph mign; 
	mign.set_name("top"); ; 
	std::unordered_map<std::string, mign_function> name_to_function;
	std::vector<std::string> output_names;
	//std::vector<std::string> wire_names; 
	std::vector<mign_function> operands;
	
    name_to_function.insert( {"1'b0", mign.get_constant( false )} );
    name_to_function.insert( {"1'b1", mign.get_constant( true )} );
		
	for (auto i = 0; i <= input; ++i)	
	{
		auto name = boost::str( boost::format( "sum_%d" ) % i );
		output_names.push_back(name); 
	}
	
	auto name = boost::str( boost::format( "MAJ0" ));
	auto sma0  = boost::str( boost::format( "xi_%d" ) % 0 );
	name_to_function.insert ({sma0, mign.create_pi(sma0)});
    
    for (auto i = 1; i <input; ++i)
    {
        auto sm0  = boost::str( boost::format( "xi_%d" ) % i );
        name_to_function.insert ({sm0, mign.create_pi(sm0)});
    }
    
	auto sma1  = boost::str( boost::format( "yi_%d" ) % 0);
	name_to_function.insert ({sma1, mign.create_pi(sma1)});
    
    for (auto i = 1; i <input; ++i)
    {
        auto sm1  = boost::str( boost::format( "yi_%d" ) % i );
        name_to_function.insert ({sm1, mign.create_pi(sm1)});
    }
    
	auto cin = boost::str( boost::format( "cin" ));
	name_to_function.insert ({cin, mign.create_pi(cin)}); 
	
	operands.push_back (find_mign(name_to_function,sma0)); 
	operands.push_back (find_mign(name_to_function,sma1)); 
	operands.push_back (find_mign(name_to_function,cin));
		
	for ( auto i = 1; i < input; ++i)
	{
		auto nand = boost::str( boost::format( "f_AND_%u" ) % i );
		auto nor = boost::str( boost::format( "f_OR_%u" ) % i );
		auto sm0  = boost::str( boost::format( "xi_%d" ) % i ); 
		auto sm1  = boost::str( boost::format( "yi_%d" ) % i ); 
		
		std::vector<mign_function> operand; 
		operand.push_back (find_mign(name_to_function,sm0)); 
		operand.push_back (find_mign(name_to_function,sm1)); 
		
		name_to_function.insert({nand, mign.create_and(operand)});
		name_to_function.insert({nor, mign.create_or(operand)});
	}
	
	name_to_function.insert({name, mign.create_maj(operands)});
	
	for (auto i = 1; i < input; ++i)	
	{
		std::vector<mign_function> operand_two; 
		for (auto j = 0; j < i; ++j)
		{
			std::vector<mign_function> operand; 
			for (auto k = j; k<= i; ++k)
			{
				if (k == j)
				{
					if (k == 0)
					{
						operand.push_back(find_mign(name_to_function,name)); 
					}
					else 
					{
						auto names = boost::str( boost::format( "f_AND_%u" ) % k );
						////std::cout << " qui entri " << std::endl; 
						operand.push_back(find_mign(name_to_function,names)); 
					}
				}
				else 
				{
					auto namet = boost::str( boost::format( "f_OR_%u" ) % k );
					operand.push_back(find_mign(name_to_function,namet)); 
				}
			}
			for (auto k = (i -j + 1); k <= (2*(i-j+1)-2); ++k)
			{
				operand.push_back(mign.get_constant(false));
			}
			auto namef = boost::str( boost::format( "f_%u_%u" ) % i % j );
			name_to_function.insert({namef, mign.create_maj(operand)});
			operand_two.push_back(find_mign(name_to_function,namef));
		}
		
		auto namee = boost::str( boost::format( "f_AND_%u" ) % i);
		operand_two.push_back(find_mign(name_to_function,namee));
		for ( auto k = i + 1; k <= (2*(i + 1)-2); ++k)
		{
			operand_two.push_back(mign.get_constant(true));
		}
		auto cout = boost::str( boost::format( "c_%u" ) % i);
		name_to_function.insert({cout, mign.create_maj(operand_two)});
	}
	
	for ( auto i = 0; i < input; ++i)
	{
		if ( i == 0) {
			operands.push_back (!find_mign(name_to_function,name)); 
			operands.push_back (!find_mign(name_to_function,name)); 
			name_to_function.insert({"out_0", mign.create_maj(operands)});
			auto name_out = boost::str( boost::format( "sum_%d" ) % i ); 
			mign.create_po(name_to_function["out_0"],name_out );
		}
	//	else if ( i == input)
		else {
			
			std::vector<mign_function> operand; 
			auto name_two = boost::str( boost::format( "aux_%d" ) % 0 );
			if ( i == 1)
			{
				auto sm0  = boost::str( boost::format( "xi_%d" ) % i );
				operand.push_back (find_mign(name_to_function,sm0));
				auto sm1  = boost::str( boost::format( "yi_%d" ) % i );
				operand.push_back (find_mign(name_to_function,sm1));
				operand.push_back (find_mign(name_to_function,name));
				name_two = boost::str( boost::format( "aux_%d" ) % i );
				name_to_function.insert({name_two, mign.create_maj(operand)});
				operand.push_back (!find_mign(name_to_function,name_two));
				operand.push_back (!find_mign(name_to_function,name_two));
				auto out = boost::str( boost::format( "out_%d" ) % i );
				name_to_function.insert({out, mign.create_maj(operand)});
				auto name_out = boost::str( boost::format( "sum_%d" ) % i ); 
				mign.create_po(name_to_function[out],name_out );
			}
			else 
			{
				auto sm0  = boost::str( boost::format( "xi_%d" ) % i );
				operand.push_back (find_mign(name_to_function,sm0));
				auto sm1  = boost::str( boost::format( "yi_%d" ) % i );
				operand.push_back (find_mign(name_to_function,sm1));
				auto sm2  = boost::str( boost::format( "c_%f" ) % (i-1) );
				operand.push_back (find_mign(name_to_function,sm2)); 
				name_two = boost::str( boost::format( "aux_%d" ) % i );
				name_to_function.insert({name_two, mign.create_maj(operand)});
				operand.push_back (!find_mign(name_to_function,name_two));
				operand.push_back (!find_mign(name_to_function,name_two));
				auto out = boost::str( boost::format( "out_%d" ) % i );
				name_to_function.insert({out, mign.create_maj(operand)});
				auto name_out = boost::str( boost::format( "sum_%d" ) % i ); 
				mign.create_po(name_to_function[out],name_out );
			}
		}
	}
	
	auto out  = boost::str( boost::format( "c_%d" ) % (input-1) );
	auto name_out = boost::str( boost::format( "sum_%d" ) % input ); 
	mign.create_po(name_to_function[out],name_out );
	
	return mign; 
}

mign_graph threshold (unsigned Nin, unsigned threshold, unsigned polarity,std::vector<unsigned> weigths)
{
	mign_graph mign; 
	mign.set_name("top"); ; 
	std::unordered_map<std::string, mign_function> name_to_function;
	std::vector<std::string> output_names;
	//std::vector<std::string> wire_names; 
	std::vector<mign_function> operands;
	
    name_to_function.insert( {"1'b0", mign.get_constant( false )} );
    name_to_function.insert( {"1'b1", mign.get_constant( true )} );
	
	for ( auto i = 0 ; i < Nin; ++i )
	{
		auto name = boost::str( boost::format( "pi_%d" ) % i ); 
		name_to_function.insert ({name, mign.create_pi(name)}); 
		operands.push_back(find_mign(name_to_function,name)); 
	}
	
	auto name = boost::str( boost::format( "Threshold" )); 
	name_to_function.insert ({name, mign.create_threshold(operands,threshold,polarity,weigths)});
	
	auto name_out = boost::str( boost::format( "out" )); 
	mign.create_po(name_to_function[name],name_out );
	
	return mign; 
}

mign_graph addm (unsigned m, unsigned n)
{
	mign_graph mign; 
	mign.set_name("top"); ; 
	std::unordered_map<std::string, mign_function> name_to_function;

	std::vector<mign_function> operands;
	
	//auto mlog = int_log2(m); 
	//std::cout << " log " << mlog << std::endl; 
	
    name_to_function.insert( {"1'b0", mign.get_constant( false )} );
    name_to_function.insert( {"1'b1", mign.get_constant( true )} );
			
	std::vector<std::vector<std::string>> names;
	names.resize(m); 

	for (auto i = 0; i < m; ++i)
	{
		for (auto j = 0; j < n; ++j)
		{
			auto name = boost::str( boost::format( "x_%d_%d" ) % i %j);
			name_to_function.insert ({name, mign.create_pi(name)});
			names[i].push_back(name); 
		}
	}
	
	for ( auto w = m; w >= 2; w++)
	{
		if ( w != 2)
		{
			
			//std::cout << " w " << w << std::endl; 
			auto results = write_addm (mign, name_to_function,names,w); 
			mign = results.mign;
			name_to_function = results.name_to_function; 
			names.erase(names.begin(),names.end()); 
			names = results.outnames; 
			w=names.size(); 
			//std::cout << " w " << w << std::endl; 
		}
		else if (w == 2)
		{
			auto result = write_add2 (mign,name_to_function,names,1) ; 
			w = 1; 
			mign = result.mign; 
			name_to_function = result.name_to_function; 
		}
		w = w- 1; 
		
	}
	// w e oi qualcosa che mi dice le iterazioni ;)
		
	return mign; 

}



mign_graph mult2 (unsigned int Nin)
{
	
	//From "Optimal-Depth Threshold Circuits for Multiplication and Related Problems" by Yeh, Varvarigos, Parhami, Lee. 
	
	mign_graph mign; 
	mign.set_name("top"); ; 
	std::unordered_map<std::string, mign_function> name_to_function;
	std::vector<std::string> output_names, outnames;
	std::vector<std::string> wire_names; 
	
    name_to_function.insert( {"1'b0", mign.get_constant( false )} );
    name_to_function.insert( {"1'b1", mign.get_constant( true )} );
		
	for (auto i = 0; i < (2*Nin); ++i)	// quanto e lungo quello di uscita? 
	{
		auto name = boost::str( boost::format( "mult_%d" ) % i );
		output_names.push_back(name); 
	}
	
	for ( auto i =0;  i < Nin; ++i)
	{
		auto name = boost::str( boost::format( "x_%d" ) % i);
		name_to_function.insert ({name, mign.create_pi(name)});
	}
	
	for ( auto i =0;  i < Nin; ++i)
	{
		auto name = boost::str( boost::format( "y_%d" ) % i);
		name_to_function.insert ({name, mign.create_pi(name)});
	}
	
	for (auto i = 0; i < Nin; ++i)
	{
		for (auto j = 0; j < Nin; ++j)
		{
			std::vector<mign_function> operands; 
			auto sm1 = boost::str( boost::format( "x_%d" ) % i );
			operands.push_back(name_to_function[sm1]); 
			auto sm2 = boost::str( boost::format( "y_%d" ) % j );
			operands.push_back(name_to_function[sm2]); 
			auto name = boost::str( boost::format( "p_%d_%d" ) % j %(i+j));
			//std::cout << " inserito = " << name << std::endl; 
			name_to_function.insert ({name, mign.create_and(operands)});
		}
	}
	
	std::vector<mign_function> operands; 
	std::vector<std::vector<std::string>> names; // una matrice di stringhe 
	names.resize(Nin); 
	
	for (auto j = 0; j < Nin; ++j)
	{
		//auto i = 0u; 
		for ( auto i = 0; i<j;++i)
		{
			names[j].push_back("1'b0"); 
		}
		for ( auto k = 0; k < Nin; ++k)
		{
			names[j].push_back(boost::str( boost::format( "p_%d_%d" ) % j %(k+j)));
		}
		for (auto h = (j+Nin); h < (2*Nin-1); ++h)
		{
			names[j].push_back("1'b0"); 
		}
	}
	
	
	auto n = 2*Nin -1; 
	//std::cout << " n = " << n << std::endl; 
	auto m = Nin; 
	//std::cout << " m = " << m << std::endl; 
	//auto mlog = int_log2(m); 
	//std::cout << " nmlog = " << mlog << std::endl; 
	
	for ( auto w = m; w >= 2; w++)
	{
		if ( w != 2)
		{
			//std::cout << " w prima del write m" << std::endl; 
			auto results = write_addm (mign, name_to_function,names,w); 
			mign = results.mign;
			name_to_function = results.name_to_function; 
			names.erase(names.begin(),names.end()); 
			names = results.outnames; 
		
			w=names.size(); 
		}
		else if (w == 2)
		{
			auto result = write_add2 (mign,name_to_function,names,Nin) ; 
			w = 1; 
			mign = result.mign; 
			name_to_function = result.name_to_function; 
		}
		w = w- 1; 
		
	}
	return mign; 	
}
}
