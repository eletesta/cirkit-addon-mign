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

#include "threshold_synthesis.hpp"

#include <iostream>
#include <string>
#include <vector>

#include <boost/assign/std/vector.hpp>
#include <boost/dynamic_bitset.hpp>
#include <classical/mign/arith_mign.hpp>

#include <core/utils/timer.hpp>
#include <classical/functions/isop.hpp>
#include <core/cube.hpp>

using namespace boost::assign;

namespace cirkit
{
	
bool is_const0(tt func)
{
	for (auto x = 0; x < func.size(); ++x)
	{
		if (func[x] == 0)
			{continue;}
		else 
			return false; 
	}
	
	return true; 
}

bool is_const1(tt func)
{
	for (auto x = 0; x < func.size(); ++x)
	{
		if (func[x] == 1)
			{continue;}
		else 
			return false; 
	}
	
	return true; 
}

signed evaluate_min_chow (std::vector<std::pair<unsigned,signed>> chow)
{
	signed min = chow[0].second; 
	
	for ( auto & e : chow)
	{
		if (e.second < min)
			min = e.second; 
	}
	return min; 
}
unsigned find_tempo(unsigned min_chow, std::vector<std::pair<unsigned, unsigned>> save)
{
	for (auto & saved : save)
	{
		if ( saved.first == min_chow)
			return saved.second; 
	}
	return 0; 
}

void vector_init (std::vector<unsigned> vettore)
{
	for (auto x = 0; x < vettore.size(); ++x)
	{
		vettore[x] = 0; 
	}
}

bool there_is (unsigned first, std::vector<std::pair<unsigned,std::pair<boost::dynamic_bitset<>,boost::dynamic_bitset<>>>> keys)
{
	auto flag = 1; 
	for (auto & x : keys)
	{
		if (x.first == first)
		{
			flag = 0; 
			//std::cout << "  e nella chiave " << std::endl; 
		}
	}
	return flag; 
	
}

void print_w (std::vector<unsigned> weigths)
{
	for ( auto & w :weigths)
	{
		std::cout << " weight : " << w << std::endl; 
	}
}

bool negative_unate (std::vector<unsigned> un, unsigned variable)
{
	if (un[variable] == 1)
	{
		//std::cout << " 	LA VARIABILE " << variable << " E NEGATIVE UNATE" << std::endl; 
		return 0; 
	}
		
	else return 1; 
}
/******************************************************************************
 * Private functions                                                          *
 ******************************************************************************/

class threshold_synthesis_manager
{
public:
  threshold_synthesis_manager( const tt& func,
                           const properties::ptr& settings,
                           const properties::ptr& statistics )
    : func( func ),
      statistics( statistics )
  {
    /* settings */
    verbose = get( settings, "verbose", false );
	very_verbose = get(settings, "very_verbose", false); 
  }

  mign_graph run()
  {
      properties_timer t( statistics );
  	  mign_graph mign; 
	
	  if ((is_const0(func)) || (is_const1(func)))
	  {
	  	std::cout << " The function is constant " << std::endl;
		return mign; 
	  }
	  
	std::vector<std::pair<unsigned,signed>> chow; 
	const auto num_vars = tt_num_vars( func );

	auto un = unate(func); 
	if (un.size() == num_vars)
	{
		//std::cout << " Unate " << std::endl; 
		chow = compute_chow(func,un); 
		mign = create_threshold(func,un,chow);
		 
	}
	else 
		{
			std::cout << " The function is not unate. Threshold gate cannot be done" << std::endl;
			return mign; 
		}
	 
	 	return mign; 
  }
  
  resulting run_T_w()
  {
    properties_timer t( statistics );
	std::vector<std::pair<unsigned,signed>> chow_tv; 
	const auto num_vars = tt_num_vars( func );
	
	std::vector<unsigned> weigths(1,0u); 
	std::vector<bool> neg_unate(1,0u); 
	
	resulting result; 
	
	 if ((is_const0(func)) || (is_const1(func)))
	{
		result.t_and_w = std::make_pair(0,weigths); 
		result.nega_un = neg_unate; 
		return result; 
	}
	
	auto un_tw = unate(func); 
	if (un_tw.size() == num_vars)
	{
		//std::cout << " Unate " << std::endl; 
		chow_tv = compute_chow(func,un_tw); 
		result = create_threshold_T_w(func,un_tw,chow_tv);
		 
	}
	else 
		{
			//std::cout << " The function is not unate. Threshold gate cannot be done" << std::endl;
			result.t_and_w = std::make_pair(0,weigths); 
			result.nega_un = neg_unate; 
		}
		
	return result; 
  }
  

private:
  std::vector<unsigned> unate(tt& table);
  std::vector<bool> nega_unate(tt& func); 
  std::vector<std::pair<unsigned,signed>> compute_chow( tt& func, std::vector<unsigned> un); 

  //flat_set<flat_set<unsigned>> find_gates_for_column( const unitized_table& table, unsigned column );
  //flat_set<unsigned> find_gate_for_table( const unitized_table& table );
  //flat_set<unsigned> find_gate_for_table_brute_force( const unitized_table& table );
  bool equation (boost::dynamic_bitset<> second_part_ineq, std::vector<unsigned> weigths, unsigned temporary);
  std::vector<std::pair<boost::dynamic_bitset<>,boost::dynamic_bitset<>>> create_ineq (tt& func, std::vector<unsigned> un); 
  std::vector<std::pair<boost::dynamic_bitset<>,boost::dynamic_bitset<>>> create_ALL_ineq (tt& func, std::vector<unsigned> un); 
  bool check_all (std::vector<unsigned> weigths, tt& func, std::vector<unsigned> un);
  unsigned compute_threshold (tt& func, std::vector<unsigned> weights, std::vector<unsigned> un);
  std::vector<std::pair<boost::dynamic_bitset<>,boost::dynamic_bitset<>>> select_ineq (tt& func, std::vector<std::pair<boost::dynamic_bitset<>,boost::dynamic_bitset<>>> ineq);
  mign_graph create_threshold( tt& func, std::vector<unsigned> pos, std::vector<std::pair<unsigned,signed>> chow);
  resulting create_threshold_T_w( tt& func, std::vector<unsigned> pos, std::vector<std::pair<unsigned,signed>> chow);
  std::vector<unsigned> weigths_assignment (std::vector<std::pair<boost::dynamic_bitset<>,boost::dynamic_bitset<>>> selected,std::vector<std::pair<unsigned,signed>> chow );
  std::vector<std::pair<unsigned,std::pair<boost::dynamic_bitset<>,boost::dynamic_bitset<>>>> create_keys (std::vector<std::pair<boost::dynamic_bitset<>,boost::dynamic_bitset<>>> selected); 
  

  

private:
  /* parameters */
  tt func;
 
  /* settings */
  bool verbose;
  bool very_verbose; 

  /* statistics */
  const properties::ptr& statistics;
};

std::vector<std::pair<boost::dynamic_bitset<>,boost::dynamic_bitset<>>> threshold_synthesis_manager::create_ALL_ineq (tt& func, std::vector<unsigned> un)
{
	const auto num_vars = tt_num_vars( func );
	std::vector<std::pair<boost::dynamic_bitset<>,boost::dynamic_bitset<>>> ineq; 
	
	std::vector< boost::dynamic_bitset<>> positive, negative; 
	
	
	for (auto x =0; x<func.size(); ++x)
	{
		
			boost::dynamic_bitset<> provv(num_vars), var(num_vars);
			for (auto j= 0; j <num_vars;++j)
			{
				auto y = j; //num_vars-j-1; 
				if (negative_unate(un,y) == 0)
				var = tt_nth_var(y).flip();
				else 
			    var =  tt_nth_var(y); 
				provv[j] = var[x]; 
			}
						 
			if (func[x] == 1)// greater part
				positive.push_back(provv); // bisogna vedere se funziona cosi 
			else 
				negative.push_back(provv); 
	}
	/*
	std::cout << " complete positive equation" << std::endl; 
	for (auto & h : positive)
		{
			std::cout << h << std::endl; 
		}
		std::cout << " complete negative equation" << std::endl; 
		for (auto & h : negative)
			{
				std::cout << h << std::endl; 
			}
		*/
	// create inequalities all together 
	for ( auto x = 0; x < positive.size(); ++x)
	{
		for ( auto y = 0; y < negative.size(); ++y)
		{
			ineq.push_back(std::make_pair(positive[x],negative[y])); 
		}
		 
	}
	
	return ineq; 
}
std::vector<std::pair<boost::dynamic_bitset<>,boost::dynamic_bitset<>>> threshold_synthesis_manager::create_ineq (tt& func, std::vector<unsigned> un)
{
	const auto num_vars = tt_num_vars( func );
	std::vector<std::pair<boost::dynamic_bitset<>,boost::dynamic_bitset<>>> ineq; 
	
	std::vector< boost::dynamic_bitset<>> positive, negative; 
	
	tt func_dc (func.size()); 
    std::vector<int> cover, cover_neg;
	
	std::vector<std::string> str, str_neg; 
    
    tt_isop( func, func, cover );
    const auto sop = cover_to_cubes( cover, num_vars );
	if (very_verbose)
	{
		std::cout << " func " << std::endl; 
	    common_pla_print( sop );
	}
	for ( const auto& c : sop )
	  {
	    str.push_back(c.to_string()); 
	  }

	
    tt_isop( ~func, ~func, cover_neg );
    const auto sop_neg = cover_to_cubes( cover_neg, num_vars );
	if (very_verbose)
	{
		std::cout << " ~func " << std::endl; 
	    common_pla_print( sop_neg );
	}
	
	for ( const auto& c : sop_neg )
	  {
	    str_neg.push_back(c.to_string()); 
	  }
		
	
	for (auto x =0; x<func.size(); ++x)
	{
		
			boost::dynamic_bitset<> provv(num_vars), var(num_vars);
			for (auto j= 0; j <num_vars;++j)
			{
				auto y = j; //num_vars - j -1; 
				if (negative_unate(un,y) == 0)
				var = tt_nth_var(y).flip(); 
				else
				var = tt_nth_var(y); 
				provv[j] = var[x]; 
			}
						 
			if (func[x] == 1)// greater part
				positive.push_back(provv); // bisogna vedere se funziona cosi 
			else 
				negative.push_back(provv); 
	}
	if (very_verbose)
	{
	std::cout << " positive " << std::endl; 
	for (auto & h : positive)
		{
			std::cout << h << std::endl; 
		}
		std::cout << " negative" << std::endl; 
		for (auto & h : negative)
			{
				std::cout << h << std::endl; 
			}
		}
	std::vector< boost::dynamic_bitset<>> new_positive, new_negative;
	
	for (auto& s: str)
	{
		//std::cout << " stiamo considerando " << s << std::endl; 
		std::vector<unsigned> salt;  
		for ( auto x = 0; x < s.size(); ++x)
		{
			if (s[x] == '-') // qui
			{
				//std::cout << " trattino alla variabile" << num_vars - x - 1 << std::endl;
				salt.push_back(x); 
			}
				
		}
		
		if (salt.size() != 0)
		{
			//std::cout << " quindi siamo nel caso salt size" << std::endl; 
		  for (auto & pos : positive)
		  {
			  auto count = 0u; 
			for ( auto j = 0; j <num_vars; ++j)
			{
				auto it = find (salt.begin(), salt.end(), j); // qui
				if (it != salt.end())
				    {
						if (pos[j] == 0)
							++count; 
					}
				else 
				
				{
					if (negative_unate(un,j) == 0)
					{
						if (!pos[j] == (unsigned)(s[j]-'0')) // qui
						++count; 
					}
					else 
					{
						if (pos[j] == (unsigned)(s[j]-'0'))
						++count;
					}
				}
			}
		
			if (count == num_vars)
			{
				//std::cout << " count = num_vars" << " per la variabile " << pos << std::endl; 
				new_positive.push_back(pos);
			}
		  }

	    }
		else 
		{
  		  for (auto & pos : positive)
  		  {
			  auto count = 0u; 
  			for ( auto j = 0; j <pos.size(); ++j)
  			{
 			   	if (negative_unate(un, j) == 0)
				{
					if (!pos[j] == (unsigned)(s[j]-'0')) // qui
					{
						++count;
					}
				}
				else 
				{
	  				if (pos[j] == (unsigned)(s[j]-'0')) // qui
	  				    {++count;}
				}
  				
  			}
  			if (count == pos.size())
  			new_positive.push_back(pos);
  		  }
		}
	}
	
	
	for (auto& s: str_neg)
	{
		std::vector<unsigned> salt;  
		for ( auto x = 0; x < s.size(); ++x)
		{
			if (s[x] == '-') // qui
			{
				salt.push_back(x);  // qui
			}
				
		}
		
		if (salt.size() != 0)
		{
		  for (auto & pos : negative)
		  {
			  auto count = 0u; 
			for ( auto j = 0; j <num_vars; ++j)
			{
				auto it = find (salt.begin(), salt.end(), (j));
				if (it != salt.end())
				    {
						 
						if (pos[j] == 1)
							++count; 
					}
				else 
				
				{
					if (negative_unate(un,j) == 0)
					{
						if (!pos[j] == (unsigned)(s[j]-'0'))
						++count; 
					}
					else 
					{
						if (pos[j] == (unsigned)(s[j]-'0'))
						++count;
					}
				}
			}
		
			if (count == num_vars)
			new_negative.push_back(pos);
		  }

	    }
		else 
		{
  		  for (auto & pos : negative)
  		  {
			  auto count = 0u; 
  			for ( auto j = 0; j <num_vars; ++j)
  			{
 
				if (negative_unate(un,j) == 0)
				{
					if (!pos[j] == (unsigned)(s[j]-'0'))
					++count; 
				}
				else 
				{
					if (pos[j] == (unsigned)(s[j]-'0'))
					++count;
				}
  			}
  			if (count == pos.size())
  			new_negative.push_back(pos);
  		  }
  		   
		}
	}
	


	// create inequalities all together 
	for ( auto x = 0; x < new_positive.size(); ++x)
	{
		for ( auto y = 0; y < new_negative.size(); ++y)
		{
			ineq.push_back(std::make_pair(new_positive[x],new_negative[y])); 
		}
		 
	}
	

	// further simplification 
	
	for (auto & inequa : ineq)
	{
		for (auto bit = 0; bit <num_vars; ++bit)
		{
			if ((inequa.first[bit] == 1) && (inequa.second[bit] == 1))
		{
			inequa.first.flip(bit); 
		    inequa.second.flip(bit);
		}
	    }
	}
	
	
	
	auto ineq_provv = ineq; 
	auto counting = 0u; 
	auto removed = 0u; 
	for (auto & inequa : ineq_provv)
	{
		
		auto count = 0u; 
		for (auto bit = 0; bit <num_vars; ++bit)
		{ if (inequa.second[bit] == 0)
		{
			++count; 
		}
		if (count == num_vars)
		{
			ineq.erase(ineq.begin() + counting - removed); 
			++removed; 
		}
	}
		++counting; 
			    
	}
	
	
	return ineq; 
}

std::vector<std::pair<boost::dynamic_bitset<>,boost::dynamic_bitset<>>> threshold_synthesis_manager::select_ineq (tt& func, std::vector<std::pair<boost::dynamic_bitset<>,boost::dynamic_bitset<>>> ineq)
{
	std::vector<std::pair<boost::dynamic_bitset<>,boost::dynamic_bitset<>>> selected_ineq; 
	const auto num_vars = tt_num_vars( func );
	
	auto ineq_provv = ineq; 
	
	for (auto & inequa : ineq)
	{
		auto count = 0u; 
		for (auto bit = 0; bit <num_vars; ++bit)
		{
		if (inequa.first[bit] == 1)
		{
			++count; 
		}
	    }
		if (count == 1)
		{
			selected_ineq.push_back(inequa); 
		}	
			    
	}
	
	
	return selected_ineq; 
}

std::vector<std::pair<unsigned,std::pair<boost::dynamic_bitset<>,boost::dynamic_bitset<>>>> threshold_synthesis_manager::create_keys (std::vector<std::pair<boost::dynamic_bitset<>,boost::dynamic_bitset<>>> selected)
{
	std::vector<std::pair<unsigned,std::pair<boost::dynamic_bitset<>,boost::dynamic_bitset<>>>> keys; 
		
	for ( auto & equa : selected)
	{
		for ( auto x = 0u; x < equa.first.size(); ++x)
		{
			if (equa.first[x] == 1)
			{
				//std::cout << " chiave " << std::endl; 
				keys.push_back(std::make_pair(x,equa)); //{x,equa}); 
				//std::cout << x << " = "<< equa.second << std::endl; 
			}
		}
	}
	
   
	
	return keys; 	
}



bool threshold_synthesis_manager::equation (boost::dynamic_bitset<> second_part_ineq, std::vector<unsigned> weigths, unsigned temporary)
{
	auto results = 0u; 
	for ( auto x = 0; x < weigths.size(); ++x)
	{
			results = results + weigths[x]*second_part_ineq[x]; 
	}
	//std::cout << " results with temporary weight " << results << "  " << temporary << std::endl; 
	if (temporary > results)
		return 0; 
	else return 1; 
	
}
std::vector<unsigned> threshold_synthesis_manager::weigths_assignment (std::vector<std::pair<boost::dynamic_bitset<>,boost::dynamic_bitset<>>> selected,std::vector<std::pair<unsigned,signed>> chow )
{
	const auto num_vars = tt_num_vars( func );
	std::vector<unsigned> weigths(num_vars,0u); 
	auto min_chow_prec = 0; 
	auto min_chow = evaluate_min_chow (chow); 
	auto temporary_w = 1u; 
	
	// bisogna creare le keys ;) 
	std::vector<bool> fatte(num_vars); 
	auto keys = create_keys (selected); 
	std::vector<std::pair<unsigned,unsigned>> save; 
	
	auto new_chow = chow;
	
	auto flag = 0; 
	while (new_chow.size() > 0)
	{
		
		//std::cout << " new chow size " << new_chow.size()  << std::endl; 
		//std::cout << " minimo chow " << min_chow << std::endl; 
		if ((min_chow == min_chow_prec) && (min_chow != 0) && (flag == 0))
		{
			//std::cout << " qui non devi entrare" << std::endl; 
			temporary_w = find_tempo(min_chow, save); 
		}
			
	for (auto & chows : chow)
	{
		//std::cout << " siamo alla variabile " << chows.first << std::endl; 
		if ((chows.second == min_chow) && (fatte[chows.first] == 0))
		{
			save.push_back({min_chow,temporary_w}); 
			std::vector<unsigned> good(num_vars,0u); 
			std::vector<unsigned> trace(num_vars,0u); 
				
			if (there_is (chows.first, keys) == 0)
			{
				
			for (auto & equas : keys)
			{
				
				if (equas.first == chows.first)
				{
					//std::cout << " keys" << std::endl; 
					++trace[equas.first]; 
					//std::cout << " blabla" << std::endl; 
					// bisogna vedere se vale equazione cerchi nelle keys chi ha la prima variabile uguale alla prima dei chow, se hai gia assegnato loro dei pesi, deve valere con i pesi, se vael va bene, altrimenti devi autmentrae 
					if (equation(equas.second.second, weigths, temporary_w) == 0) 
						{
							//std::cout << " blabla2" << std::endl;
							++good[equas.first]; 
							//std::cout << " quindi il probblemae quiii" << std::endl; 
							//continue; 
						}
					else 
						{
							//std::cout << " blabla33" << std::endl;
							flag = 1; 
							goto end; 
						}		
				}
			}
		
				
			if (good == trace)
			{
				//std::cout << " goodd == trace" << std::endl; 
				weigths[chows.first] = temporary_w; 
				fatte[chows.first] = 1; 
				const auto it = std::find(new_chow.begin(), new_chow.end(), chows); 
				new_chow.erase(it);
				flag = 0; 
				goto end; 
			}
			
		}
		else 
		{
			weigths[chows.first] = temporary_w; 
			fatte[chows.first] = 1; 
			//std::cout << " diamo alla variabile " << chows.first << " il valore di weight" << temporary_w << std::endl; 
			const auto it = std::find(new_chow.begin(), new_chow.end(), chows); 
			new_chow.erase(it);
		}
			
			
		}
	}
	
	end: 
	//std::cout << " size" << new_chow.size() << std::endl; 
	++temporary_w; 
	//min_chow_prec = min_chow; 
	min_chow_prec = min_chow;
	if (new_chow.size() > 0)
	min_chow = evaluate_min_chow(new_chow);
   }

   //std::cout << " fine di w" << std::endl; 
	return weigths; 
}
std::vector<unsigned> threshold_synthesis_manager::unate(tt& func)
{
  if ( verbose )
  {
    std::cout << "[i] check unateness" << std::endl;
  }
  properties_timer t( statistics, "unateness_runtime" );

  const auto num_vars = tt_num_vars( func );
  std::vector<unsigned> unate; // says unateness 

  for ( auto i = 0u; i < num_vars; ++i )
  {
	  auto y = i; //num_vars - i - 1; 
	  auto zero = tt_cof0(func,y); 
	  auto one = tt_cof1(func,y); 
	  
	  if  (zero.is_subset_of(one)) 
		 unate.push_back(0); 
	 	
   	  else if (one.is_subset_of(zero))
   		 {
			 unate.push_back(1);
			 //std::cout << " variable " << i << " is negative" << std::endl; 
		 }
	  
   	  else break;   
  }

  return unate;
}

std::vector<bool> threshold_synthesis_manager::nega_unate(tt& func)
{
  
  const auto num_vars = tt_num_vars( func );
  std::vector<bool> neg_unate; // says unateness 

  for ( auto i = 0u; i < num_vars; ++i )
  {
	  auto y = i; //num_vars - i - 1; 
	  auto zero = tt_cof0(func,y); 
	  auto one = tt_cof1(func,y); 
	  
	  if  (zero.is_subset_of(one)) 
		 neg_unate.push_back(0); 
	 	
   	  else if (one.is_subset_of(zero))
   		 {
			 neg_unate.push_back(1);
			 //std::cout << " variable " << i << " is negative" << std::endl; 
		 }
	  
   	  else break;   
  }

  return neg_unate;
}

std::vector<std::pair<unsigned,signed>> threshold_synthesis_manager::compute_chow (tt& func, std::vector<unsigned> un)
{
	std::vector<std::pair<unsigned,signed>> chow; 
	const auto num_vars = tt_num_vars( func );
	const auto t_size = func.size(); 
	boost::dynamic_bitset<> var(num_vars);
	
	for (auto x = 0; x<num_vars; ++x)
	{
		signed m = 0; 
		signed n = 0; 
		if (negative_unate(un,x) == 0)
		var = tt_nth_var(x).flip(); 
		else 
		var = tt_nth_var(x); 

		//std::cout << var << std::endl; 
		//std::cout << var[0] << std::endl; 
		for (auto i = 0; i <t_size; ++i)
		{
			if ((var[i] == 1) && (func[i] == 1))
				++m; 
			if ((var[i] == 0) && (func[i] == 1))
				++n;
		}
		signed p = 2 * (m - n); 
		chow.push_back({x,p}); 
		//std::cout << "variable " << x << " parameter " << p <<std::endl; 
	}
	
	std::vector<std::pair<unsigned,signed>> order_chow; 
	signed min = 0; 
	signed max = 0;
	
	for (auto &c : chow)
	{
		if (c.second < min)
			min = c.second; 
		if (c.second > max)
			max = c.second; 
	}
	for (auto x = min; x <= max; ++x)
	{
		for (auto &c : chow)
		{
			if (c.second == x)
				order_chow.push_back({c.first,c.second}); 
		}
	}
	return order_chow; 
}

bool threshold_synthesis_manager::check_all (std::vector<unsigned> weigths, tt& func, std::vector<unsigned> un)
{
	const auto num_vars = tt_num_vars( func );
	auto ineq = create_ALL_ineq (func, un); 
	for (auto & equ : ineq)
	{
		auto bigger = 0; 
		auto less = 0; 
		for ( auto x = 0; x <num_vars; ++x)
		{
			bigger = bigger + weigths[x]*equ.first[x]; 
			less = less + weigths[x]*equ.second[x]; 
		}
		if (bigger > less)
			continue; 
		else 
			return 1; 
	}
	return 0; 
}

unsigned threshold_synthesis_manager::compute_threshold (tt& func, std::vector<unsigned> weights, std::vector<unsigned> un)
{
	unsigned T = 0u; 
	std::vector< boost::dynamic_bitset<>> positive;
    const auto num_vars = tt_num_vars( func );

	for (auto x =0; x<func.size(); ++x)
	{
		
			boost::dynamic_bitset<> provv(num_vars), var(num_vars);
			for (auto j= 0; j <num_vars;++j)
			{
				auto y = j; //num_vars - j - 1; 
				if (negative_unate(un,y) == 0)
				var = tt_nth_var(y).flip(); 
				else
				var = tt_nth_var(y); 
				provv[j] = var[x]; 
			}
						 
			if (func[x] == 1)// greater part
				positive.push_back(provv); // bisogna vedere se funziona cosi 
	}
	
	auto result = 0u; 
	for (auto x = 0; x < num_vars; ++x)
	{
		result = result + positive[0][x]*weights[x]; 
	}
	T = result; 
	for (auto y = 1; y<positive.size(); ++y)
	{
		auto partial = 0u; 
		for (auto x = 0; x < num_vars; ++x)
		{
			partial = partial + positive[y][x]*weights[x]; 
		}
		if (partial < T)
			T = partial; 
	}
	
		//std::cout << "Threshold value = " << T << std::endl;
	
	
	return T; 
}


mign_graph threshold_synthesis_manager::create_threshold( tt& func, std::vector<unsigned> pos, std::vector<std::pair<unsigned,signed>> chow)
{
	auto c = create_ineq(func,pos); /// qui in input serve anche un, perche qundo ho quelli negati devo tenere conto di questo ;) 
	auto selected = select_ineq(func,c); 
	
		for (auto & d :c)
	{
		d.first.clear(); 
		d.second.clear(); 
	}
	
	if (verbose)
		std::cout << " [i] evaluating weights" << std::endl; 
	
	auto weigths = weigths_assignment(selected, chow); 
	
	for (auto & s :selected)
{
	s.first.clear(); 
	s.second.clear(); 
}

	mign_graph mign; 
	
	if (check_all(weigths,func,pos) == 0)
		
	{
		auto T = compute_threshold(func,weigths,pos); 
		std::cout << "Threshold value = " << T << std::endl;
		const auto num_vars = tt_num_vars( func );
		print_w (weigths); 
		//std::cout << " assertion failed qui" << std::endl; 
		mign = threshold (num_vars, T, 1,weigths);
	}
	else 
		{ 
			auto T = 0; 
			const auto num_vars = tt_num_vars( func );
			mign = threshold (num_vars, T, 1,weigths);
			std::cout << " Threshold function cannot be realized. " << std::endl; 
		}
	
	return mign; 
}

resulting threshold_synthesis_manager::create_threshold_T_w( tt& func, std::vector<unsigned> pos, std::vector<std::pair<unsigned,signed>> chow)
{
	//std::cout << " create ineq" << std::endl; 
	auto c = create_ineq(func,pos); /// qui in input serve anche un, perche qundo ho quelli negati devo tenere conto di questo ;) 
	
	//std::cout << " create selectd ineq" << std::endl; 
	auto selected = select_ineq(func,c); 
	
	//if (verbose)
	//std::cout << " [i] evaluating weights" << std::endl;
	 
	 //std::cout << " weights" << std::endl; 
	auto weigths = weigths_assignment(selected, chow); 

	resulting result; 
	auto T = 0; 
	
	//std::cout << "check" << std::endl; 
	if (check_all(weigths,func,pos) == 0)
		
	{
		T = compute_threshold(func,weigths,pos); 
	}
	else 
	{ 
	    T = 0; 
	    //std::cout << " Threshold function cannot be realized. " << std::endl; 
	}
	
	//neg_unate = nega_unate(func); 
	result.t_and_w = std::make_pair(T, weigths); 
	result.nega_un = nega_unate(func); 
	return result; 
}

/******************************************************************************
 * Public functions                                                           *
 ******************************************************************************/

mign_graph threshold_synthesis( const tt& func,
                           const properties::ptr& settings,
                           const properties::ptr& statistics )
{
	/* check unateness*/ 
    threshold_synthesis_manager mgr( func, settings, statistics );
    return mgr.run();
}


resulting compute_T_w ( tt& func)
{
	auto settings = std::make_shared<properties>();
	auto statistics = std::make_shared<properties>();
    threshold_synthesis_manager mgr( func, settings, statistics );
    return mgr.run_T_w();
	
}

}

// Local Variables:
// c-basic-offset: 2
// eval: (c-set-offset 'substatement-open 0)
// eval: (c-set-offset 'innamespace 0)
// End:
