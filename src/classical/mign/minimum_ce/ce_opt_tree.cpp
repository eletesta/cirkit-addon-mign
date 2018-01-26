#include "ce_opt_tree.hpp"

#include <core/utils/timer.hpp>

#include <boost/assign/std/vector.hpp>
#include <boost/format.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/graph/topological_sort.hpp>
#include <boost/range/algorithm.hpp>
#include <boost/range/iterator_range.hpp>
#include <boost/regex.hpp>

using namespace boost::assign;

using boost::format;

namespace cirkit
{

	bool e_output (mig_graph old_mig, mig_node node)
	{
		bool output_ce = false; 
		for ( const auto& output : mig_info( old_mig ).outputs )
		{ 
			if (output.first.node == node)
			{
				output_ce = true; 
				break; 
			}
		}
		return output_ce; 
	}
		
	unsigned number_indegree (mig_node node,mig_graph mig)
		{
			auto count = 0u; 
			for (const auto& node_bis : boost::make_iterator_range(boost::vertices(mig)))
			{
				if ( !boost::out_degree( node_bis, mig ) ) { continue; }
				auto children = get_children(mig,node_bis); 
				for ( auto x = 0; x<3; ++x)
				{
					if (children[x].node == node)
						count++; 
				}
			}
			
			for (const auto& output : mig_info(mig).outputs)
			{
				if (output.first.node == node)
					count++; 
			}
			return count; 
		}

	
std::map<mig_node, std::vector<std::pair<mig_function,unsigned>>> init_visited_table_ce( const mig_graph& mig, mig_graph& mig_new )
{
 std::map<mig_node, std::vector<std::pair<mig_function,unsigned>>> old_to_new;

 old_to_new[0].push_back(std::make_pair(mig_get_constant( mig_new, false ),0));
for ( const auto& pi : mig_info( mig ).inputs )
	 {
	   old_to_new[pi].push_back(std::make_pair(mig_create_pi( mig_new, mig_info( mig ).node_names.at( pi ) ),0));
	 }

 return old_to_new;
	}
	
	unsigned sum_ce(mig_function& a, mig_function& b, mig_function& c, bool outcompl)
	{
		unsigned sum = 0u; 
		if (a.node != 0)
		{
			if (b.node != 0)
			{
				if (c.node != 0)
					sum = a.complemented + b.complemented + c.complemented ; //+ outcompl;
				else 
					sum = a.complemented + b.complemented ;// + outcompl;
			}
			else if (b.node == 0)
			{
					sum = a.complemented + c.complemented;// + outcompl;
			}
		}
		else 
		sum = b.complemented + c.complemented ; //+ outcompl;
			
		return sum; 
	}

	unsigned sum_ce_changed (mig_function& a, mig_function& b, mig_function& c, bool outcompl)
	{
		unsigned sum = 0u; 
		if (a.node != 0)
		{
			if (b.node != 0)
			{
				if (c.node != 0)
					sum = !a.complemented + !b.complemented + !c.complemented ; //+ !outcompl;
				else 
					sum = !a.complemented + !b.complemented ;// + !outcompl;
			}
			else if (b.node == 0)
			{
					sum = !a.complemented + !c.complemented; // + !outcompl;
			}
		}
		else 
		sum = !b.complemented + !c.complemented; // + !outcompl;
			
		return sum; 
	}
	
	
std::vector<std::pair<mig_function,unsigned>> comple (std::vector<std::pair<mig_function,unsigned>> ciao, bool complemented)
		{
				auto size = ciao.size(); 
				std::vector<std::pair<mig_function,unsigned>> schifo(size);
				for (auto i = 0; i< size; ++i)
				{
					schifo[i].first.node = ciao[i].first.node; 
					schifo[i].second = ciao[i].second; 
					schifo[i].first.complemented = ciao[i].first.complemented != complemented; 
				}
			
				return schifo; 
			}
			
			
			std::vector<std::pair<mig_function,unsigned>> a_double_last (mig_graph& mig_new, std::vector<std::pair<mig_function,unsigned>> a, std::vector<std::pair<mig_function,unsigned>> b, std::vector<std::pair<mig_function,unsigned>> c, bool outcompl)
			{
				std::vector<std::pair<mig_function,unsigned>> node; 
	
				auto sum_start = sum_ce(a[0].first,b[0].first,c[0].first,outcompl) + a[0].second + b[0].second + c[0].second + outcompl; 
				auto sum_final = sum_ce_changed (a[0].first,b[0].first,c[0].first,outcompl) + a[0].second + b[0].second + c[0].second + !outcompl; 
				auto sum_double = sum_ce(a[1].first,b[0].first,c[0].first,outcompl) + a[1].second + b[0].second + c[0].second + outcompl; 
				auto sum_change_double = sum_ce_changed (a[1].first,b[0].first,c[0].first,outcompl) + a[1].second + b[0].second + c[0].second + !outcompl; 

	
				auto best_pos = sum_start;
				auto best_neg = sum_final; 
				std::vector<unsigned> flag(2);    
				if (sum_start < sum_double)
				{
					best_pos = sum_start;
					flag[0] = 1; 
				}
				else 
				{ 
					best_pos = sum_double;
					flag[0] = 2;  
				}
	
				if (sum_final < sum_change_double)
				{
					best_neg = sum_final; 
					flag[1] = 1; 
				}
				else 
				{ 
					best_neg = sum_change_double;
					flag[1] = 2; 
				}
	
				if ( best_pos < best_neg)
				{
						if ( flag[0] == 1)
							{
										auto s = best_pos ; //- outcompl; 
												node.push_back(std::make_pair(mig_create_maj(mig_new,a[0].first,b[0].first,c[0].first),s)); 
													}
														if (flag[0] == 2)
															{
																		auto s = best_pos- outcompl; 
																				node.push_back(std::make_pair(mig_create_maj(mig_new,a[1].first,b[0].first,c[0].first),s)); 
																			} 
																		}
																		else if (best_pos > best_neg)
																		{
																				if ( flag[1] == 1)
																					{
																								auto s = best_neg - !outcompl; 
																										node.push_back(std::make_pair(mig_create_maj(mig_new,!a[0].first,!b[0].first,!c[0].first),s)); 
																												node[0].first.complemented = 1; 
																													}
																														if (flag[1] == 2)
																															{
																																		auto s = best_neg - !outcompl; 
																																				node.push_back(std::make_pair(mig_create_maj(mig_new,!a[1].first,!b[0].first,!c[0].first),s)); 
																																						node[0].first.complemented = 1; 
																																							}
																																						} 
																																						else if ( best_pos == best_neg)
																																						{
					if ( flag[0] == 1)
					{
						auto s = best_pos - outcompl; 
						node.push_back(std::make_pair(mig_create_maj(mig_new,a[0].first,b[0].first,c[0].first),s)); 
					} 
					else if (flag[0] == 2)
					{
						auto s = best_pos - outcompl; 
						node.push_back(std::make_pair(mig_create_maj(mig_new,a[1].first,b[0].first,c[0].first),s)); 
					} 
					if ( flag[1] == 1)
					{
						auto s = best_neg - !outcompl; 
						node.push_back(std::make_pair(mig_create_maj(mig_new,!a[0].first,!b[0].first,!c[0].first),s)); 
						node[1].first.complemented = 1; 
					}
					else if (flag[1] == 2)
					{
						auto s = best_neg - !outcompl; 
						node.push_back(std::make_pair(mig_create_maj(mig_new,!a[1].first,!b[0].first,!c[0].first),s)); 
						node[1].first.complemented = 1; 
					}
		
				}
				return node; 
	
			}

			std::vector<std::pair<mig_function,unsigned>> b_double_last (mig_graph& mig_new, std::vector<std::pair<mig_function,unsigned>> a, std::vector<std::pair<mig_function,unsigned>> b, std::vector<std::pair<mig_function,unsigned>> c, bool outcompl)
			{
				std::vector<std::pair<mig_function,unsigned>> node; 
	
				auto sum_start = sum_ce(a[0].first,b[0].first,c[0].first,outcompl) + a[0].second + b[0].second + c[0].second + outcompl; 
				auto sum_final = sum_ce_changed (a[0].first,b[0].first,c[0].first,outcompl) + a[0].second + b[0].second + c[0].second + !outcompl; 
				auto sum_double = sum_ce(a[0].first,b[1].first,c[0].first,outcompl) + a[0].second + b[1].second + c[0].second + outcompl; 
				auto sum_change_double = sum_ce_changed (a[0].first,b[1].first,c[0].first,outcompl) + a[0].second + b[1].second + c[0].second + !outcompl; 

	
				auto best_pos = sum_start;
				auto best_neg = sum_final; 
				std::vector<unsigned> flag(2);    
				if (sum_start < sum_double)
				{
					best_pos = sum_start;
					flag[0] = 1; 
				}
				else 
				{ 
					best_pos = sum_double;
					flag[0] = 2;  
				}
	
				if (sum_final < sum_change_double)
				{
					best_neg = sum_final; 
					flag[1] = 1; 
				}
				else 
				{ 
					best_neg = sum_change_double;
					flag[1] = 2; 
				}
	
				if ( best_pos < best_neg)
				{
					if ( flag[0] == 1)
					{
						auto s = best_pos - outcompl; 
						node.push_back(std::make_pair(mig_create_maj(mig_new,a[0].first,b[0].first,c[0].first),s)); 
					}
					if (flag[0] == 2)
					{
						auto s = best_pos - outcompl; 
						node.push_back(std::make_pair(mig_create_maj(mig_new,a[0].first,b[1].first,c[0].first),s)); 
					} 
				}
				else if (best_pos > best_neg)
				{
					if ( flag[1] == 1)
					{
						auto s = best_neg - !outcompl; 
						node.push_back(std::make_pair(mig_create_maj(mig_new,!a[0].first,!b[0].first,!c[0].first),s)); 
						node[0].first.complemented = 1; 
					}
					if (flag[1] == 2)
					{
						auto s = best_neg - !outcompl; 
						node.push_back(std::make_pair(mig_create_maj(mig_new,!a[0].first,!b[1].first,!c[0].first),s)); 
						node[0].first.complemented = 1; 
					}
				} 
				else if ( best_pos == best_neg)
				{
					if ( flag[0] == 1)
					{
						auto s = best_pos - outcompl; 
						node.push_back(std::make_pair(mig_create_maj(mig_new,a[0].first,b[0].first,c[0].first),s)); 
					} 
					else if (flag[0] == 2)
					{
						auto s = best_pos - outcompl; 
						node.push_back(std::make_pair(mig_create_maj(mig_new,a[0].first,b[1].first,c[0].first),s)); 
					} 
					if ( flag[1] == 1)
					{
						auto s = best_neg - !outcompl; 
						node.push_back(std::make_pair(mig_create_maj(mig_new,!a[0].first,!b[0].first,!c[0].first),s)); 
						node[1].first.complemented = 1; 
					}
					else if (flag[1] == 2)
					{
						auto s = best_neg - !outcompl; 
						node.push_back(std::make_pair(mig_create_maj(mig_new,!a[0].first,!b[1].first,!c[0].first),s)); 
						node[1].first.complemented = 1; 
					}
		
				}
				return node; 
	
			}

			std::vector<std::pair<mig_function,unsigned>> c_double_last (mig_graph& mig_new, std::vector<std::pair<mig_function,unsigned>> a, std::vector<std::pair<mig_function,unsigned>> b, std::vector<std::pair<mig_function,unsigned>> c, bool outcompl)
			{
				std::vector<std::pair<mig_function,unsigned>> node; 
	
				auto sum_start = sum_ce(a[0].first,b[0].first,c[0].first,outcompl) + a[0].second + b[0].second + c[0].second + outcompl; 
				auto sum_final = sum_ce_changed (a[0].first,b[0].first,c[0].first,outcompl) + a[0].second + b[0].second + c[0].second + !outcompl; 
				auto sum_double = sum_ce(a[0].first,b[0].first,c[1].first,outcompl) + a[0].second + b[0].second + c[1].second + outcompl; 
				auto sum_change_double = sum_ce_changed (a[0].first,b[0].first,c[1].first,outcompl) + a[0].second + b[0].second + c[1].second + !outcompl; 

	
				auto best_pos = sum_start;
				auto best_neg = sum_final; 
				std::vector<unsigned> flag(2);    
				if (sum_start < sum_double)
				{
					best_pos = sum_start;
					flag[0] = 1; 
				}
				else 
				{ 
					best_pos = sum_double;
					flag[0] = 2;  
				}
	
				if (sum_final < sum_change_double)
				{
					best_neg = sum_final; 
					flag[1] = 1; 
				}
				else 
				{ 
					best_neg = sum_change_double;
					flag[1] = 2; 
				}
	
				if ( best_pos < best_neg)
				{
					if ( flag[0] == 1)
					{
						auto s = best_pos - outcompl; 
						node.push_back(std::make_pair(mig_create_maj(mig_new,a[0].first,b[0].first,c[0].first),s)); 
					}
					if (flag[0] == 2)
					{
						auto s = best_pos - outcompl; 
						node.push_back(std::make_pair(mig_create_maj(mig_new,a[0].first,b[0].first,c[1].first),s)); 
					} 
				}
				else if (best_pos > best_neg)
				{
					if ( flag[1] == 1)
					{
						auto s = best_neg - !outcompl; 
						node.push_back(std::make_pair(mig_create_maj(mig_new,!a[0].first,!b[0].first,!c[0].first),s)); 
						node[0].first.complemented = 1; 
					}
					if (flag[1] == 2)
					{
						auto s = best_neg - !outcompl; 
						node.push_back(std::make_pair(mig_create_maj(mig_new,!a[0].first,!b[0].first,!c[1].first),s)); 
						node[0].first.complemented = 1; 
					}
				} 
				else if ( best_pos == best_neg)
				{
					if ( flag[0] == 1)
					{
						auto s = best_pos - outcompl; 
						node.push_back(std::make_pair(mig_create_maj(mig_new,a[0].first,b[0].first,c[0].first),s)); 
					} 
					else if (flag[0] == 2)
					{
						auto s = best_pos - outcompl; 
						node.push_back(std::make_pair(mig_create_maj(mig_new,a[0].first,b[0].first,c[1].first),s)); 
					} 
					if ( flag[1] == 1)
					{
						auto s = best_neg - !outcompl; 
						node.push_back(std::make_pair(mig_create_maj(mig_new,!a[0].first,!b[0].first,!c[0].first),s)); 
						node[1].first.complemented = 1; 
					}
					else if (flag[1] == 2)
					{
						auto s = best_neg - !outcompl; 
						node.push_back(std::make_pair(mig_create_maj(mig_new,!a[0].first,!b[0].first,!c[1].first),s)); 
						node[1].first.complemented = 1; 
					}
		
				}
				return node; 
			}

			std::vector<std::pair<mig_function,unsigned>> ab_double_last (mig_graph& mig_new, std::vector<std::pair<mig_function,unsigned>> a, std::vector<std::pair<mig_function,unsigned>> b, std::vector<std::pair<mig_function,unsigned>> c, bool outcompl)
			{
				std::vector<std::pair<mig_function,unsigned>> node; 
	
				auto sum_start = sum_ce(a[0].first,b[0].first,c[0].first,outcompl) + a[0].second + b[0].second + c[0].second + outcompl; // caso 1
				auto sum_final = sum_ce_changed(a[0].first,b[0].first,c[0].first,outcompl) + a[0].second + b[0].second + c[0].second + !outcompl; 
				auto sum_double_a = sum_ce(a[1].first,b[0].first,c[0].first, outcompl) + a[1].second + b[0].second + c[0].second + outcompl; // caso 2
				auto sum_change_a = sum_ce_changed(a[1].first,b[0].first,c[0].first,outcompl) + a[1].second + b[0].second + c[0].second + !outcompl; 
				auto sum_double_b = sum_ce(a[0].first,b[1].first,c[0].first, outcompl)+ a[0].second + b[1].second + c[0].second + outcompl; //caso 3
				auto sum_change_b = sum_ce_changed(a[0].first,b[1].first,c[0].first,outcompl) + a[0].second + b[1].second + c[0].second + !outcompl; 
				auto sum_double_ab = sum_ce(a[1].first,b[1].first,c[0].first, outcompl)+ a[1].second + b[1].second + c[0].second + outcompl; // caso 4
				auto sum_change_ab = sum_ce_changed(a[1].first,b[1].first,c[0].first,outcompl)+ a[1].second + b[1].second + c[0].second + !outcompl; 
	
				auto best_pos = sum_start;
				auto best_neg = sum_final; 
				std::vector<unsigned> flag(2);  
				flag[0] = 1; 
				flag[1] = 1;   
				if (sum_double_a < best_pos)
				{
					best_pos = sum_double_a;
					flag[0] = 2; 
				}
				if (sum_double_b < best_pos)
				{
					best_pos = sum_double_b;
					flag[0] = 3; 
				}
				if (sum_double_ab < best_pos)
				{
					best_pos = sum_double_ab;
					flag[0] = 4; 
				}
	
				if (sum_change_a < best_neg)
				{
					best_neg = sum_change_a;
					flag[1] = 2; 
				}
				if (sum_change_b < best_neg)
				{
					best_neg = sum_change_b;
					flag[1] = 3; 
				}
				if (sum_change_ab < best_neg)
				{
					best_neg = sum_change_ab;
					flag[1] = 4; 
				}
	
				if ( best_pos < best_neg)
				{
					if ( flag[0] == 1)
					{
						auto s = best_pos - outcompl; 
						node.push_back(std::make_pair(mig_create_maj(mig_new,a[0].first,b[0].first,c[0].first),s)); 
					}
		
					if (flag[0] == 2)
					{
						auto s = best_pos - outcompl; 
						node.push_back(std::make_pair(mig_create_maj(mig_new,a[1].first,b[0].first,c[0].first),s)); 
					}
					if (flag[0] == 3)
					{
						auto s = best_pos - outcompl; 
						node.push_back(std::make_pair(mig_create_maj(mig_new,a[0].first,b[1].first,c[0].first),s)); 
					} 
					if (flag[0] == 4)
					{
						auto s = best_pos - outcompl; 
						node.push_back(std::make_pair(mig_create_maj(mig_new,a[1].first,b[1].first,c[0].first),s)); 
					} 
				}
				else if (best_pos > best_neg)
				{
					if ( flag[1] == 1)
					{
						auto s = best_neg - !outcompl; 
						node.push_back(std::make_pair(mig_create_maj(mig_new,!a[0].first,!b[0].first,!c[0].first),s)); 
						node[0].first.complemented = 1; 
					}
					if (flag[1] == 2)
					{
						auto s = best_neg - !outcompl; 
						node.push_back(std::make_pair(mig_create_maj(mig_new,!a[1].first,!b[0].first,!c[0].first),s)); 
						node[0].first.complemented = 1; 
					}
					if ( flag[1] == 3)
					{
						auto s = best_neg - !outcompl; 
						node.push_back(std::make_pair(mig_create_maj(mig_new,!a[0].first,!b[1].first,!c[0].first),s)); 
						node[0].first.complemented = 1; 
					}
					if (flag[1] == 4)
					{
						auto s = best_neg - !outcompl; 
						node.push_back(std::make_pair(mig_create_maj(mig_new,!a[1].first,!b[1].first,!c[0].first),s)); 
						node[0].first.complemented = 1; 
					}
				} 
				else if ( best_pos == best_neg)
				{
					if ( flag[0] == 1)
					{
						auto s = best_pos - outcompl; 
						node.push_back(std::make_pair(mig_create_maj(mig_new,a[0].first,b[0].first,c[0].first),s)); 
					}
		
					else if (flag[0] == 2)
					{
						auto s = best_pos - outcompl; 
						node.push_back(std::make_pair(mig_create_maj(mig_new,a[1].first,b[0].first,c[0].first),s)); 
					}
					else if (flag[0] == 3)
					{
						auto s = best_pos - outcompl; 
						node.push_back(std::make_pair(mig_create_maj(mig_new,a[0].first,b[1].first,c[0].first),s)); 
					} 
					else if (flag[0] == 4)
					{
						auto s = best_pos - outcompl; 
						node.push_back(std::make_pair(mig_create_maj(mig_new,a[1].first,b[1].first,c[0].first),s)); 
					} 
					if ( flag[1] == 1)
					{
						auto s = best_neg - !outcompl; 
						node.push_back(std::make_pair(mig_create_maj(mig_new,!a[0].first,!b[0].first,!c[0].first),s)); 
						node[1].first.complemented = 1; 
					}
					else if (flag[1] == 2)
					{
						auto s = best_neg - !outcompl; 
						node.push_back(std::make_pair(mig_create_maj(mig_new,!a[1].first,!b[0].first,!c[0].first),s)); 
						node[1].first.complemented = 1; 
					}
					else if ( flag[1] == 3)
					{
						auto s = best_neg - !outcompl; 
						node.push_back(std::make_pair(mig_create_maj(mig_new,!a[0].first,!b[1].first,!c[0].first),s)); 
						node[1].first.complemented = 1; 
					}
					else if (flag[1] == 4)
					{
						auto s = best_neg - !outcompl; 
						node.push_back(std::make_pair(mig_create_maj(mig_new,!a[1].first,!b[1].first,!c[0].first),s)); 
						node[1].first.complemented = 1; 
					}
		
				}
				return node; 
			}

			std::vector<std::pair<mig_function,unsigned>> ac_double_last (mig_graph& mig_new, std::vector<std::pair<mig_function,unsigned>> a, std::vector<std::pair<mig_function,unsigned>> b, std::vector<std::pair<mig_function,unsigned>> c, bool outcompl)
			{
				std::vector<std::pair<mig_function,unsigned>> node; 
	
				auto sum_start = sum_ce(a[0].first,b[0].first,c[0].first,outcompl) + a[0].second + b[0].second + c[0].second + outcompl; // caso 1
				auto sum_final = sum_ce_changed(a[0].first,b[0].first,c[0].first,outcompl) + a[0].second + b[0].second + c[0].second + !outcompl; 
				auto sum_double_a = sum_ce(a[1].first,b[0].first,c[0].first, outcompl) + a[1].second + b[0].second + c[0].second + outcompl; // caso 2
				auto sum_change_a = sum_ce_changed(a[1].first,b[0].first,c[0].first,outcompl) + a[1].second + b[0].second + c[0].second + !outcompl; 
				auto sum_double_b = sum_ce(a[0].first,b[0].first,c[1].first, outcompl)+ a[0].second + b[0].second + c[1].second + outcompl; //caso 3
				auto sum_change_b = sum_ce_changed(a[0].first,b[0].first,c[1].first,outcompl) + a[0].second + b[0].second + c[1].second + !outcompl; 
				auto sum_double_ab = sum_ce(a[1].first,b[0].first,c[1].first, outcompl)+ a[1].second + b[0].second + c[1].second + outcompl; // caso 4
				auto sum_change_ab = sum_ce_changed(a[1].first,b[0].first,c[1].first,outcompl)+ a[1].second + b[0].second + c[1].second + !outcompl; 
	
				auto best_pos = sum_start;
				auto best_neg = sum_final; 
				std::vector<unsigned> flag(2);  
				flag[0] = 1; 
				flag[1] = 1;   
				if (sum_double_a < best_pos)
				{
					best_pos = sum_double_a;
					flag[0] = 2; 
				}
				if (sum_double_b < best_pos)
				{
					best_pos = sum_double_b;
					flag[0] = 3; 
				}
				if (sum_double_ab < best_pos)
				{
					best_pos = sum_double_ab;
					flag[0] = 4; 
				}
	
				if (sum_change_a < best_neg)
				{
					best_neg = sum_change_a;
					flag[1] = 2; 
				}
				if (sum_change_b < best_neg)
				{
					best_neg = sum_change_b;
					flag[1] = 3; 
				}
				if (sum_change_ab < best_neg)
				{
					best_neg = sum_change_ab;
					flag[1] = 4; 
				}
	
				if ( best_pos < best_neg)
				{
					if ( flag[0] == 1)
					{
						auto s = best_pos - outcompl; 
						node.push_back(std::make_pair(mig_create_maj(mig_new,a[0].first,b[0].first,c[0].first),s)); 
					}
		
					if (flag[0] == 2)
					{
						auto s = best_pos - outcompl; 
						node.push_back(std::make_pair(mig_create_maj(mig_new,a[1].first,b[0].first,c[0].first),s)); 
					}
					if (flag[0] == 3)
					{
						auto s = best_pos - outcompl; 
						node.push_back(std::make_pair(mig_create_maj(mig_new,a[0].first,b[0].first,c[1].first),s)); 
					} 
					if (flag[0] == 4)
					{
						auto s = best_pos - outcompl; 
						node.push_back(std::make_pair(mig_create_maj(mig_new,a[1].first,b[0].first,c[1].first),s)); 
					} 
				}
				else if (best_pos > best_neg)
				{
					if ( flag[1] == 1)
					{
						auto s = best_neg - !outcompl; 
						node.push_back(std::make_pair(mig_create_maj(mig_new,!a[0].first,!b[0].first,!c[0].first),s)); 
						node[0].first.complemented = 1; 
					}
					if (flag[1] == 2)
					{
						auto s = best_neg - !outcompl; 
						node.push_back(std::make_pair(mig_create_maj(mig_new,!a[1].first,!b[0].first,!c[0].first),s)); 
						node[0].first.complemented = 1; 
					}
					if ( flag[1] == 3)
					{
						auto s = best_neg - !outcompl; 
						node.push_back(std::make_pair(mig_create_maj(mig_new,!a[0].first,!b[0].first,!c[0].first),s)); 
						node[0].first.complemented = 1; 
					}
					if (flag[1] == 4)
					{
						auto s = best_neg - !outcompl; 
						node.push_back(std::make_pair(mig_create_maj(mig_new,!a[1].first,!b[0].first,!c[1].first),s)); 
						node[0].first.complemented = 1; 
					}
				} 
				else if ( best_pos == best_neg)
				{
					if ( flag[0] == 1)
					{
						auto s = best_pos - outcompl; 
						node.push_back(std::make_pair(mig_create_maj(mig_new,a[0].first,b[0].first,c[0].first),s)); 
					}
		
					else if (flag[0] == 2)
					{
						auto s = best_pos - outcompl; 
						node.push_back(std::make_pair(mig_create_maj(mig_new,a[1].first,b[0].first,c[0].first),s)); 
					}
					else if (flag[0] == 3)
					{
						auto s = best_pos - outcompl; 
						node.push_back(std::make_pair(mig_create_maj(mig_new,a[0].first,b[0].first,c[1].first),s)); 
					} 
					else if (flag[0] == 4)
					{
						auto s = best_pos - outcompl; 
						node.push_back(std::make_pair(mig_create_maj(mig_new,a[1].first,b[0].first,c[1].first),s)); 
					} 
					if ( flag[1] == 1)
					{
						auto s = best_neg - !outcompl; 
						node.push_back(std::make_pair(mig_create_maj(mig_new,!a[0].first,!b[0].first,!c[0].first),s)); 
						node[1].first.complemented = 1; 
					}
					else if (flag[1] == 2)
					{
						auto s = best_neg - !outcompl; 
						node.push_back(std::make_pair(mig_create_maj(mig_new,!a[1].first,!b[0].first,!c[0].first),s)); 
						node[1].first.complemented = 1; 
					}
					else if ( flag[1] == 3)
					{
						auto s = best_neg - !outcompl; 
						node.push_back(std::make_pair(mig_create_maj(mig_new,!a[0].first,!b[0].first,!c[1].first),s)); 
						node[1].first.complemented = 1; 
					}
					else if (flag[1] == 4)
					{
						auto s = best_neg - !outcompl; 
						node.push_back(std::make_pair(mig_create_maj(mig_new,!a[1].first,!b[0].first,!c[1].first),s)); 
						node[1].first.complemented = 1; 
					}
		
				}
				return node; 
			}

			std::vector<std::pair<mig_function,unsigned>> bc_double_last (mig_graph& mig_new, std::vector<std::pair<mig_function,unsigned>> a, std::vector<std::pair<mig_function,unsigned>> b, std::vector<std::pair<mig_function,unsigned>> c, bool outcompl)
			{
				std::vector<std::pair<mig_function,unsigned>> node; 
				
				//std::cout << " siamo dentro bc double" << std::endl;
	
				auto sum_start = sum_ce(a[0].first,b[0].first,c[0].first,outcompl) + a[0].second + b[0].second + c[0].second + outcompl; // caso 1
				auto sum_final = sum_ce_changed(a[0].first,b[0].first,c[0].first,outcompl) + a[0].second + b[0].second + c[0].second+ !outcompl; 
				auto sum_double_a = sum_ce(a[0].first,b[0].first,c[1].first, outcompl) + a[0].second + b[0].second + c[1].second+ outcompl; // caso 2
				auto sum_change_a = sum_ce_changed(a[0].first,b[0].first,c[1].first,outcompl) + a[0].second + b[0].second + c[1].second+ !outcompl; 
				auto sum_double_b = sum_ce(a[0].first,b[1].first,c[0].first, outcompl)+ a[0].second + b[1].second + c[0].second+ outcompl; //caso 3
				auto sum_change_b = sum_ce_changed(a[0].first,b[1].first,c[0].first,outcompl) + a[0].second + b[1].second + c[0].second+ !outcompl; 
				auto sum_double_ab = sum_ce(a[0].first,b[1].first,c[1].first, outcompl)+ a[0].second + b[1].second + c[1].second+ outcompl; // caso 4
				auto sum_change_ab = sum_ce_changed(a[0].first,b[1].first,c[1].first,outcompl)+ a[0].second + b[1].second + c[1].second+ !outcompl; 
				
				//std::cout << " sum caso bph " << sum_start << std::endl; 
				//std::cout << " sum caso bph " << sum_final << std::endl; 
				//std::cout << " sum caso bph " << sum_double_a << std::endl; 
				//std::cout << " sum caso bph " << sum_change_a << std::endl; 
				//std::cout << " sum caso bph " << sum_double_b << std::endl; 
				//std::cout << " sum caso bph " << sum_change_b << std::endl; 
				//std::cout << " sum caso bph " << sum_double_ab << std::endl; 
				//std::cout << " sum caso bph " << sum_change_ab << std::endl; 
				
				//std::cout << " a0 . second = " << a[0].second; 
				//std::cout << "b[0] second " << b[0].second; 
				//std::cout << "b[1] second " << b[1].second; 
				//std::cout << "c[0] second " << c[0].second; 
				//std::cout << "c[1] second " << c[1].second; 
				
				auto best_pos = sum_start;
				auto best_neg = sum_final; 
				std::vector<unsigned> flag(2);  
				flag[0] = 1; 
				flag[1] = 1;   
				if (sum_double_a < best_pos)
				{
					best_pos = sum_double_a;
					flag[0] = 2; 
				}
				if (sum_double_b < best_pos)
				{
					best_pos = sum_double_b;
					flag[0] = 3; 
				}
				if (sum_double_ab < best_pos)
				{
					best_pos = sum_double_ab;
					flag[0] = 4; 
				}
	
				if (sum_change_a < best_neg)
				{
					best_neg = sum_change_a;
					flag[1] = 2; 
				}
				if (sum_change_b < best_neg)
				{
					best_neg = sum_change_b;
					flag[1] = 3; 
				}
				if (sum_change_ab < best_neg)
				{
					best_neg = sum_change_ab;
					flag[1] = 4; 
				}
	
	
				//std::cout << " best pos =  " << best_pos << std::endl; 
				//std::cout << " best neg =  " << best_neg << std::endl; 
				
				if ( best_pos < best_neg)
				{
					if ( flag[0] == 1)
					{
						auto s = best_pos - outcompl; 
						node.push_back(std::make_pair(mig_create_maj(mig_new,a[0].first,b[0].first,c[0].first),s)); 
					}
		
					if (flag[0] == 2)
					{
						auto s = best_pos - outcompl; 
						node.push_back(std::make_pair(mig_create_maj(mig_new,a[0].first,b[0].first,c[1].first),s)); 
					}
					if (flag[0] == 3)
					{
						auto s = best_pos - outcompl; 
						node.push_back(std::make_pair(mig_create_maj(mig_new,a[0].first,b[1].first,c[0].first),s)); 
					} 
					if (flag[0] == 4)
					{
						auto s = best_pos - outcompl; 
						node.push_back(std::make_pair(mig_create_maj(mig_new,a[0].first,b[1].first,c[1].first),s)); 
					} 
				}
				else if (best_pos > best_neg)
				{
					if ( flag[1] == 1)
					{
						auto s = best_neg - !outcompl; 
						node.push_back(std::make_pair(mig_create_maj(mig_new,!a[0].first,!b[0].first,!c[0].first),s)); 
						node[0].first.complemented = 1; 
					}
					if (flag[1] == 2)
					{
						auto s = best_neg - !outcompl; 
						node.push_back(std::make_pair(mig_create_maj(mig_new,!a[0].first,!b[0].first,!c[1].first),s)); 
						node[0].first.complemented = 1; 
					}
					if ( flag[1] == 3)
					{
						auto s = best_neg - !outcompl; 
						node.push_back(std::make_pair(mig_create_maj(mig_new,!a[0].first,!b[1].first,!c[0].first),s)); 
						node[0].first.complemented = 1; 
					}
					if (flag[1] == 4)
					{
						auto s = best_neg - !outcompl; 
						node.push_back(std::make_pair(mig_create_maj(mig_new,!a[0].first,!b[1].first,!c[1].first),s)); 
						node[0].first.complemented = 1; 
					}
				} 
				else if ( best_pos == best_neg)
				{
					if ( flag[0] == 1)
					{
						auto s = best_pos - outcompl; 
						node.push_back(std::make_pair(mig_create_maj(mig_new,a[0].first,b[0].first,c[0].first),s)); 
					}
		
					else if (flag[0] == 2)
					{
						auto s = best_pos - outcompl; 
						node.push_back(std::make_pair(mig_create_maj(mig_new,a[0].first,b[0].first,c[1].first),s)); 
					}
					else if (flag[0] == 3)
					{
						auto s = best_pos - outcompl; 
						node.push_back(std::make_pair(mig_create_maj(mig_new,a[0].first,b[1].first,c[0].first),s)); 
					} 
					else if (flag[0] == 4)
					{
						auto s = best_pos - outcompl; 
						node.push_back(std::make_pair(mig_create_maj(mig_new,a[0].first,b[1].first,c[1].first),s)); 
					} 
					if ( flag[1] == 1)
					{
						auto s = best_neg - !outcompl; 
						node.push_back(std::make_pair(mig_create_maj(mig_new,!a[0].first,!b[0].first,!c[0].first),s)); 
						node[1].first.complemented = 1; 
					}
					else if (flag[1] == 2)
					{
						auto s = best_neg - !outcompl; 
						node.push_back(std::make_pair(mig_create_maj(mig_new,!a[0].first,!b[0].first,!c[1].first),s)); 
						node[1].first.complemented = 1; 
					}
					else if ( flag[1] == 3)
					{
						auto s = best_neg - !outcompl; 
						node.push_back(std::make_pair(mig_create_maj(mig_new,!a[0].first,!b[1].first,!c[0].first),s)); 
						node[1].first.complemented = 1; 
					}
					else if (flag[1] == 4)
					{
						auto s = best_neg - !outcompl; 
						node.push_back(std::make_pair(mig_create_maj(mig_new,!a[0].first,!b[1].first,!c[1].first),s)); 
						node[1].first.complemented = 1; 
					}
		
				}
				return node; 
			}

			std::vector<std::pair<mig_function,unsigned>> abc_double_last (mig_graph& mig_new, std::vector<std::pair<mig_function,unsigned>> a, std::vector<std::pair<mig_function,unsigned>> b, std::vector<std::pair<mig_function,unsigned>> c, bool outcompl)
			{
				std::vector<std::pair<mig_function,unsigned>> node; 
	
				auto sum_start = sum_ce(a[0].first,b[0].first,c[0].first,outcompl) + a[0].second + b[0].second + c[0].second + outcompl; // caso 1
				auto sum_final = sum_ce_changed(a[0].first,b[0].first,c[0].first,outcompl) + a[0].second + b[0].second + c[0].second + !outcompl; 
				auto sum_double_a = sum_ce(a[1].first,b[0].first,c[0].first, outcompl) + a[1].second + b[0].second + c[0].second+ outcompl; // caso 2
				auto sum_change_a = sum_ce_changed(a[1].first,b[0].first,c[0].first,outcompl) + a[1].second + b[0].second + c[0].second+ !outcompl; 
				auto sum_double_b = sum_ce(a[0].first,b[1].first,c[0].first, outcompl)+ a[0].second + b[1].second + c[0].second+ outcompl; //caso 3
				auto sum_change_b = sum_ce_changed(a[0].first,b[1].first,c[0].first,outcompl) + a[0].second + b[1].second + c[0].second+ !outcompl; 
				auto sum_double_bc = sum_ce(a[0].first,b[1].first,c[1].first, outcompl)+ a[0].second + b[1].second + c[1].second+ outcompl; // caso 4
				auto sum_change_bc = sum_ce_changed(a[0].first,b[1].first,c[1].first,outcompl)+ a[0].second + b[1].second + c[1].second+ !outcompl; 
				auto sum_double_c = sum_ce(a[0].first,b[0].first,c[1].first, outcompl)+ a[0].second + b[0].second + c[1].second+ outcompl; // caso 5
				auto sum_change_c = sum_ce_changed(a[0].first,b[0].first,c[1].first,outcompl)+ a[0].second + b[0].second + c[1].second+ !outcompl; 
				auto sum_double_ac = sum_ce(a[1].first,b[0].first,c[1].first, outcompl)+ a[1].second + b[0].second + c[1].second+ outcompl; // caso 6
				auto sum_change_ac = sum_ce_changed(a[1].first,b[0].first,c[1].first,outcompl)+ a[1].second + b[0].second + c[1].second+ !outcompl; 
				auto sum_double_ab = sum_ce(a[1].first,b[1].first,c[0].first, outcompl)+ a[1].second + b[1].second + c[0].second+ outcompl; // caso 7
				auto sum_change_ab = sum_ce_changed(a[1].first,b[1].first,c[0].first,outcompl)+ a[1].second + b[1].second + c[0].second+ !outcompl; 
				auto sum_double_abc = sum_ce(a[1].first,b[1].first,c[1].first, outcompl)+ a[1].second + b[1].second + c[1].second+ outcompl; // caso 8
				auto sum_change_abc = sum_ce_changed(a[1].first,b[1].first,c[1].first,outcompl)+ a[1].second + b[1].second + c[1].second+ !outcompl; 
	
				auto best_pos = sum_start;
				auto best_neg = sum_final; 
				std::vector<unsigned> flag(2);  
				flag[0] = 1; 
				flag[1] = 1;   
				if (sum_double_a < best_pos)
				{
					best_pos = sum_double_a;
					flag[0] = 2; 
				}
				if (sum_double_b < best_pos)
				{
					best_pos = sum_double_b;
					flag[0] = 3; 
				}
				if (sum_double_ab < best_pos)
				{
					best_pos = sum_double_ab;
					flag[0] = 7; 
				}
				if (sum_double_bc < best_pos)
				{
					best_pos = sum_double_bc;
					flag[0] = 4; 
				}
				if (sum_double_ac < best_pos)
				{
					best_pos = sum_double_ac;
					flag[0] = 6; 
				}
				if (sum_double_c < best_pos)
				{
					best_pos = sum_double_ac;
					flag[0] = 5; 
				}
				if (sum_double_abc < best_pos)
				{
					best_pos = sum_double_ac;
					flag[0] = 8; 
				}
	
				if (sum_change_a < best_neg)
				{
					best_neg = sum_change_a;
					flag[1] = 2; 
				}
				if (sum_change_b < best_neg)
				{
					best_neg = sum_change_b;
					flag[1] = 3; 
				}
				if (sum_change_ab < best_neg)
				{
					best_neg = sum_change_ab;
					flag[1] = 7; 
				}
				if (sum_change_ac < best_neg)
				{
					best_neg = sum_change_ac;
					flag[1] = 6; 
				}
				if (sum_change_bc < best_neg)
				{
					best_neg = sum_change_bc;
					flag[1] = 4; 
				}
				if (sum_change_c < best_neg)
				{
					best_neg = sum_change_c;
					flag[1] = 5; 
				}
				if (sum_change_abc < best_neg)
				{
					best_neg = sum_change_abc;
					flag[1] = 8; 
				}
				
				if (best_pos <= best_neg)
				{
					if ( flag[0] == 1)
					{
						auto s = best_pos;
						node.push_back(std::make_pair(mig_create_maj(mig_new,a[0].first,b[0].first,c[0].first),s)); 
					}
		
					else if (flag[0] == 2)
					{
						auto s = best_pos ;
						node.push_back(std::make_pair(mig_create_maj(mig_new,a[1].first,b[0].first,c[0].first),s)); 
					}
					else if (flag[0] == 3)
					{
						auto s = best_pos;
						node.push_back(std::make_pair(mig_create_maj(mig_new,a[0].first,b[1].first,c[0].first),s)); 
					} 
					else if (flag[0] == 4)
					{
						auto s = best_pos ;
						node.push_back(std::make_pair(mig_create_maj(mig_new,a[0].first,b[1].first,c[1].first),s)); 
					} 
					if ( flag[0] == 5)
					{
						auto s = best_pos ;
						node.push_back(std::make_pair(mig_create_maj(mig_new,a[0].first,b[0].first,c[1].first),s)); 
					}
		
					else if (flag[0] == 6)
					{
						auto s = best_pos ;
						node.push_back(std::make_pair(mig_create_maj(mig_new,a[1].first,b[0].first,c[1].first),s)); 
					}
					else if (flag[0] == 7)
					{
						auto s = best_pos ;
						node.push_back(std::make_pair(mig_create_maj(mig_new,a[1].first,b[1].first,c[0].first),s)); 
					} 
					else if (flag[0] == 8)
					{
						auto s = best_pos;
						node.push_back(std::make_pair(mig_create_maj(mig_new,a[1].first,b[1].first,c[1].first),s)); 
					} 
				}
				else if (best_pos >best_neg)
				{
					if ( flag[1] == 1)
					{
						auto s = best_neg ;
						node.push_back(std::make_pair(mig_create_maj(mig_new,!a[0].first,!b[0].first,!c[0].first),s)); 
						node[1].first.complemented = 1; 
					}
					else if (flag[1] == 2)
					{
						auto s = best_neg ;
						node.push_back(std::make_pair(mig_create_maj(mig_new,!a[1].first,!b[0].first,!c[0].first),s)); 
						node[1].first.complemented = 1; 
					}
					else if ( flag[1] == 3)
					{
						auto s = best_neg ;
						node.push_back(std::make_pair(mig_create_maj(mig_new,!a[0].first,!b[1].first,!c[0].first),s)); 
						node[1].first.complemented = 1; 
					}
					else if (flag[1] == 4)
					{
						auto s = best_neg ;
						node.push_back(std::make_pair(mig_create_maj(mig_new,!a[0].first,!b[1].first,!c[1].first),s)); 
						node[1].first.complemented = 1; 
					}
					if ( flag[1] == 5)
					{
						auto s = best_neg ;
						node.push_back(std::make_pair(mig_create_maj(mig_new,!a[0].first,!b[0].first,!c[1].first),s)); 
						node[1].first.complemented = 1; 
					}
					else if (flag[1] == 6)
					{
						auto s = best_neg ;
						node.push_back(std::make_pair(mig_create_maj(mig_new,!a[1].first,!b[0].first,!c[1].first),s)); 
						node[1].first.complemented = 1; 
					}
					else if ( flag[1] == 7)
					{
						auto s = best_neg;
						node.push_back(std::make_pair(mig_create_maj(mig_new,!a[1].first,!b[1].first,!c[0].first),s)); 
						node[1].first.complemented = 1; 
					}
					else if (flag[1] == 8)
					{
						auto s = best_neg ;
						node.push_back(std::make_pair(mig_create_maj(mig_new,!a[1].first,!b[1].first,!c[1].first),s)); 
						node[1].first.complemented = 1; 
					}
	
				}
						
				return node; 
			
			}


				

std::vector<std::pair<mig_function,unsigned>> a_double (mig_graph& mig_new, std::vector<std::pair<mig_function,unsigned>> a, std::vector<std::pair<mig_function,unsigned>> b, std::vector<std::pair<mig_function,unsigned>> c, bool outcompl)
{
	std::vector<std::pair<mig_function,unsigned>> node; 
	
	auto sum_start = sum_ce(a[0].first,b[0].first,c[0].first,outcompl) + a[0].second + b[0].second + c[0].second; 
	auto sum_final = sum_ce_changed (a[0].first,b[0].first,c[0].first,outcompl) + a[0].second + b[0].second + c[0].second; 
	auto sum_double = sum_ce(a[1].first,b[0].first,c[0].first,outcompl) + a[1].second + b[0].second + c[0].second; 
	auto sum_change_double = sum_ce_changed (a[1].first,b[0].first,c[0].first,outcompl) + a[1].second + b[0].second + c[0].second; 

	
	auto best_pos = sum_start;
	auto best_neg = sum_final; 
	std::vector<unsigned> flag(2);    
	if (sum_start < sum_double)
	{
		best_pos = sum_start;
		flag[0] = 1; 
	}
	else 
	{ 
		best_pos = sum_double;
		flag[0] = 2;  
	}
	
	if (sum_final < sum_change_double)
	{
		best_neg = sum_final; 
		flag[1] = 1; 
	}
	else 
	{ 
		best_neg = sum_change_double;
		flag[1] = 2; 
	}
	
	//if ( best_pos < best_neg)
	//{
	//	if ( flag[0] == 1)
	//	{
	//		auto s = best_pos - outcompl; 
	//		node.push_back(std::make_pair(mig_create_maj(mig_new,a[0].first,b[0].first,c[0].first),s)); 
	//	}
	//	if (flag[0] == 2)
	//	{
	//		auto s = best_pos - outcompl; 
	//		node.push_back(std::make_pair(mig_create_maj(mig_new,a[1].first,b[0].first,c[0].first),s)); 
		//} 
	//}
	//else if (best_pos > best_neg)
	//{
	//	if ( flag[1] == 1)
	//	{
	//		auto s = best_neg - !outcompl; 
	//		node.push_back(std::make_pair(mig_create_maj(mig_new,!a[0].first,!b[0].first,!c[0].first),s)); 
	//		node[0].first.complemented = 1; 
	//	}
	//	if (flag[1] == 2)
	//	{
	//		auto s = best_neg - !outcompl; 
	//		node.push_back(std::make_pair(mig_create_maj(mig_new,!a[1].first,!b[0].first,!c[0].first),s)); 
	//		node[0].first.complemented = 1; 
	//	}
	//} 
	//else if ( best_pos == best_neg)
	//{
		if ( flag[0] == 1)
		{
			auto s = best_pos;
			node.push_back(std::make_pair(mig_create_maj(mig_new,a[0].first,b[0].first,c[0].first),s)); 
		} 
		else if (flag[0] == 2)
		{
			auto s = best_pos;
			node.push_back(std::make_pair(mig_create_maj(mig_new,a[1].first,b[0].first,c[0].first),s)); 
		} 
		if ( flag[1] == 1)
		{
			auto s = best_neg;
			node.push_back(std::make_pair(mig_create_maj(mig_new,!a[0].first,!b[0].first,!c[0].first),s)); 
			node[1].first.complemented = 1; 
		}
		else if (flag[1] == 2)
		{
			auto s = best_neg ;
			node.push_back(std::make_pair(mig_create_maj(mig_new,!a[1].first,!b[0].first,!c[0].first),s)); 
			node[1].first.complemented = 1; 
		}
		
		//}
	return node; 
	
}

std::vector<std::pair<mig_function,unsigned>> b_double (mig_graph& mig_new, std::vector<std::pair<mig_function,unsigned>> a, std::vector<std::pair<mig_function,unsigned>> b, std::vector<std::pair<mig_function,unsigned>> c, bool outcompl)
{
	std::vector<std::pair<mig_function,unsigned>> node; 
	
	auto sum_start = sum_ce(a[0].first,b[0].first,c[0].first,outcompl) + a[0].second + b[0].second + c[0].second; 
	auto sum_final = sum_ce_changed (a[0].first,b[0].first,c[0].first,outcompl) + a[0].second + b[0].second + c[0].second; 
	auto sum_double = sum_ce(a[0].first,b[1].first,c[0].first,outcompl) + a[0].second + b[1].second + c[0].second; 
	auto sum_change_double = sum_ce_changed (a[0].first,b[1].first,c[0].first,outcompl) + a[0].second + b[1].second + c[0].second; 

	
	auto best_pos = sum_start;
	auto best_neg = sum_final; 
	std::vector<unsigned> flag(2);    
	if (sum_start < sum_double)
	{
		best_pos = sum_start;
		flag[0] = 1; 
	}
	else 
	{ 
		best_pos = sum_double;
		flag[0] = 2;  
	}
	
	if (sum_final < sum_change_double)
	{
		best_neg = sum_final; 
		flag[1] = 1; 
	}
	else 
	{ 
		best_neg = sum_change_double;
		flag[1] = 2; 
	}
	
		if ( flag[0] == 1)
		{
			auto s = best_pos;
			node.push_back(std::make_pair(mig_create_maj(mig_new,a[0].first,b[0].first,c[0].first),s)); 
		} 
		else if (flag[0] == 2)
		{
			auto s = best_pos;
			node.push_back(std::make_pair(mig_create_maj(mig_new,a[0].first,b[1].first,c[0].first),s)); 
		} 
		if ( flag[1] == 1)
		{
			auto s = best_neg ;
			node.push_back(std::make_pair(mig_create_maj(mig_new,!a[0].first,!b[0].first,!c[0].first),s)); 
			node[1].first.complemented = 1; 
		}
		else if (flag[1] == 2)
		{
			auto s = best_neg;
			node.push_back(std::make_pair(mig_create_maj(mig_new,!a[0].first,!b[1].first,!c[0].first),s)); 
			node[1].first.complemented = 1; 
		}
		
	return node; 
	
}

std::vector<std::pair<mig_function,unsigned>> c_double (mig_graph& mig_new, std::vector<std::pair<mig_function,unsigned>> a, std::vector<std::pair<mig_function,unsigned>> b, std::vector<std::pair<mig_function,unsigned>> c, bool outcompl)
{
	std::vector<std::pair<mig_function,unsigned>> node; 
	
	auto sum_start = sum_ce(a[0].first,b[0].first,c[0].first,outcompl) + a[0].second + b[0].second + c[0].second; 
	auto sum_final = sum_ce_changed (a[0].first,b[0].first,c[0].first,outcompl) + a[0].second + b[0].second + c[0].second; 
	auto sum_double = sum_ce(a[0].first,b[0].first,c[1].first,outcompl) + a[0].second + b[0].second + c[1].second; 
	auto sum_change_double = sum_ce_changed (a[0].first,b[0].first,c[1].first,outcompl) + a[0].second + b[0].second + c[1].second; 

	
	auto best_pos = sum_start;
	auto best_neg = sum_final; 
	std::vector<unsigned> flag(2);    
	if (sum_start < sum_double)
	{
		best_pos = sum_start;
		flag[0] = 1; 
	}
	else 
	{ 
		best_pos = sum_double;
		flag[0] = 2;  
	}
	
	if (sum_final < sum_change_double)
	{
		best_neg = sum_final; 
		flag[1] = 1; 
	}
	else 
	{ 
		best_neg = sum_change_double;
		flag[1] = 2; 
	}
	
		if ( flag[0] == 1)
		{
			auto s = best_pos;
			node.push_back(std::make_pair(mig_create_maj(mig_new,a[0].first,b[0].first,c[0].first),s)); 
		} 
		else if (flag[0] == 2)
		{
			auto s = best_pos;
			node.push_back(std::make_pair(mig_create_maj(mig_new,a[0].first,b[0].first,c[1].first),s)); 
		} 
		if ( flag[1] == 1)
		{
			auto s = best_neg;
			node.push_back(std::make_pair(mig_create_maj(mig_new,!a[0].first,!b[0].first,!c[0].first),s)); 
			node[1].first.complemented = 1; 
		}
		else if (flag[1] == 2)
		{
			auto s = best_neg;
			node.push_back(std::make_pair(mig_create_maj(mig_new,!a[0].first,!b[0].first,!c[1].first),s)); 
			node[1].first.complemented = 1; 
		}
		
	
	return node; 
}

std::vector<std::pair<mig_function,unsigned>> ab_double (mig_graph& mig_new, std::vector<std::pair<mig_function,unsigned>> a, std::vector<std::pair<mig_function,unsigned>> b, std::vector<std::pair<mig_function,unsigned>> c, bool outcompl)
{
	std::vector<std::pair<mig_function,unsigned>> node; 
	
	auto sum_start = sum_ce(a[0].first,b[0].first,c[0].first,outcompl) + a[0].second + b[0].second + c[0].second; // caso 1
	auto sum_final = sum_ce_changed(a[0].first,b[0].first,c[0].first,outcompl) + a[0].second + b[0].second + c[0].second; 
	auto sum_double_a = sum_ce(a[1].first,b[0].first,c[0].first, outcompl) + a[1].second + b[0].second + c[0].second; // caso 2
	auto sum_change_a = sum_ce_changed(a[1].first,b[0].first,c[0].first,outcompl) + a[1].second + b[0].second + c[0].second; 
	auto sum_double_b = sum_ce(a[0].first,b[1].first,c[0].first, outcompl)+ a[0].second + b[1].second + c[0].second; //caso 3
	auto sum_change_b = sum_ce_changed(a[0].first,b[1].first,c[0].first,outcompl) + a[0].second + b[1].second + c[0].second; 
	auto sum_double_ab = sum_ce(a[1].first,b[1].first,c[0].first, outcompl)+ a[1].second + b[1].second + c[0].second; // caso 4
	auto sum_change_ab = sum_ce_changed(a[1].first,b[1].first,c[0].first,outcompl)+ a[1].second + b[1].second + c[0].second; 
	
	auto best_pos = sum_start;
	auto best_neg = sum_final; 
	std::vector<unsigned> flag(2);  
	flag[0] = 1; 
	flag[1] = 1;   
	if (sum_double_a < best_pos)
	{
		best_pos = sum_double_a;
		flag[0] = 2; 
	}
	if (sum_double_b < best_pos)
	{
		best_pos = sum_double_b;
		flag[0] = 3; 
	}
	if (sum_double_ab < best_pos)
	{
		best_pos = sum_double_ab;
		flag[0] = 4; 
	}
	
	if (sum_change_a < best_neg)
	{
		best_neg = sum_change_a;
		flag[1] = 2; 
	}
	if (sum_change_b < best_neg)
	{
		best_neg = sum_change_b;
		flag[1] = 3; 
	}
	if (sum_change_ab < best_neg)
	{
		best_neg = sum_change_ab;
		flag[1] = 4; 
	}
	
		if ( flag[0] == 1)
		{
			auto s = best_pos;
			node.push_back(std::make_pair(mig_create_maj(mig_new,a[0].first,b[0].first,c[0].first),s)); 
		}
		
		else if (flag[0] == 2)
		{
			auto s = best_pos ;
			node.push_back(std::make_pair(mig_create_maj(mig_new,a[1].first,b[0].first,c[0].first),s)); 
		}
		else if (flag[0] == 3)
		{
			auto s = best_pos;
			node.push_back(std::make_pair(mig_create_maj(mig_new,a[0].first,b[1].first,c[0].first),s)); 
		} 
		else if (flag[0] == 4)
		{
			auto s = best_pos ;
			node.push_back(std::make_pair(mig_create_maj(mig_new,a[1].first,b[1].first,c[0].first),s)); 
		} 
		if ( flag[1] == 1)
		{
			auto s = best_neg;
			node.push_back(std::make_pair(mig_create_maj(mig_new,!a[0].first,!b[0].first,!c[0].first),s)); 
			node[1].first.complemented = 1; 
		}
		else if (flag[1] == 2)
		{
			auto s = best_neg ;
			node.push_back(std::make_pair(mig_create_maj(mig_new,!a[1].first,!b[0].first,!c[0].first),s)); 
			node[1].first.complemented = 1; 
		}
		else if ( flag[1] == 3)
		{
			auto s = best_neg ;
			node.push_back(std::make_pair(mig_create_maj(mig_new,!a[0].first,!b[1].first,!c[0].first),s)); 
			node[1].first.complemented = 1; 
		}
		else if (flag[1] == 4)
		{
			auto s = best_neg ;
			node.push_back(std::make_pair(mig_create_maj(mig_new,!a[1].first,!b[1].first,!c[0].first),s)); 
			node[1].first.complemented = 1; 
		}
		
	
	return node; 
}

std::vector<std::pair<mig_function,unsigned>> ac_double (mig_graph& mig_new, std::vector<std::pair<mig_function,unsigned>> a, std::vector<std::pair<mig_function,unsigned>> b, std::vector<std::pair<mig_function,unsigned>> c, bool outcompl)
{
	std::vector<std::pair<mig_function,unsigned>> node; 
	
	auto sum_start = sum_ce(a[0].first,b[0].first,c[0].first,outcompl) + a[0].second + b[0].second + c[0].second; // caso 1
	auto sum_final = sum_ce_changed(a[0].first,b[0].first,c[0].first,outcompl) + a[0].second + b[0].second + c[0].second; 
	auto sum_double_a = sum_ce(a[1].first,b[0].first,c[0].first, outcompl) + a[1].second + b[0].second + c[0].second; // caso 2
	auto sum_change_a = sum_ce_changed(a[1].first,b[0].first,c[0].first,outcompl) + a[1].second + b[0].second + c[0].second; 
	auto sum_double_b = sum_ce(a[0].first,b[0].first,c[1].first, outcompl)+ a[0].second + b[0].second + c[1].second; //caso 3
	auto sum_change_b = sum_ce_changed(a[0].first,b[0].first,c[1].first,outcompl) + a[0].second + b[0].second + c[1].second; 
	auto sum_double_ab = sum_ce(a[1].first,b[0].first,c[1].first, outcompl)+ a[1].second + b[0].second + c[1].second; // caso 4
	auto sum_change_ab = sum_ce_changed(a[1].first,b[0].first,c[1].first,outcompl)+ a[1].second + b[0].second + c[1].second; 
	
	auto best_pos = sum_start;
	auto best_neg = sum_final; 
	std::vector<unsigned> flag(2);  
	flag[0] = 1; 
	flag[1] = 1;   
	if (sum_double_a < best_pos)
	{
		best_pos = sum_double_a;
		flag[0] = 2; 
	}
	if (sum_double_b < best_pos)
	{
		best_pos = sum_double_b;
		flag[0] = 3; 
	}
	if (sum_double_ab < best_pos)
	{
		best_pos = sum_double_ab;
		flag[0] = 4; 
	}
	
	if (sum_change_a < best_neg)
	{
		best_neg = sum_change_a;
		flag[1] = 2; 
	}
	if (sum_change_b < best_neg)
	{
		best_neg = sum_change_b;
		flag[1] = 3; 
	}
	if (sum_change_ab < best_neg)
	{
		best_neg = sum_change_ab;
		flag[1] = 4; 
	}
	
		if ( flag[0] == 1)
		{
			auto s = best_pos ;
			node.push_back(std::make_pair(mig_create_maj(mig_new,a[0].first,b[0].first,c[0].first),s)); 
		}
		
		else if (flag[0] == 2)
		{
			auto s = best_pos ;
			node.push_back(std::make_pair(mig_create_maj(mig_new,a[1].first,b[0].first,c[0].first),s)); 
		}
		else if (flag[0] == 3)
		{
			auto s = best_pos ;
			node.push_back(std::make_pair(mig_create_maj(mig_new,a[0].first,b[0].first,c[1].first),s)); 
		} 
		else if (flag[0] == 4)
		{
			auto s = best_pos ;
			node.push_back(std::make_pair(mig_create_maj(mig_new,a[1].first,b[0].first,c[1].first),s)); 
		} 
		if ( flag[1] == 1)
		{
			auto s = best_neg ;
			node.push_back(std::make_pair(mig_create_maj(mig_new,!a[0].first,!b[0].first,!c[0].first),s)); 
			node[1].first.complemented = 1; 
		}
		else if (flag[1] == 2)
		{
			auto s = best_neg; 
			node.push_back(std::make_pair(mig_create_maj(mig_new,!a[1].first,!b[0].first,!c[0].first),s)); 
			node[1].first.complemented = 1; 
		}
		else if ( flag[1] == 3)
		{
			auto s = best_neg ;
			node.push_back(std::make_pair(mig_create_maj(mig_new,!a[0].first,!b[0].first,!c[1].first),s)); 
			node[1].first.complemented = 1; 
		}
		else if (flag[1] == 4)
		{
			auto s = best_neg ; 
			node.push_back(std::make_pair(mig_create_maj(mig_new,!a[1].first,!b[0].first,!c[1].first),s)); 
			node[1].first.complemented = 1; 
		}
		
	
	return node; 
}

std::vector<std::pair<mig_function,unsigned>> bc_double (mig_graph& mig_new, std::vector<std::pair<mig_function,unsigned>> a, std::vector<std::pair<mig_function,unsigned>> b, std::vector<std::pair<mig_function,unsigned>> c, bool outcompl)
{
	std::vector<std::pair<mig_function,unsigned>> node; 
	
	//std::cout << " siamo dentro bc " << std::endl; 
	
	auto sum_start = sum_ce(a[0].first,b[0].first,c[0].first,outcompl) + a[0].second + b[0].second + c[0].second; // caso 1
	auto sum_final = sum_ce_changed(a[0].first,b[0].first,c[0].first,outcompl) + a[0].second + b[0].second + c[0].second; 
	auto sum_double_a = sum_ce(a[0].first,b[0].first,c[1].first, outcompl) + a[0].second + b[0].second + c[1].second; // caso 2
	auto sum_change_a = sum_ce_changed(a[0].first,b[0].first,c[1].first,outcompl) + a[0].second + b[0].second + c[1].second; 
	auto sum_double_b = sum_ce(a[0].first,b[1].first,c[0].first, outcompl)+ a[0].second + b[1].second + c[0].second; //caso 3
	auto sum_change_b = sum_ce_changed(a[0].first,b[1].first,c[0].first,outcompl) + a[0].second + b[1].second + c[0].second; 
	auto sum_double_ab = sum_ce(a[0].first,b[1].first,c[1].first, outcompl)+ a[0].second + b[1].second + c[1].second; // caso 4
	auto sum_change_ab = sum_ce_changed(a[0].first,b[1].first,c[1].first,outcompl)+ a[0].second + b[1].second + c[1].second; 
	
	auto best_pos = sum_start;
	auto best_neg = sum_final; 
	std::vector<unsigned> flag(2);  
	flag[0] = 1; 
	flag[1] = 1;   
	if (sum_double_a < best_pos)
	{
		best_pos = sum_double_a;
		flag[0] = 2; 
	}
	if (sum_double_b < best_pos)
	{
		best_pos = sum_double_b;
		flag[0] = 3; 
	}
	if (sum_double_ab < best_pos)
	{
		best_pos = sum_double_ab;
		flag[0] = 4; 
	}
	
	if (sum_change_a < best_neg)
	{
		best_neg = sum_change_a;
		flag[1] = 2; 
	}
	if (sum_change_b < best_neg)
	{
		best_neg = sum_change_b;
		flag[1] = 3; 
	}
	if (sum_change_ab < best_neg)
	{
		best_neg = sum_change_ab;
		flag[1] = 4; 
	}
	
	//std::cout << " best pos = " << best_pos << std::endl;
	//std::cout << " best neg = " << best_neg << std::endl;
		if ( flag[0] == 1)
		{
			auto s = best_pos ;
			node.push_back(std::make_pair(mig_create_maj(mig_new,a[0].first,b[0].first,c[0].first),s)); 
		}
		
		else if (flag[0] == 2)
		{
			auto s = best_pos;
			node.push_back(std::make_pair(mig_create_maj(mig_new,a[0].first,b[0].first,c[1].first),s)); 
		}
		else if (flag[0] == 3)
		{
			auto s = best_pos ;
			node.push_back(std::make_pair(mig_create_maj(mig_new,a[0].first,b[1].first,c[0].first),s)); 
		} 
		else if (flag[0] == 4)
		{
			auto s = best_pos;
			node.push_back(std::make_pair(mig_create_maj(mig_new,a[0].first,b[1].first,c[1].first),s)); 
		} 
		if ( flag[1] == 1)
		{
			auto s = best_neg ;
			node.push_back(std::make_pair(mig_create_maj(mig_new,!a[0].first,!b[0].first,!c[0].first),s)); 
			node[1].first.complemented = 1; 
		}
		else if (flag[1] == 2)
		{
			auto s = best_neg ;
			node.push_back(std::make_pair(mig_create_maj(mig_new,!a[0].first,!b[0].first,!c[1].first),s)); 
			node[1].first.complemented = 1; 
		}
		else if ( flag[1] == 3)
		{
			auto s = best_neg ;
			node.push_back(std::make_pair(mig_create_maj(mig_new,!a[0].first,!b[1].first,!c[0].first),s)); 
			node[1].first.complemented = 1; 
		}
		else if (flag[1] == 4)
		{
			auto s = best_neg ;
			node.push_back(std::make_pair(mig_create_maj(mig_new,!a[0].first,!b[1].first,!c[1].first),s)); 
			node[1].first.complemented = 1; 
		}
		
	
	return node; 
}

std::vector<std::pair<mig_function,unsigned>> abc_double (mig_graph& mig_new, std::vector<std::pair<mig_function,unsigned>> a, std::vector<std::pair<mig_function,unsigned>> b, std::vector<std::pair<mig_function,unsigned>> c, bool outcompl)
{
	std::vector<std::pair<mig_function,unsigned>> node; 
	
	auto sum_start = sum_ce(a[0].first,b[0].first,c[0].first,outcompl) + a[0].second + b[0].second + c[0].second; // caso 1
	auto sum_final = sum_ce_changed(a[0].first,b[0].first,c[0].first,outcompl) + a[0].second + b[0].second + c[0].second; 
	auto sum_double_a = sum_ce(a[1].first,b[0].first,c[0].first, outcompl) + a[0].second + b[0].second + c[0].second; // caso 2
	auto sum_change_a = sum_ce_changed(a[1].first,b[0].first,c[0].first,outcompl) + a[1].second + b[0].second + c[0].second; 
	auto sum_double_b = sum_ce(a[0].first,b[1].first,c[0].first, outcompl)+ a[0].second + b[1].second + c[0].second; //caso 3
	auto sum_change_b = sum_ce_changed(a[0].first,b[1].first,c[0].first,outcompl) + a[0].second + b[1].second + c[0].second; 
	auto sum_double_bc = sum_ce(a[0].first,b[1].first,c[1].first, outcompl)+ a[1].second + b[0].second + c[1].second; // caso 4
	auto sum_change_bc = sum_ce_changed(a[0].first,b[1].first,c[1].first,outcompl)+ a[0].second + b[1].second + c[1].second; 
	auto sum_double_c = sum_ce(a[0].first,b[0].first,c[1].first, outcompl)+ a[0].second + b[0].second + c[1].second; // caso 5
	auto sum_change_c = sum_ce_changed(a[0].first,b[0].first,c[1].first,outcompl)+ a[0].second + b[0].second + c[1].second; 
	auto sum_double_ac = sum_ce(a[1].first,b[0].first,c[1].first, outcompl)+ a[1].second + b[0].second + c[1].second; // caso 6
	auto sum_change_ac = sum_ce_changed(a[1].first,b[0].first,c[1].first,outcompl)+ a[1].second + b[0].second + c[1].second; 
	auto sum_double_ab = sum_ce(a[1].first,b[1].first,c[0].first, outcompl)+ a[1].second + b[1].second + c[0].second; // caso 7
	auto sum_change_ab = sum_ce_changed(a[1].first,b[1].first,c[0].first,outcompl)+ a[1].second + b[1].second + c[0].second; 
	auto sum_double_abc = sum_ce(a[1].first,b[1].first,c[1].first, outcompl)+ a[1].second + b[1].second + c[1].second; // caso 8
	auto sum_change_abc = sum_ce_changed(a[1].first,b[1].first,c[1].first,outcompl)+ a[1].second + b[1].second + c[1].second; 
	
	auto best_pos = sum_start;
	auto best_neg = sum_final; 
	std::vector<unsigned> flag(2);  
	flag[0] = 1; 
	flag[1] = 1;   
	if (sum_double_a < best_pos)
	{
		best_pos = sum_double_a;
		flag[0] = 2; 
	}
	if (sum_double_b < best_pos)
	{
		best_pos = sum_double_b;
		flag[0] = 3; 
	}
	if (sum_double_ab < best_pos)
	{
		best_pos = sum_double_ab;
		flag[0] = 7; 
	}
	if (sum_double_bc < best_pos)
	{
		best_pos = sum_double_bc;
		flag[0] = 4; 
	}
	if (sum_double_ac < best_pos)
	{
		best_pos = sum_double_ac;
		flag[0] = 6; 
	}
	if (sum_double_c < best_pos)
	{
		best_pos = sum_double_ac;
		flag[0] = 5; 
	}
	if (sum_double_abc < best_pos)
	{
		best_pos = sum_double_ac;
		flag[0] = 8; 
	}
	
	if (sum_change_a < best_neg)
	{
		best_neg = sum_change_a;
		flag[1] = 2; 
	}
	if (sum_change_b < best_neg)
	{
		best_neg = sum_change_b;
		flag[1] = 3; 
	}
	if (sum_change_ab < best_neg)
	{
		best_neg = sum_change_ab;
		flag[1] = 7; 
	}
	if (sum_change_ac < best_neg)
	{
		best_neg = sum_change_ac;
		flag[1] = 6; 
	}
	if (sum_change_bc < best_neg)
	{
		best_neg = sum_change_bc;
		flag[1] = 4; 
	}
	if (sum_change_c < best_neg)
	{
		best_neg = sum_change_c;
		flag[1] = 5; 
	}
	if (sum_change_abc < best_neg)
	{
		best_neg = sum_change_abc;
		flag[1] = 8; 
	}
	
		if ( flag[0] == 1)
		{
			auto s = best_pos;
			node.push_back(std::make_pair(mig_create_maj(mig_new,a[0].first,b[0].first,c[0].first),s)); 
		}
		
		else if (flag[0] == 2)
		{
			auto s = best_pos ;
			node.push_back(std::make_pair(mig_create_maj(mig_new,a[1].first,b[0].first,c[0].first),s)); 
		}
		else if (flag[0] == 3)
		{
			auto s = best_pos;
			node.push_back(std::make_pair(mig_create_maj(mig_new,a[0].first,b[1].first,c[0].first),s)); 
		} 
		else if (flag[0] == 4)
		{
			auto s = best_pos ;
			node.push_back(std::make_pair(mig_create_maj(mig_new,a[0].first,b[1].first,c[1].first),s)); 
		} 
		if ( flag[0] == 5)
		{
			auto s = best_pos ;
			node.push_back(std::make_pair(mig_create_maj(mig_new,a[0].first,b[0].first,c[1].first),s)); 
		}
		
		else if (flag[0] == 6)
		{
			auto s = best_pos ;
			node.push_back(std::make_pair(mig_create_maj(mig_new,a[1].first,b[0].first,c[1].first),s)); 
		}
		else if (flag[0] == 7)
		{
			auto s = best_pos ;
			node.push_back(std::make_pair(mig_create_maj(mig_new,a[1].first,b[1].first,c[0].first),s)); 
		} 
		else if (flag[0] == 8)
		{
			auto s = best_pos;
			node.push_back(std::make_pair(mig_create_maj(mig_new,a[1].first,b[1].first,c[1].first),s)); 
		} 
		
		if ( flag[1] == 1)
		{
			auto s = best_neg ;
			node.push_back(std::make_pair(mig_create_maj(mig_new,!a[0].first,!b[0].first,!c[0].first),s)); 
			node[1].first.complemented = 1; 
		}
		else if (flag[1] == 2)
		{
			auto s = best_neg ;
			node.push_back(std::make_pair(mig_create_maj(mig_new,!a[1].first,!b[0].first,!c[0].first),s)); 
			node[1].first.complemented = 1; 
		}
		else if ( flag[1] == 3)
		{
			auto s = best_neg ;
			node.push_back(std::make_pair(mig_create_maj(mig_new,!a[0].first,!b[1].first,!c[0].first),s)); 
			node[1].first.complemented = 1; 
		}
		else if (flag[1] == 4)
		{
			auto s = best_neg ;
			node.push_back(std::make_pair(mig_create_maj(mig_new,!a[0].first,!b[1].first,!c[1].first),s)); 
			node[1].first.complemented = 1; 
		}
		if ( flag[1] == 5)
		{
			auto s = best_neg ;
			node.push_back(std::make_pair(mig_create_maj(mig_new,!a[0].first,!b[0].first,!c[1].first),s)); 
			node[1].first.complemented = 1; 
		}
		else if (flag[1] == 6)
		{
			auto s = best_neg ;
			node.push_back(std::make_pair(mig_create_maj(mig_new,!a[1].first,!b[0].first,!c[1].first),s)); 
			node[1].first.complemented = 1; 
		}
		else if ( flag[1] == 7)
		{
			auto s = best_neg;
			node.push_back(std::make_pair(mig_create_maj(mig_new,!a[1].first,!b[1].first,!c[0].first),s)); 
			node[1].first.complemented = 1; 
		}
		else if (flag[1] == 8)
		{
			auto s = best_neg ;
			node.push_back(std::make_pair(mig_create_maj(mig_new,!a[1].first,!b[1].first,!c[1].first),s)); 
			node[1].first.complemented = 1; 
		}
	
	return node; 
}

std::vector<std::pair<mig_function,unsigned>>  change_ce (mig_graph& mig_new, std::vector<std::pair<mig_function,unsigned>> a, std::vector<std::pair<mig_function,unsigned>> b, std::vector<std::pair<mig_function,unsigned>> c, bool outcompl, mig_graph mig, mig_node last)
{
	std::vector<std::pair<mig_function,unsigned>> node; 
	
	//std::cout << " il valore del nodo e " << last << std::endl; 
	//std::cout << " il valore del output e " << mig_info( mig ).outputs[0].first.node << std::endl; 
	//std::cout << " la size di a é " << a.size() << std::endl; 
	//std::cout << " la size di b é " << b.size() << std::endl; 
	//std::cout << " la size di c é " << c.size() << std::endl; 
	if ((e_output (mig,last) == false)  && (number_indegree(last,mig) == 1))

	//if ((last < mig_info( mig ).outputs[0].first.node) && (number_indegree(last,mig) == 1))
	{
	if ((a.size() == 1) && (b.size() == 1) && (c.size() == 1))
	{
	auto sum_start = sum_ce(a[0].first,b[0].first,c[0].first,outcompl) + a[0].second + b[0].second + c[0].second; 
	auto sum_final = sum_ce_changed (a[0].first,b[0].first,c[0].first,outcompl) + a[0].second + b[0].second + c[0].second; 
	
//	std::cout << " sum_start = " << sum_start << std::endl; 
	//std::cout << " sum_final = " << sum_final << std::endl; 
	//if (sum_start > sum_final) 
	//{		
		//auto s = sum_final - !outcompl;
	//node.push_back(std::make_pair(mig_create_maj(mig_new,!a[0].first,!b[0].first,!c[0].first),s)); 
	//node[0].first.complemented = 1; 
	//}
	//else if (sum_start < sum_final)
	//{
	//	auto s = sum_start - outcompl;
	//node.push_back(std::make_pair(mig_create_maj(mig_new,a[0].first,b[0].first,c[0].first),s));
	
	//}
	//else if (sum_start == sum_final)
	//{
		auto s = sum_start; //- outcompl;
	node.push_back(std::make_pair(mig_create_maj(mig_new,a[0].first,b[0].first,c[0].first),s));
 //else 
	
	s = sum_final; // - !outcompl;
		//{
	node.push_back(std::make_pair(mig_create_maj(mig_new,!a[0].first,!b[0].first,!c[0].first),s));
	node[1].first.complemented = 1; 
	//}
	//}
	
	return node; 
	}
	
	else if ((a.size() > 1) && (b.size() == 1) && (c.size() == 1))
				node = a_double (mig_new, a,b,c,outcompl); 
	else if ((a.size() == 1) && (b.size() > 1) && (c.size() == 1))
				node = b_double (mig_new,a,b,c,outcompl); 
	else if ((a.size() == 1) && (b.size() == 1) && (c.size() > 1))
				node = c_double (mig_new,a,b,c,outcompl); 
	else if ((a.size() > 1) && (b.size() > 1) && (c.size() == 1))
				node = ab_double (mig_new,a,b,c,outcompl); 
	else if ((a.size() > 1) && (b.size() == 1) && (c.size() > 1))
				node = ac_double (mig_new,a,b,c,outcompl); 
	else if ((a.size() == 1) && (b.size() > 1) && (c.size() > 1))
				node = bc_double (mig_new,a,b,c,outcompl); 
	else if ((a.size() > 1) && (b.size() > 1) && (c.size() > 1))
				node = abc_double (mig_new,a,b,c,outcompl); 
	//else 
	//{
	//	auto sum_start = sum_ce(a[0].first,b[0].first,c[0].first,outcompl) + a[0].second + b[0].second + c[0].second; 
	//auto s = sum_start - outcompl;
	//node.push_back(std::make_pair(mig_create_maj(mig_new,a[0].first,b[0].first,c[0].first),s));
	//}
	return node; 
}

else if ((e_output (mig,last) == true) && (number_indegree(last,mig) == 1))
{
	if ((a.size() == 1) && (b.size() == 1) && (c.size() == 1))
	{
	auto sum_start = sum_ce(a[0].first,b[0].first,c[0].first,outcompl) + a[0].second + b[0].second + c[0].second; 
	auto sum_final = sum_ce_changed (a[0].first,b[0].first,c[0].first,outcompl) + a[0].second + b[0].second + c[0].second; 
	auto s = sum_start; //- outcompl;
	if (sum_start <= sum_final) 
		node.push_back(std::make_pair(mig_create_maj(mig_new,a[0].first,b[0].first,c[0].first),s));
	else 
	{
		s = sum_final; 
	node.push_back(std::make_pair(mig_create_maj(mig_new,!a[0].first,!b[0].first,!c[0].first),s));
	node[0].first.complemented = 1;
}
//s = sum_final; // - !outcompl;
	//{
//node.push_back(std::make_pair(mig_create_maj(mig_new,!a[0].first,!b[0].first,!c[0].first),s));
//node[1].first.complemented = 1; 
}
	else if ((a.size() > 1) && (b.size() == 1) && (c.size() == 1))
				node = a_double_last (mig_new, a,b,c,outcompl); 
	else if ((a.size() == 1) && (b.size() > 1) && (c.size() == 1))
				node = b_double_last (mig_new,a,b,c,outcompl); 
	else if ((a.size() == 1) && (b.size() == 1) && (c.size() > 1))
				node = c_double_last (mig_new,a,b,c,outcompl); 
	else if ((a.size() > 1) && (b.size() > 1) && (c.size() == 1))
				node = ab_double_last (mig_new,a,b,c,outcompl); 
	else if ((a.size() > 1) && (b.size() == 1) && (c.size() > 1))
				node = ac_double_last (mig_new,a,b,c,outcompl); 
	else if ((a.size() == 1) && (b.size() > 1) && (c.size() > 1))
				node = bc_double_last (mig_new,a,b,c,outcompl); 
	else if ((a.size() > 1) && (b.size() > 1) && (c.size() > 1))
				node = abc_double_last (mig_new,a,b,c,outcompl); 
}

else if ((e_output (mig,last) == false) && (number_indegree(last,mig) > 1))
	node.push_back(std::make_pair(mig_create_maj(mig_new,a[0].first,b[0].first,c[0].first),0));	
else if ((e_output (mig,last) == true) && (number_indegree(last,mig) > 1))
	node.push_back(std::make_pair(mig_create_maj(mig_new,a[0].first,b[0].first,c[0].first),0));

return node; 
}
	 

std::vector<std::pair<mig_function,unsigned>> mig_rewrite_ce_rec( const mig_graph& old_mig, mig_node node,
                                       mig_graph& new_mig,
                                       std::map<mig_node, std::vector<std::pair<mig_function,unsigned>>> & old_to_new, bool outcompl )
{
  /* visited */
	//bool outcompl = false; 
    const auto it = old_to_new.find( node );
    if ( it != old_to_new.end() )
    {
      return it->second;
    }

    const auto c = get_children( old_mig, node );
    const auto f = change_ce( new_mig,comple(
                           mig_rewrite_ce_rec( old_mig, c[0].node, new_mig, old_to_new, c[0].complemented ),c[0].complemented),
                           comple(mig_rewrite_ce_rec( old_mig, c[1].node, new_mig,old_to_new, c[1].complemented ),c[1].complemented),
                           comple(mig_rewrite_ce_rec( old_mig, c[2].node, new_mig,old_to_new, c[2].complemented ) ,c[2].complemented), outcompl, old_mig, node);

    old_to_new.insert( {node, f} );
    return f;
 

}

/******************************************************************************
 * Public functions                                                           *
 ******************************************************************************/

mig_graph ce_opt_tree( const mig_graph mig,
                                const properties::ptr& settings,
                                const properties::ptr& statistics )
{
  mig_graph new_mig; 
  mig_initialize( new_mig, mig_info( mig ).model_name );

  /* create constant and PIs */
  auto old_to_new = init_visited_table_ce( mig, new_mig);

  /* map nodes */
  for ( const auto& output : mig_info( mig ).outputs )
  {
      mig_create_po( new_mig, mig_rewrite_ce_rec( mig, output.first.node, new_mig, old_to_new ,output.first.complemented)[0].first ^ output.first.complemented, output.second );
  }
  
   new_mig = mig_strash(new_mig); 

  return new_mig;
}

}

// Local Variables:
// c-basic-offset: 2
// eval: (c-set-offset 'substatement-open 0)
// eval: (c-set-offset 'innamespace 0)
// End: