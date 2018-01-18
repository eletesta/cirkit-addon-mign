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

#include "mign.hpp"

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/topological_sort.hpp>
#include <boost/range/iterator_range.hpp>

#include <classical/mign/mign_utils.hpp> 
#include <classical/mign/mign_cover.hpp> 
#include <classical/mign/mign_bitmarks.hpp>

namespace cirkit
{

/******************************************************************************
 * Types                                                                      *
 ******************************************************************************/

mign_graph::mign_graph (/*const std::string& name, */const boost::optional<unsigned> nin)
	: constant(add_vertex(g)), 
//_name (name),
	_nin (nin),
	    _bitmarks( std::make_shared<mign_bitmarks>() )
{
	assert (constant == 0); 
	
}

void mign_graph::compute_fanout()
{
	fanout.update( [this]() { return precompute_in_degrees( g ); } );
}

void mign_graph::compute_parents()
{
  parentss.update( [this]() { return precompute_ingoing_vertices( g ); } );
}

void mign_graph::compute_levels()
{
  levels.update( [this]() { return mign_compute_levels( *this ); } );
}

mign_function mign_graph::get_constant(bool value) const
{
	return mign_function (constant, value); 
}

mign_function mign_graph::create_pi (const std::string& name)
{
	const auto node = add_vertex(g); 
	_input_to_id.insert ({node, _inputs.size()}); // from 0 
		/*we dont care about the compl since the node will store it*/ 
	_inputs.push_back({node,name}); 
	return mign_function (node); 
}

void mign_graph::create_po (const mign_function& f, const std::string name)
{
	_outputs.push_back({f,name}); // info about compl is in the mign_function
}

mign_function mign_graph::create_maj_fo_restricted (std::vector<mign_function> operands, const unsigned int fanout) // works for MAJORity 3 
{
	assert (operands.size() == 3); // By now it work for MAJ-3 for Odysseas project 
	
    if ( operands[0] == operands[1] )  { return operands[0]; }
    if ( operands[0] == operands[2] )  { return operands[0]; }
    if ( operands[1] == operands[2] )  { return operands[1]; }
	
	
    /* structural hashing */
	std::vector<mign_function> children; 
	for (auto i = 0; i<operands.size(); ++i)
	{
		children.push_back(operands[i]); 
	}
	
    std::sort(children.begin(),children.begin() + children.size());

    auto node_complement = false;
    
    auto key = children; 

    const auto it = maj_strash_fo_restricted.find( key );
    if ( it != maj_strash_fo_restricted.end() )
    {
		if ( it -> second.second < fanout)
		{
			maj_strash_fo_restricted.at(key).second++; 
			return mign_function( it->second.first, node_complement );
		}
		else 
		{
			//std::cout << " ERASER" << std::endl; 
			maj_strash_fo_restricted.erase(key); 
		}
    }

    /* insert node */
    const auto node = add_vertex( g );
    ++_num_maj;

    const auto ea = add_edge( node, children[0].node, g ).first;
    const auto eb = add_edge( node, children[1].node, g ).first;
    const auto ec = add_edge( node, children[2].node, g ).first;

    complement[ea] = children[0].complemented;
    complement[eb] = children[1].complemented;
    complement[ec] = children[2].complemented;

   
    maj_strash_fo_restricted[key].first = node;
	maj_strash_fo_restricted[key].second = 0;
	
    return mign_function( node, node_complement );
  }
  

mign_function mign_graph::create_maj_general (const std::vector<mign_function> operands)
{
	auto node_complement = false; //
	
	if (operands.size() == 1)
		return operands[0]; 
	// special cases not valid if nin is imposted 
	if ((_nin) && (*_nin != 3))
	{
		std::vector<mign_function> children; 
		for (auto i = 0; i<operands.size(); ++i)
		{
			children.push_back(operands[i]); 
		}
		std::sort(children.begin(),children.begin() + children.size()); // in order 
	
		auto key = children; 
		if (!(_disable))
		{
			const auto it = maj_strash.find(key); 
		    if ( it != maj_strash.end() )
		    {
				 return mign_function( it->second, node_complement );
		    }
		}
		
		const auto node = add_vertex(g); 
		++_num_maj; 
	
		for ( auto i =0 ;i<operands.size(); ++i)
		{
			const auto x = add_edge(node, children[i].node, g).first; 
			complement[x] = children[i].complemented; 
		}
	
		maj_strash[key] = node; 
		return mign_function(node, node_complement); 
	}
	else 
	{
		for (auto i = 0; i <operands.size(); ++i)
		{
			auto count = 0u; 
			for ( auto j = 0; j <operands.size(); ++j)
			{
				if (operands[i] == operands[j])
					++count; 
			}
			if (count >= operands.size()/2 + 1)
				return operands[i]; 
		}
	
		for (auto i = 0; i <operands.size(); ++i)
		{
			for ( auto j = 0; j <operands.size(); ++j)
			{
				if (operands[i] == !operands[j])
				{
					std::vector<mign_function> new_operands; 
					for (auto h = 0; h < operands.size(); ++h)
					{
						if ( h != i && h != j)
							new_operands.push_back(operands[h]); 
					}
					return create_maj_general(new_operands); 
				}
			}
		}

		std::vector<mign_function> children; 
		for (auto i = 0; i<operands.size(); ++i)
		{
			children.push_back(operands[i]); 
		}
		std::sort(children.begin(),children.begin() + children.size()); // in order 
	
		auto key = children; 
	
		if (!(_disable))
		{
			const auto it = maj_strash.find(key); 
		    if ( it != maj_strash.end() )
		    {
				 return mign_function( it->second, node_complement );
		    }
		}
	
		const auto node = add_vertex(g); 
		++_num_maj; 
	
		for ( auto i =0 ;i<operands.size(); ++i)
		{
			const auto x = add_edge(node, children[i].node, g).first; 
			complement[x] = children[i].complemented; 
		}
	
		maj_strash[key] = node; 
		return mign_function(node, node_complement); 
	}
}

mign_function mign_graph::create_maj_general_10 (const std::vector<mign_function> operands)
{
	auto node_complement = false; //
	
	if (operands.size() == 1)
		return operands[0]; 
	// special cases not valid if nin is imposted 
	if ((_nin) && (*_nin != 3))
	{
		std::vector<mign_function> children; 
		for (auto i = 0; i<operands.size(); ++i)
		{
			children.push_back(operands[i]); 
		}
		std::sort(children.begin(),children.begin() + children.size()); // in order 
	
		auto key = children; 
		if (!(_disable))
		{
			const auto it = maj_strash.find(key); 
		    if ( it != maj_strash.end() )
		    {
				 return mign_function( it->second, node_complement );
		    }
		}
		
		const auto node = add_vertex(g); 
		++_num_maj; 
	
		for ( auto i =0 ;i<operands.size(); ++i)
		{
			const auto x = add_edge(node, children[i].node, g).first; 
			complement[x] = children[i].complemented; 
		}
	
		maj_strash[key] = node; 
		return mign_function(node, node_complement); 
	}
	else 
	{

		std::vector<mign_function> children; 
		for (auto i = 0; i<operands.size(); ++i)
		{
			children.push_back(operands[i]); 
		}
		std::sort(children.begin(),children.begin() + children.size()); // in order 
	
		auto key = children; 
	
		if (!(_disable))
		{
			const auto it = maj_strash.find(key); 
		    if ( it != maj_strash.end() )
		    {
				 return mign_function( it->second, node_complement );
		    }
		}
	
		const auto node = add_vertex(g); 
		++_num_maj; 
	
		for ( auto i =0 ;i<operands.size(); ++i)
		{
			const auto x = add_edge(node, children[i].node, g).first; 
			complement[x] = children[i].complemented; 
		}
	
		maj_strash[key] = node; 
		return mign_function(node, node_complement); 
	}
}
	
mign_function mign_graph::create_maj ( std::vector<mign_function> operands)
{
	/*TODO Special cases*/
	
	const auto size = operands.size(); 
	
	if ((_nin))
	{
		num_input = *_nin; 
		
		if (size == num_input)	
		{
		return create_maj_general(operands); 
		}
		
		else if (size < num_input)
		{
			auto volte = (num_input - size )/2; 
			for ( auto i = 0; i < volte;++i )
			{
				//std::cout << " i = " << i << std::endl; 
				operands.push_back(get_constant(false));
				operands.push_back(get_constant(true));
			}
			
			return create_maj_general(operands); 
		}
		
		else //(size > num_input)
		{
			std::cout << " size = " << size << std::endl; 
			std::cout << "Cannot be reduced " << std::endl; 
			return create_maj_general(operands); 
		}
	}
	else 
		return create_maj_general(operands); 

}	

mign_function mign_graph::create_maj_10 ( std::vector<mign_function> operands)
{
	/*TODO Special cases*/
	
	//auto node_complement = false; // could be useful 
	const auto size = operands.size(); 
	
	// TODO devo aggiungere anche l'opzione non solo omogeneo, ma heterogeneo con al massimo nodi da 11 per esempio. In quel caso pero devo aggiungere un opzione in piu, tipo un comando vero e falso che controlli questa opzione 
	if ((_nin))
	{
		num_input = *_nin; 
		
		//std::cout << " numero di input = " << num_input << std::endl; 
		if (size == num_input)	
		{
			//std::cout << " size 1 = " << size << std::endl; 
		return create_maj_general_10(operands); 
		}
		
		else if (size < num_input)
		{
			//std::cout << " size 2 = " << size << std::endl; 
			auto volte = (num_input - size )/2; 
			for ( auto i = 0; i < volte;++i )
			{
				//std::cout << " i = " << i << std::endl; 
				operands.push_back(get_constant(false));
				operands.push_back(get_constant(true));
			}
			
			return create_maj_general_10(operands); 
		}
		
		else //(size > num_input)
		{
			std::cout << " size = " << size << std::endl; 
			std::cout << "Cannot be reduced " << std::endl; 
			return create_maj_general_10(operands); 
		}
	}
	else 
		return create_maj_general_10(operands); 

}	

mign_function mign_graph::create_and (std::vector<mign_function> operands)
{
	const auto num = operands.size(); 
	
	// un AND can be realized using n-1 number of 0 as input, where n is the number of input 
	for (auto i = 0; i < (num-1); ++i)
	{
		operands.push_back(get_constant(false)); 
	}

	return create_maj(operands); 

}

mign_function mign_graph::create_or (std::vector<mign_function> operands)
{
	//std::cout<< " CReate or generale" << std::endl; 
	const auto num = operands.size(); 
	
	// un AND can be realized using n-1 number of 0 as input, where n is the number of input 
	for (auto i = 0; i < (num - 1) ; ++i)
	{
		operands.push_back(get_constant(true)); 
	}
	
	return create_maj(operands); 

}

// function to creat the threshold function. polarity = 0 means >= polarity = 1 means < =. x1 + x2 + .... + xn >= threshold. This is possible with one majority of unlimited size. 

mign_function mign_graph::create_threshold (std::vector<mign_function> operands, unsigned T, unsigned polarity, std::vector<unsigned> weigths)
{
	
	auto count = 0u; 
	//const auto s = operands.size() - 1; 
	
		for (auto& w :weigths)
		{
			if ( w != 1)
			{
				for ( auto x = 2; x<=w; ++x)
				{
					mign_function doppia; 
					doppia = operands[count]; 
					operands.push_back(doppia); 
				}
			}
			++count; 
		}
		const auto num = operands.size(); 
	
		auto th = 0u; 
		if (polarity == 1)
		{
			th = num - T; 
		}
		else 
			th = T; 
		
	

	if (th > num/2)
	{
		auto narity = 2* th -1; 
		auto filling = narity - num; 
		if (polarity == 1)
		{
			for ( auto i = 0; i <num; ++ i)
			{
				operands[i].complemented = !operands[i].complemented; 
			} 
		}
		if (filling > 0)
		{
			for (auto i = 0; i < filling ; ++ i)
			{
				operands.push_back(get_constant(false));
			}
		}
		return create_maj(operands);
		
	}
	else {
		auto narity = 2*(num-th) + 1; 
		auto filling = narity - num; 
		if (polarity == 1)
		{
			for ( auto i = 0; i <num; ++ i)
			{
				operands[i].complemented = !operands[i].complemented; 
			} 
		}
		if (filling > 0)
		{
			for (auto i = 0; i < filling ; ++ i)
			{
				operands.push_back(get_constant(true));
			}
		}
		return create_maj(operands);
	}
		/*if (polarity == 0) // case x1 + x2 + x3 ... > t 
		{
			 
			const signed a = 2*th - 1 - num; 
			const signed b = num + 1 - 2*th; 
			//std::cout << " a = " << a << std::endl;
			//std::cout << " b = " << b << std::endl;
			if ( a > 0){
				
			for (auto i = 0; i < a; ++i)
			{
				operands.push_back(get_constant(false)); 
			}
		 	}
			if (b > 0){
				//std::cout << " ERROR 3" << std::endl; 
				
			for (auto i = 0; i < b; ++i)
			{
				operands.push_back(get_constant(true)); 
			}
			}
			
			return create_maj(operands); 
		}
		else //if (polarity == 1)
		{
			const signed c = num - 1 - 2*th; 
			const signed d = 2*th - 1 - num;
			
			for ( auto i = 0; i <num; ++ i)
			{
				operands[i].complemented = !operands[i].complemented; 
			} 
			
			if (c > 0)
			{
			for (auto i = 0; i < c; ++ i)
			{
				operands.push_back(get_constant(false)); 
			}
			}
			if ( d > 0)
			{
			for (auto i = 0; i < d; ++ i)
			{
				operands.push_back(get_constant(true)); 
			}
			}
			return create_maj(operands);
		}*/
	}
	
unsigned mign_graph::fanin_count( node_t n ) const
{
  return boost::out_degree( n, g );
}

unsigned mign_graph::fanout_count( node_t n ) const
{
	
  return (*fanout)[n];
}

const std::vector<mign_graph::node_t>& mign_graph::parents( node_t n ) const
{
  return (*parentss)[n];
}

unsigned mign_graph::level( node_t n ) const
{
  return (*levels)[n];
}


bool mign_graph::is_input( node_t n ) const
{
  return fanin_count( n ) == 0u;
}

bool mign_graph::is_output(node_t n) const
{
	bool flag = false; 
	for (auto& output : outputs())
	{
		if (n == output.first.node)
		flag = true; 	
	}
	return flag; 
}
//const std::string& mign_graph::name() const
//{
//  return _name;
//}

boost::optional<unsigned> mign_graph::nin() const
{
	if (_nin)
	{
		return _nin; 
	}
	else 
	{std::cout << "Heterogeneous" << std::endl;  
	return _nin; 
	}
}

mign_bitmarks& mign_graph::bitmarks()
{
  return *_bitmarks;
}

const mign_bitmarks& mign_graph::bitmarks() const
{
  return *_bitmarks;
}

void mign_graph::init_refs()
{
  compute_fanout();
  ref_count = *fanout;
}

unsigned mign_graph::get_ref( mign_node n ) const
{
  return ref_count[n];
}

unsigned mign_graph::inc_ref( mign_node n )
{
  return ref_count[n]++;
}

unsigned mign_graph::dec_ref( mign_node n )
{
  assert( ref_count[n] > 0 );
  return --ref_count[n];
}

void mign_graph::inc_output_refs()
{
  for ( const auto& output : outputs() )
  {
    ++ref_count[output.first.node];
  }
}

std::size_t mign_graph::size() const
{
  return num_vertices( g );
}

unsigned mign_graph::num_gates() const
{
  return _num_maj;
}

const mign_graph::graph_t& mign_graph::graph() const
{
  return g;
}



mign_graph::graph_t& mign_graph::graph()
{
  return g;
}

const mign_graph::input_vec_t& mign_graph::inputs() const
{
  return _inputs;
}

const mign_graph::output_vec_t& mign_graph::outputs() const
{
  return _outputs;
}

const std::string& mign_graph::input_name( mign_node n ) const
{
  return _inputs[_input_to_id.at( n )].second;
}

std::string mign_graph::output_name( mign_node n ) const
{
	for ( const auto& x : _outputs)
	{
		if ( x.first.node == n )
    {
			return x.second;
    }
	}

  return std::string();
}

const std::string& mign_graph::name() const
{
  return _name;
}

void mign_graph::set_name( const std::string& name )
{
  _name = name;
}

const unsigned mign_graph::input_index( mign_node n ) const
{
  return _input_to_id.at( n );
}


mign_graph::vertex_range_t mign_graph::vertices() const
{
  return boost::make_iterator_range( boost::vertices( g ) );
}

std::vector<mign_function> mign_graph::children( mign_node n ) const
{
  std::vector<mign_function> c;
  for ( const auto& e : boost::make_iterator_range( boost::out_edges( n, g ) ) )
  {
    c.push_back( mign_function( boost::target( e, g ), complement[e] ) );
  }
  return c;
}


std::vector<mign_graph::node_t> mign_graph::topological_nodes() const
{
  std::vector<node_t> top( num_vertices( g ) );
  boost::topological_sort( g, top.begin() );
  return top;
}


bool mign_graph::has_cover() const
{
  return (bool)_cover;
}

bool mign_graph::has_multi_cover() const
{
	if (_multi_cover.size() > 0)
		return true; 
	else return false; 
}


const mign_cover& mign_graph::cover() const
{
  return *_cover;
}

std::vector<mign_cover> mign_graph::multi_cover() const
{
	std::vector<mign_cover> provv; 
	for ( auto y = 0; y < _multi_cover.size(); ++y)
	{
		provv.push_back(*_multi_cover[y]); 
	}
  return provv; 
}

void mign_graph::set_cover( const mign_cover& other )
{
  if ( !_cover )
  {
    _cover = std::make_shared<mign_cover>( 0u, *this );
  }
  *_cover = other;
}



void mign_graph::set_multi_cover( const std::vector<mign_cover>& other )
{
	//unsigned y = 0; 
	//std::cout << " set multi cover " << std::endl; 
	_multi_cover.resize(other.size()); 
	//std::cout << " set multi cover dopo resize " << other.size() << std::endl; 
	for ( auto i = 0; i < other.size(); ++i)
	{
	    if ( !_multi_cover[i] )
	    {
	      _multi_cover[i] = std::make_shared<mign_cover>( 0u, *this );
	    }
		//std::cout << i << std::endl; 
		*_multi_cover[i] = other[i]; 
	}

}

unsigned mign_graph::size_multi_cover()
{
	return _multi_cover.size();
}

mign_cover& mign_graph::take_one_cover(unsigned i) 
{
	return *_multi_cover[i]; 
}

/******************************************************************************
 * mign_fuction                                                            *
 *********************************************************************/
mign_function::mign_function( mign_node node, bool complemented )
  : node( node ),
    complemented( complemented )
{
}

bool mign_function::operator==( const mign_function& other ) const
{
  return node == other.node && complemented == other.complemented;
}

bool mign_function::operator!=( const mign_function& other ) const
{
  return !operator==( other );
}

bool mign_function::operator<( const mign_function& other ) const
{
  if ( node < other.node )
  {
    return true;
  }
  else if ( node == other.node )
  {
    return !complemented && other.complemented;
  }
  else
  {
    return false;
  }
}

bool mign_function::operator>( const mign_function& other ) const
{
  return ( other < *this );
}

mign_function mign_function::operator!() const
{
  return mign_function( node, !complemented );
}

mign_function mign_function::operator^( bool value ) const
{
  return mign_function( node, complemented != value );
}


}

// Local Variables:
// c-basic-offset: 2
// eval: (c-set-offset 'substatement-open 0)
// eval: (c-set-offset 'innamespace 0)
// End:
