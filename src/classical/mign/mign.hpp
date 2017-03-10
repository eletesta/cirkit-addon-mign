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

/**
 * @file mign.hpp
 *
 * @brief MIG-n graph
 *
 * @author Eleonora Testa
 * @since  2.3
 */

#ifndef MIGN_HPP
#define MIGN_HPP

#include <functional>
#include <string>
#include <vector>
#include <unordered_map>

#include <boost/graph/adjacency_list.hpp>
#include <boost/range/iterator_range.hpp>
#include <boost/optional.hpp>
#include <boost/functional/hash.hpp> 

#include <core/utils/graph_utils.hpp>
#include <core/utils/hash_utils.hpp>
#include <core/utils/dirty.hpp>
//#include "hash_utils.hpp"

namespace cirkit
{
	// if you want to stroe more info about verteces
//<boost::property<boost::vertex_name_t, my_vertex_info> = to store the information about the type of nodes with want in term of number of inputs. 
//struct my_vertex_info
//{

	// my_vertex_info() {}
	// my_vertex_info( const std::string& name ) : name( name ) {}

	// std::string name;
	// unsigned        ninput = 3; // default value of the number of input 
	//};

	// hashing part 

	
using mign_graph_t = digraph_t<boost::no_property, boost::property<boost::edge_complement_t,bool>>; 
using mign_node    = vertex_t<mign_graph_t>; 

class mign_function 
{
public: 
	mign_function (mign_node node = 0, bool complemented = false); 
	
	bool operator==(const mign_function& other) const; 
	bool operator!=(const mign_function& other) const; 
	bool operator<(const mign_function& other) const; 
	bool operator>(const mign_function& other) const; 
	
	mign_function operator!() const; 
	mign_function operator^(bool value) const; 
	
public: 
	mign_node node; 
	bool complemented; 
};
}

namespace std
{
	template<>
	struct hash<cirkit::mign_function>
	{
		inline std::size_t operator() (const cirkit::mign_function& f) const
		{
			return (f.node << 1u) + static_cast<int> (f.complemented); 
		}
	};
}

namespace cirkit
{

class mign_cover; 
	
class mign_graph 
{
public: 
	using graph_t = mign_graph_t; 
	using node_t = vertex_t<graph_t>; 
	
	using input_vec_t = std::vector<std::pair<node_t, std::string>>; 
	using output_vec_t = std::vector<std::pair<mign_function, std::string>>; 
	using vertex_range_t = boost::iterator_range<boost::graph_traits<mign_graph_t>::vertex_iterator>;
	
public: 
	mign_graph(/*const std::string& name = std::string(), */boost::optional<unsigned> nin = boost::optional<unsigned>()); 
	
    void compute_fanout();
    void compute_parents();
    void compute_levels();
	
	//void disable_strashing(); 
	
	unsigned fanin_count (node_t n) const; 
	unsigned fanout_count (node_t n) const; 
	
	mign_function get_constant(bool value) const; 
	mign_function create_pi (const std::string& name); 
	void create_po (const mign_function& f, const std::string name); 
	mign_function create_maj_general (const std::vector<mign_function>);
	mign_function create_maj(const std::vector<mign_function>);
	mign_function create_maj_general_10 (const std::vector<mign_function>);
	mign_function create_maj_10(const std::vector<mign_function>); 
	mign_function create_maj_fo_restricted(const std::vector<mign_function>, const unsigned int); 
	mign_function create_and(const std::vector<mign_function>); 
	mign_function create_or(const std::vector<mign_function>);
	mign_function create_threshold (std::vector<mign_function> operands, unsigned th, unsigned polarity, std::vector<unsigned> weigths); 

    const std::vector<node_t>& parents( node_t n ) const;
    unsigned level( node_t n ) const;
	

    bool is_input( node_t n ) const;
    bool is_maj( node_t n ) const;
	bool is_output (node_t n) const; 
	//bool is_and (node_t n) const; 
	//bool is_or (node_t n) const; 
	
    const std::string& name() const;
    void set_name( const std::string& name );
    std::size_t size() const;
    unsigned num_gates() const;
    //unsigned num_maj() const;
	//unsigned num_input() const; 
	
	boost::optional<unsigned> nin() const; 
    const graph_t& graph() const;
    graph_t& graph();
    const input_vec_t& inputs() const;
    const output_vec_t& outputs() const;
    const std::string& input_name( mign_node n ) const;
	std::string output_name( mign_node n ) const;
    const unsigned input_index( mign_node n ) const;
    std::vector<mign_function> children( mign_node n ) const;
    vertex_range_t vertices() const;
    std::vector<node_t> topological_nodes() const;
	
    bool has_cover() const;
    const mign_cover& cover() const;
    bool has_multi_cover() const; 
    std::vector<mign_cover> multi_cover() const;
	unsigned size_multi_cover(); 
	mign_cover& take_one_cover(unsigned i); 
    void set_cover( const mign_cover& other );
	void set_multi_cover (const std::vector<mign_cover>& other); 
	
    inline void structural_hashing( bool disable ) { _disable = disable; }
    inline bool has_structural_hashing() const         { return !(_disable); }
	
private: 
    graph_t g;
    node_t  constant;

	std::string  _name;
    input_vec_t  _inputs;
    output_vec_t _outputs;
	boost::optional<unsigned> _nin; 
	boost::optional<bool> _enable_strashing; 
	
    std::unordered_map<mign_node, unsigned> _input_to_id; // difference with map is in the way they are accessed 

    boost::property_map<graph_t, boost::edge_complement_t>::type         complement;
	
    dirty<std::vector<unsigned>>                                         fanout;
	//dirty<std::vector<unsigned>                                                fanout;
    dirty<std::vector<std::vector<node_t>>>                              parentss;
    dirty<std::vector<unsigned>>                                         levels;

	
	  unsigned                                                           _num_maj = 0u;
	  unsigned 															 num_input = 0u; 
    //boost::property_map<graph_t, boost::vertex_ninput_t>::type         ninput;
	using maj_strash_key_t = std::vector<mign_function>;   
	std::unordered_map<maj_strash_key_t, node_t, hash<maj_strash_key_t> > maj_strash; // store and access it faster 
	
    std::unordered_map<maj_strash_key_t, std::pair<node_t,unsigned>, hash<maj_strash_key_t> > maj_strash_fo_restricted;

    std::shared_ptr<mign_cover>                                           _cover = nullptr;
	std::vector<std::shared_ptr<mign_cover>>                              _multi_cover ;
	
	bool                                                                  _disable = false; 
    
	
}; 
}

#endif

// Local Variables:
// c-basic-offset: 2
// eval: (c-set-offset 'substatement-open 0)
// eval: (c-set-offset 'innamespace 0)
// End:
