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
 * @file mign_simulate.hpp
 *
 * @brief MIGN simulation
 *
 * @author Eleonora Testa 
 * @since  2.3
 */

#ifndef MIGN_SIMULATE_HPP
#define MIGN_SIMULATE_HPP

#include <map>
#include <unordered_map>

#include <core/properties.hpp>
#include <core/utils/timer.hpp>
#include <classical/utils/truth_table_utils.hpp>
#include <classical/mign/mign.hpp>
#include <classical/mign/mign_utils.hpp>
#include <classical/mign/mign_dfs.hpp>


namespace cirkit
{

/******************************************************************************
 * Abstract class for simulators                                              *
 ******************************************************************************/

template<typename T>
class mign_simulator
{
public:
  /**
   * @brief Simulator routine when input is encountered
   *
   * @param node MIGN node reference in the `mig' graph
   * @param name Name of that node
   * @param pos  Position of that node (usually wrt. to the input vector)
   * @param mign  MIGN
   */
  virtual T get_input( const mign_node& node, const mign_graph& mign ) const = 0;
  virtual T get_constant() const = 0;
  virtual T invert( const T& v ) const = 0;
  virtual T maj_op( const mign_node& node, const std::vector<T>& v) const = 0;

  virtual bool terminate( const mign_node& node, const mign_graph& mign ) const
  {
    return false;
  }
};

/******************************************************************************
 * Several simulator implementations                                          *
 ******************************************************************************/

class mign_tt_simulator : public mign_simulator<tt>
{
public:
  tt get_input( const mign_node& node, const mign_graph& mign ) const;
  tt get_constant() const;
  tt invert( const tt& v ) const;
  tt maj_op( const mign_node& node,const std::vector<tt>& v) const;
};

template<typename T>
class mign_partial_node_assignment_simulator : public mign_simulator<T>
{
public:
  mign_partial_node_assignment_simulator( const mign_simulator<T>& total_simulator,
                                         const std::map<mign_node, T>& assignment,
                                         const T& default_value )
    : total_simulator( total_simulator ),
      assignment( assignment ),
      default_value( default_value ) {}

  T get_input( const mign_node& node, const mign_graph& mign ) const
  {
    auto it = assignment.find( node );
    if ( it == assignment.end() )
    {
      std::cout << "[w] no assignment given for '" << node << "', assume default" << std::endl;
      return default_value;
    }
    else
    {
      return it->second;
    }
  }

  T get_constant() const { return total_simulator.get_constant(); }
  T invert( const T& v ) const { return total_simulator.invert( v ); }
  
  T maj_op( const mign_node& node, const std::vector<T>& v) const
  {
    auto it = assignment.find( node );
    if ( it == assignment.end() )
    {
      return total_simulator.maj_op( node, v);
    }
    else
    {
      return it->second;
    }
  }

  bool terminate( const mign_node& node, const mign_graph& mign ) const
  {
    return assignment.find( node ) != assignment.end();
  }

private:
  const mign_simulator<T>& total_simulator;
  const std::map<mign_node, T>& assignment;
  T default_value;
};

/******************************************************************************
 * DFS visitor for actual simulation                                          *
 ******************************************************************************/

using mign_node_color_map = std::map<mign_node, boost::default_color_type>;

template<typename T>
struct simulate_mign_node_visitor : public mign_dfs_visitor
{
public:
  simulate_mign_node_visitor( const mign_graph& mign, const mign_simulator<T>& simulator, std::map<mign_node, T>& node_values )
    : mign_dfs_visitor( mign ),
      simulator( simulator ),
      node_values( node_values ) {}

  void finish_input( const mign_node& node, const mign_graph& mign )
  {
    node_values[node] = simulator.get_input( node, mign );
  }

  void finish_constant( const mign_node& node, const mign_graph& mign )
  {
    node_values[node] = simulator.get_constant();
  }

  void finish_maj_node( const mign_node& node, const std::vector<mign_function>& ciao, const mign_graph& mign )
  {
	  std::vector<T> tc(ciao.size()); 
	  //std::cout << ciao.size() << std::endl; 
	  
	  for (auto x = 0; x < ciao.size(); ++x)
	  {
		 
		  tc[x] = node_values[ciao[x].node]; 
		 
		  if ( ciao[x].complemented)
			 {
				   
			 	tc[x] = simulator.invert(tc[x]);
				
			 }  
	  }
    
    node_values[node] = simulator.maj_op( node,tc );
  }

private:
  const mign_simulator<T>& simulator;
  std::map<mign_node, T>& node_values;
};

/******************************************************************************
 * Methods to trigger simulation                                              *
 ******************************************************************************/

template<typename T>
T simulate_mign_node( const mign_graph& mign, const mign_node& node,
                     const mign_simulator<T>& simulator,
                     mign_node_color_map& colors,
                     std::map<mign_node, T>& node_values )
{
  boost::depth_first_visit( mign.graph(), node,
                            simulate_mign_node_visitor<T>( mign, simulator, node_values ),
                            boost::make_assoc_property_map( colors ),
                            [&simulator, &mign]( const mign_node& node, const mign_graph::graph_t& g ) { return simulator.terminate( node, mign ); } );
					
  return node_values[node];
}

template<typename T>
T simulate_mign_node( const mign_graph& mign, const mign_node& node,
                     const mign_simulator<T>& simulator )
{
  mign_node_color_map colors;
  std::map<mign_node, T> node_values;


  return simulate_mign_node<T>( mign, node, simulator, colors, node_values );
}

template<typename T>
T simulate_mign_function( const mign_graph& mign, const mign_function& f,
                         const mign_simulator<T>& simulator,
                         mign_node_color_map& colors,
                         std::map<mign_node, T>& node_values )
{
  T value = simulate_mign_node<T>( mign, f.node, simulator, colors, node_values );
  return f.complemented ? simulator.invert( value ) : value;
}

template<typename T>
T simulate_mign_function( const mign_graph& mign, const mign_function& f,
                         const mign_simulator<T>& simulator )
{
	
  T value = simulate_mign_node<T>( mign, f.node, simulator );
  return f.complemented ? simulator.invert( value ) : value;
}

template<typename T>
std::map<mign_function, T> simulate_mign( const mign_graph& mign, const mign_simulator<T>& simulator,
                                        const properties::ptr& settings = properties::ptr(),
                                        const properties::ptr& statistics = properties::ptr() )
{
  /* settings */
  const auto verbose = get( settings, "verbose", false );

  /* timer */
  properties_timer t( statistics );

  mign_node_color_map colors;
  std::map<mign_node, T> node_values;

  std::map<mign_function, T> results;

  for ( const auto& o : mign.outputs() )
  {
    if ( verbose )
    {
      std::cout << "[i] simulate '" << o.second << "'" << std::endl;
    }
    T value = simulate_mign_node<T>( mign, o.first.node, simulator, colors, node_values );

    /* value may need to be inverted */
    results[o.first] = o.first.complemented ? simulator.invert( value ) : value;
  }

  set( statistics, "node_values", node_values );

  return results;
}

}

#endif

// Local Variables:
// c-basic-offset: 2
// eval: (c-set-offset 'substatement-open 0)
// eval: (c-set-offset 'innamespace 0)
// End:
