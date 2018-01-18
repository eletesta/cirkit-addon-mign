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

#include "mign_cuts_paged.hpp"

#include <map>
#include <mutex>
#include <stack>
#include <type_traits>

#include <core/utils/bitset_utils.hpp>
#include <core/utils/range_utils.hpp>
#include <core/utils/timer.hpp>
#include <core/utils/program_options.hpp>
#include <classical/utils/truth_table_utils.hpp>
#include <classical/mig/mig.hpp>
#include <classical/mig/mig_utils.hpp>
#include <classical/mig/mig_from_string.hpp>

#include <classical/mign/mig_to_mign.hpp>
#include <classical/mign/mign_from_string.hpp>
#include <classical/mign/mign_simulate.hpp>
#include <classical/mign/mign_utils.hpp>
#include <classical/mign/mign_utils_majtt.hpp>

#include <formal/synthesis/exact_mig.hpp>

#include <boost/assign/std/vector.hpp>
#include <boost/optional.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <boost/graph/topological_sort.hpp>
#include <boost/range/adaptors.hpp>
#include <boost/range/algorithm.hpp>
#include <boost/range/combine.hpp>
#include <boost/range/numeric.hpp>

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

std::vector<std::pair<unsigned, unsigned>> compute_level_ranges( mign_graph& mign, unsigned& max_level )
{
  mign.compute_levels();
  mign.compute_parents();

  max_level = 0u;
  for ( const auto& o : mign.outputs() )
  {
    max_level = std::max( max_level, mign.level( o.first.node ) );
  }

  std::vector<unsigned> levels( mign.size() );
  for ( const auto& v : mign.vertices() )
  {
    levels[v] = mign.level( v );
  }

  std::vector<std::pair<unsigned, unsigned>> level_ranges( mign.size() );

  const auto top = mign.topological_nodes();
  for ( const auto& v : top )
  {
    const auto& parents = mign.parents( v );

    if ( parents.empty() )
    {
      level_ranges[v] = {levels[v], max_level};
      levels[v] = max_level;
      continue;
    }

    const auto min_parent = *boost::min_element( parents, [&levels]( const mign_node& p1, const mign_node& p2 ) {
        return levels[p1] < levels[p2];
      } );
    auto to_l = levels[min_parent] - 1u;
    level_ranges[v] = {levels[v], to_l};
    levels[v] = to_l;
  }

  return level_ranges;
}

mign_cuts_paged::mign_cuts_paged( mign_graph& mign, unsigned k, const properties::ptr& settings )
  : _mign( mign ),
    _k( k ),
    _priority( get( settings, "priority_cut" , 0u) ),
	_almost( get( settings, "almost" , 0u) ),
    _extra( get( settings, "extra", 0u ) ),
    data( _mign.size(), 2u + _extra ),
    cones( _mign.size() ), 
	almost_info( _mign.size() )
{
  unsigned max_level;
  _levels = compute_level_ranges( mign, max_level );
  enumerate();
}

const mign_graph& mign_cuts_paged::mign() const
{
  return _mign;
}

unsigned mign_cuts_paged::total_cut_count() const
{
  return data.sets_count();
}

double mign_cuts_paged::enumeration_time() const
{
  return _enumeration_time;
}

unsigned mign_cuts_paged::memory() const
{
  return data.memory() + cones.memory();
}

unsigned mign_cuts_paged::count( mign_node node ) const
{
  return data.count( node );
}

boost::iterator_range<paged_memory::iterator> mign_cuts_paged::cuts( mign_node node )
{
  return data.sets( node );
}

boost::iterator_range<paged_memory::iterator> mign_cuts_paged::cut_cones( mign_node node )
{
  return cones.sets( node );
}

tt mign_cuts_paged::simulate( mign_node node, const mign_cuts_paged::cut& c ) const
{
  std::vector<mign_node> leafs;
  for ( auto child : c )
  {
    leafs.push_back( child );
  }
  return mign_simulate_cut( _mign, node, leafs );
}

unsigned mign_cuts_paged::depth( mign_node node, const mign_cuts_paged::cut& c ) const
{
  return c.extra( 0u );
}

unsigned mign_cuts_paged::size( mign_node node, const mign_cuts_paged::cut& c ) const
{
  return c.extra( 1u );
}

tt mign_cuts_paged::tt_func( mign_node node, const int cut) const
{
	auto prova = almost_info[node][cut]; 
    return std::get<0>(prova);
}

int mign_cuts_paged::is_maj_or_almost( mign_node node, const int cut) const  // = 0 é una majority, = 1 é un almost majority , =2 é una lmost majority < 
{
	auto prova = almost_info[node][cut]; 
    return std::get<1>(prova);
}

tt mign_cuts_paged::tt_reminder( mign_node node, const int cut) const
{
  return std::get<2>(almost_info[node][cut]);
}

tt mign_cuts_paged::tt_comp_in( mign_node node, const int cut) const
{
  return std::get<3>(almost_info[node][cut]);
}

std::string mign_cuts_paged::tt_exact( mign_node node, const int cut) const
{
  return std::get<4>(almost_info[node][cut]);
}


unsigned mign_cuts_paged::index( const mign_cuts_paged::cut& c ) const
{
  return data.index( c );
}

mign_cuts_paged::cut mign_cuts_paged::from_address( unsigned address )
{
  return data.from_address( address );
}

void mign_cuts_paged::foreach_cut( const std::function<void(mign_node, cut&)>& func )
{
  for ( const auto& n : _mign.vertices() )
  {
    for ( auto cut : cuts( n ) )
    {
      func( n, cut );
    }
  }
}

void mign_cuts_paged::enumerate()
{
  reference_timer t( &_enumeration_time );

  /* topsort */
  const auto top = _mign.topological_nodes();

  /* loop */
  _top_index = 0u;
  for ( auto n : top )
  {
    if ( _mign.is_input( n ) )
    {
      /* constant */
      if ( n == 0u )
      {
        data.assign_empty( 0u, get_extra( 0u, 0u ) );
        cones.assign_empty( 0u );
      }
      /* PI */
      else
      {
        data.assign_singleton( n, n, get_extra( 0u, 1u ) );
        cones.assign_singleton( n, n );
      }
    }
    else
    {
      data.append_begin( n );
      cones.append_begin( n );

      std::vector<mign_node> cns;
      for ( const auto& c : _mign.children( n ) )
      {
        cns.push_back( c.node );
      }
	  
      enumerate_node_with_bitsets( n, cns );

      data.append_singleton( n, n, get_extra( 0u, 1u ) );
      cones.append_singleton( n, n );
    }

    _top_index++;
  }
}

void mign_cuts_paged::merge_cut( local_cut_vec_t& local_cuts, const boost::dynamic_bitset<>& new_cut, unsigned min_level, const boost::dynamic_bitset<>& new_cone , mign_node parent) const
{
  /* too large? */
  if ( new_cut.count() > _k ) { return; }
  
  auto first_subsume = true;
  auto add = true;
  
  auto l = 0u;
  while ( l < local_cuts.size() )
  {
    auto cut = std::get<0>( local_cuts[l] );
		
    /* same cut */
    if ( cut == new_cut ) { add = false; break; }

    /* cut subsumes new_cut */  
    else if ( ( cut & new_cut ) == cut ) {  add = false; break; }

    /* new_cut subsumes cut */
    else if ( ( cut & new_cut ) == new_cut )  
    {
		
	  add = false;
      if ( first_subsume )
      {
        local_cuts[l] = std::make_tuple( new_cut, min_level, new_cone);
        first_subsume = false;
      }
      else
      {
        local_cuts[l] = local_cuts.back();
        local_cuts.pop_back();
      }
    }

    ++l;
  }

  if ( add )
  {
    local_cuts.push_back( std::make_tuple( new_cut, min_level, new_cone));
  }
  
}

mign_cuts_paged::local_cut_vec_t mign_cuts_paged::enumerate_local_cuts( mign_node n1, mign_node n2, unsigned max_cut_size, mign_node parent )
{
  local_cut_vec_t local_cuts;

  for ( const auto& c1 : boost::combine( cuts( n1 ), cut_cones( n1 ) ) )
  {
    for ( const auto& c2 : boost::combine( cuts( n2 ), cut_cones( n2 ) ) )
    {
      
        auto min_level = std::numeric_limits<unsigned>::max();
        boost::dynamic_bitset<> new_cut( max_cut_size );
        auto f = [&new_cut, &min_level, this]( unsigned pos ) {
          new_cut.set( pos );
          min_level = std::min( min_level, this->_levels[pos].second );
        };
        std::for_each( boost::get<0>( c1 ).begin(), boost::get<0>( c1 ).end(), f );
        std::for_each( boost::get<0>( c2 ).begin(), boost::get<0>( c2 ).end(), f );
       

        boost::dynamic_bitset<> new_cone( max_cut_size );
        auto f2 = [&new_cone]( unsigned pos ) {
          new_cone.set( pos );
        };
        std::for_each( boost::get<1>( c1 ).begin(), boost::get<1>( c1 ).end(), f2 );
        std::for_each( boost::get<1>( c2 ).begin(), boost::get<1>( c2 ).end(), f2 );
       
        merge_cut( local_cuts, new_cut, min_level, new_cone, parent );
      
    }
  }

  return local_cuts;
}
mign_cuts_paged::local_cut_vec_t mign_cuts_paged::enumerate_local_cuts( mign_node n1, mign_node n2, mign_node n3, unsigned max_cut_size , mign_node parent)
{
  local_cut_vec_t local_cuts;
  
  for ( const auto& c1 : boost::combine( cuts( n1 ), cut_cones( n1 ) ) )
  {
    for ( const auto& c2 : boost::combine( cuts( n2 ), cut_cones( n2 ) ) )
    {
      for ( const auto& c3 : boost::combine( cuts( n3 ), cut_cones( n3 ) ) )
      {
        auto min_level = std::numeric_limits<unsigned>::max();
        boost::dynamic_bitset<> new_cut( max_cut_size );
        auto f = [&new_cut, &min_level, this]( unsigned pos ) {
			
          new_cut.set( pos );
          min_level = std::min( min_level, this->_levels[pos].second );
        };
        std::for_each( boost::get<0>( c1 ).begin(), boost::get<0>( c1 ).end(), f );
        std::for_each( boost::get<0>( c2 ).begin(), boost::get<0>( c2 ).end(), f );
        std::for_each( boost::get<0>( c3 ).begin(), boost::get<0>( c3 ).end(), f );

        boost::dynamic_bitset<> new_cone( max_cut_size );
        auto f2 = [&new_cone]( unsigned pos ) {
          new_cone.set( pos );
        };
        std::for_each( boost::get<1>( c1 ).begin(), boost::get<1>( c1 ).end(), f2 );
        std::for_each( boost::get<1>( c2 ).begin(), boost::get<1>( c2 ).end(), f2 );
        std::for_each( boost::get<1>( c3 ).begin(), boost::get<1>( c3 ).end(), f2 );

        merge_cut( local_cuts, new_cut, min_level, new_cone , parent);
      }
    }
  }
   
  return local_cuts;
}

mign_cuts_paged::local_cut_vec_t mign_cuts_paged::enumerate_local_cuts( mign_node n1, mign_node n2, mign_node n3, mign_node n4, mign_node n5, unsigned max_cut_size , mign_node parent)
{
  local_cut_vec_t local_cuts;
 
  for ( const auto& c1 : boost::combine( cuts( n1 ), cut_cones( n1 ) ) )
  {
    for ( const auto& c2 : boost::combine( cuts( n2 ), cut_cones( n2 ) ) )
    {
      for ( const auto& c3 : boost::combine( cuts( n3 ), cut_cones( n3 ) ) )
      {
		  for ( const auto& c4 : boost::combine( cuts( n4 ), cut_cones( n4 ) ))
		  {
			  for (  const auto& c5 : boost::combine( cuts( n5 ), cut_cones( n5 ) ))
			  {
		          auto min_level = std::numeric_limits<unsigned>::max();
		          boost::dynamic_bitset<> new_cut( max_cut_size );
		          auto f = [&new_cut, &min_level, this]( unsigned pos ) {
		            new_cut.set( pos );
				
		            min_level = std::min( min_level, this->_levels[pos].second );
		          };
		          std::for_each( boost::get<0>( c1 ).begin(), boost::get<0>( c1 ).end(), f );
		          std::for_each( boost::get<0>( c2 ).begin(), boost::get<0>( c2 ).end(), f );
		          std::for_each( boost::get<0>( c3 ).begin(), boost::get<0>( c3 ).end(), f );
				  std::for_each( boost::get<0>( c4 ).begin(), boost::get<0>( c4 ).end(), f );
				  std::for_each( boost::get<0>( c5 ).begin(), boost::get<0>( c5 ).end(), f );

		          boost::dynamic_bitset<> new_cone( max_cut_size );
		          auto f2 = [&new_cone]( unsigned pos ) {
		            new_cone.set( pos );
		          };
		          std::for_each( boost::get<1>( c1 ).begin(), boost::get<1>( c1 ).end(), f2 );
		          std::for_each( boost::get<1>( c2 ).begin(), boost::get<1>( c2 ).end(), f2 );
		          std::for_each( boost::get<1>( c3 ).begin(), boost::get<1>( c3 ).end(), f2 );
				  std::for_each( boost::get<1>( c4 ).begin(), boost::get<1>( c4 ).end(), f2 );
				  std::for_each( boost::get<1>( c5 ).begin(), boost::get<1>( c5 ).end(), f2 );

		          merge_cut( local_cuts, new_cut, min_level, new_cone , parent);
			  }
		  }  
      }
    }
  }

  return local_cuts;
}

mign_cuts_paged::local_cut_vec_t mign_cuts_paged::enumerate_local_cuts( const std::vector<mign_node>& ns, unsigned max_cut_size, mign_node n )
{
  local_cut_vec_t local_cuts, double_cuts;
  
  if (ns.size() == 2u)
  {
  	local_cuts = enumerate_local_cuts( ns[0u], ns[1u], max_cut_size , n);
  }
  if ( ns.size() == 3u )
  {
    local_cuts = enumerate_local_cuts( ns[0u], ns[1u], ns[2u], max_cut_size, n );
  }
  else 
  {
	  const auto num_children = ns.size(); 
	  
	  std::vector<unsigned> a(num_children+1,0u); 
	  std::vector<unsigned> m(num_children+1,2u); 
	  std::vector<unsigned> prova(num_children+1,0u);
	   
	  auto tot = 0; 
	  for (auto y = 1; y < num_children+1; ++y)
	  {
		  m[y] = count(ns[tot]); 
		  tot++;  
	  }  
	   std::vector<paged_memory::iterator> as;
	   std::vector<paged_memory::iterator> as_c;
	 
	  for ( auto i = 0u; i < num_children; ++i )
	  {
			as.push_back(cuts( ns[i] ).begin());
			as_c.push_back(cut_cones(ns[i]).begin()); 
	  }

	  while ( true )
	  {
		   
          auto min_level = std::numeric_limits<unsigned>::max();
          boost::dynamic_bitset<> new_cut( max_cut_size );
          auto f = [&new_cut, &min_level, this]( unsigned pos ) {
          new_cut.set( pos );
          min_level = std::min( min_level, this->_levels[pos].second );
          };
		  
		  for (auto x = 0; x < num_children; ++x)
		  {
		  	std::for_each( (*as[x]).begin(), (*as[x]).end(), f );
		  }
		  
          boost::dynamic_bitset<> new_cone( max_cut_size );
          auto f2 = [&new_cone]( unsigned pos ) {
            new_cone.set( pos );
          };
		  for (auto x = 0; x < num_children; ++x)
		  {
		  	std::for_each( (*as_c[x]).begin(), (*as_c[x]).end(), f2 );
		  }
          merge_cut( local_cuts, new_cut, min_level, new_cone , n);
		 
		  auto flag = 0; 
		  for ( auto i = 1u; i <= num_children; ++i )
		  {
		    if (a[i] == m[i] -1) 
		    	++flag; 
		  }
		  if (flag == num_children)
		  break; 

	      auto j = a.size()-1;    
	      auto c = as.size() - 1;  
	      while ( a[j] == m[j]-1 )
	      {
		    auto x = j--;
		    auto y = c--; 
	        as[y] = cuts( ns[y] ).begin();
		    a[x] = prova[0]; 
		    as_c[y] = cut_cones(ns[y]).begin();   
	      }
		
	      if (!j) {break;}
	      as[c]++; 
	      a[j]++;  
	      as_c[c]++; 
	   }
   }
   
   boost::sort( local_cuts, []( const std::tuple<boost::dynamic_bitset<>, unsigned, boost::dynamic_bitset<>>& e1,
                                const std::tuple<boost::dynamic_bitset<>, unsigned, boost::dynamic_bitset<>>& e2 ) {
                  return ( std::get<1>( e1 ) > std::get<1>( e2 ) ) || ( std::get<1>( e1 ) == std::get<1>( e2 ) && std::get<0>( e1 ).count() < std::get<0>( e2 ).count() ); }
 				  );
				  
   std::tuple<boost::dynamic_bitset<>, unsigned, boost::dynamic_bitset<>> first_local; 
   double_cuts.push_back(local_cuts[0]); 
   local_cuts.erase(local_cuts.begin() + 0); 
   
   boost::sort( local_cuts, []( const std::tuple<boost::dynamic_bitset<>, unsigned, boost::dynamic_bitset<>>& e1,
                               const std::tuple<boost::dynamic_bitset<>, unsigned, boost::dynamic_bitset<>>& e2 ) {
                 return ( std::get<1>( e1 ) < std::get<1>( e2 ) ) || ( std::get<1>( e1 ) == std::get<1>( e2 ) && std::get<0>( e1 ).count() < std::get<0>( e2 ).count() ); }
				  );
  for (auto x = 0; x < local_cuts.size(); ++x)
  {
	  double_cuts.push_back(local_cuts[x]); 
  }	
  if (_priority > 0)  
  {
	 if ( double_cuts.size() > _priority )
		double_cuts.resize( _priority );
  }

  return double_cuts;
}

void mign_cuts_paged::enumerate_node_with_bitsets( mign_node n, const std::vector<mign_node>& ns )
{
  for ( const auto& cut : enumerate_local_cuts( ns, _top_index, n ) ) 
  {

    auto area = std::get<2>( cut );
    area.resize( n + 1 );
    area.set( n );
    const auto extra = get_extra( _levels[n].first - std::get<1>( cut ), static_cast<unsigned int>( area.count() ) );
    data.append_set( n, get_index_vector( std::get<0>( cut ) ), extra );
    cones.append_set( n, get_index_vector( area ) );

	if (_almost == 1)
	{
		auto new_cut = std::get<0>( cut ); 
	    std::vector<mign_node> leafs; 
	    for (auto y = 0; y < new_cut.size(); y++)
	    {
	     if (new_cut[y] == 1)
	  	  leafs.push_back(y); 
	    }
	    auto func = mign_simulate_cut( _mign, n, leafs ); 
	    const auto num_vars = tt_num_vars( func );
	
		
		unsigned flag = 0u; 
		std::string str; 
	
		boost::dynamic_bitset<> reminder(func.size(),0); 
		boost::dynamic_bitset<> inputs_c(num_vars,0); 
		auto is_maj_p = is_maj(func, leafs, _mign, n); 
		if (is_maj_p.first > 0) 
		{	
			flag = 1; 
			inputs_c = is_maj_p.second; 
			almost_info[n].push_back(std::make_tuple(func,is_maj_p.first - 1,reminder,inputs_c, str));
		}
		else if (flag == 0) // qua ci va 0
		{
			auto is_almost_maj_p = is_almost_maj(func); 
			if ((std::get<0>(is_almost_maj_p) > 0) && (flag == 0)) 
			{
				inputs_c = std::get<1>(is_almost_maj_p); 
				reminder = std::get<2>(is_almost_maj_p);
				auto statistics = std::make_shared<properties>();
				auto settings = std::make_shared<properties>();
				boost::optional<unsigned> timeout = 5; 
				signed stop = static_cast<unsigned int>( area.count() ) - 3 - leafs.size(); 
				if (stop <= 0) 
				{
					almost_info[n].push_back(std::make_tuple(func,-1,reminder,inputs_c, str));
				}
				
				else 
				{
					settings->set( "timeout", timeout);
					settings->set( "stop", (unsigned)stop);
					auto mig = exact_mig_with_sat ( reminder,settings,statistics);  
					if (( (bool)mig ) )
					{ 
						auto mign = mig_to_mign (*mig);
					    str = mign_to_string( mign, mign.outputs()[0].first,settings,statistics);	
					    almost_info[n].push_back(std::make_tuple(func,std::get<0>(is_almost_maj_p),reminder,inputs_c, str));
						assert (mign.size() - mign.inputs().size() - 1 <= stop); 
							//}
					}
					else 
						almost_info[n].push_back(std::make_tuple(func,-1,reminder,inputs_c, str));
				}
				
			}
			else 
			{
				almost_info[n].push_back(std::make_tuple(func,-1,reminder,inputs_c, str)); 
			}
		}
		else 
		{
			almost_info[n].push_back(std::make_tuple(func,-1,reminder,inputs_c, str)); 
		}
		
	}
    
  }

}

std::vector<unsigned> mign_cuts_paged::get_extra( unsigned depth, unsigned size ) const
{
  std::vector<unsigned> v( 2u + _extra, 0u );
  v[0u] = depth;
  v[1u] = size;
  return v;
}

}

// Local Variables:
// c-basic-offset: 2
// eval: (c-set-offset 'substatement-open 0)
// eval: (c-set-offset 'innamespace 0)
// End:
