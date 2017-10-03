#include "minimum_ce.hpp"

#include <cmath>
#include <vector>

#include <boost/assign/std/vector.hpp>
#include <boost/format.hpp>

#include <core/utils/timer.hpp>
#include <classical/mign/mign_utils.hpp>

#ifdef ADDON_FORMAL
#include <formal/utils/z3_utils.hpp>
#include <z3++.h>
#endif

using namespace boost::assign;

using boost::format;
using boost::str;
//using boost::numeric::ublas; 

namespace cirkit
{

/******************************************************************************
 * Types                                                                      *
 ******************************************************************************/
#ifdef ADDON_FORMAL
struct exact_mig_instance_min
{
  exact_mig_instance_min( unsigned num_level, unsigned num_vars, unsigned nr) :
  num_level (num_level), 
    num_vars( num_vars ),
	nr(nr),
    solver( ctx ),
    out( 1u << num_vars ),
    in( 3u, std::vector<std::vector<z3::expr>>( 1u << num_vars ) ),
    sel( 3u ),
	neg( 3u ),
	s (nr),
    out_neg( ctx.bool_const( "out_neg" ))
  {}

	  void add_level()
	    {
	      const auto bw = (unsigned)ceil( log( 2u + num_vars + 7u ) / log( 2u ) );

	      unsigned level = sel[0u].size();
	        //std::cout<<"valore di level in quello buono" << level << std::endl;

	      for ( auto x = 0u; x < 3u; ++x )
	      {
	        sel[x] += ctx.bv_const( str( format( "sel%d_%d" ) % x % level ).c_str(), bw );
	        neg[x] += ctx.bool_const( str( format( "neg%d_%d" ) % x % level ).c_str() );

	        /* assertion for sel[level] <= n + k */
	        solver.add( z3::ule( sel[x][level], ctx.bv_val( num_vars + level, bw ) ) );
	      }

	     // solver.add( sel[0][level] < sel[1][level] );
	      //solver.add( sel[1][level] < sel[2][level] );

	      for ( auto i = 0u; i < level; ++i )
	      {
	        solver.add( sel[0][level] != sel[0][i] || sel[1][level] != sel[1][i] || sel[2][level] != sel[2][i] ||
	                    neg[0][level] != neg[0][i] || neg[1][level] != neg[1][i] || neg[2][level] != neg[2][i] );
	      }

	      for ( auto j = 0u; j < ( 1u << num_vars ); ++j )
	      {
	        out[j] += ctx.bool_const( str( format( "out_%d_%d" ) % j % level ).c_str() );
	        in[0u][j] += ctx.bool_const( str( format( "in1_%d_%d" ) % j % level ).c_str() );
	        in[1u][j] += ctx.bool_const( str( format( "in2_%d_%d" ) % j % level ).c_str() );
	        in[2u][j] += ctx.bool_const( str( format( "in3_%d_%d" ) % j % level ).c_str() );

	        /* assertion for out[j][level] = M(in1[j][level],in2[j][level],in3[j][level] */
	        solver.add( out[j][level] == ( ( in[0u][j][level] && in[1u][j][level] ) || ( in[0u][j][level] && in[2u][j][level] ) || ( in[1u][j][level] && in[2u][j][level] ) ) );

	        /* assertions for in[x][j][level] = neg[level] ^ ite( sel[level], ... ) */
	        boost::dynamic_bitset<> val( num_vars, j );
	        for ( auto x = 0u; x < 3u; ++x )
	        {
	          solver.add( implies( sel[x][level] == ctx.bv_val( 0u, bw ),
	                               in[x][j][level] == ( neg[x][level] != ctx.bool_val( false ) ) ) );

	          for ( auto l = 0u; l < num_vars; ++l )
	          {
	            solver.add( implies( sel[x][level] == ctx.bv_val( l + 1u, bw ),
	                                 in[x][j][level] == ( neg[x][level] != ctx.bool_val( val[l] ) ) ) );
	          }
	          for ( auto l = 0u; l < level; ++l )
	          {
	            solver.add( implies( sel[x][level] == ctx.bv_val( l + num_vars + 1u, bw ),
	                                 in[x][j][level] == ( neg[x][level] != out[j][l] ) ) );
	          }
	        }
	      }
	    }

	    void constrain( const tt& spec )
	    {
	      for ( auto j = 0u; j < ( 1u << num_vars ); ++j )
	      {
	        solver.add( (out[j].back() != out_neg) == ctx.bool_val( spec[j] ) );
	      }
	    }

 
  void extract_sel (mign_graph mig) 
 	  {
 		  unsigned n = mig.inputs().size(); 
 		  //unsigned k = num_vertices(mig) - n - 1u;
 		  unsigned count = 0u; 
 		  const auto bw = (unsigned)ceil( log( 2u + num_vars + 7u ) / log( 2u ) );
		
		  
 		for ( unsigned node=(n+1); node< mig.num_gates();node++)
 		{
			 
 		    const auto children = mig.children(node);
 			//std::cout << node << std::endl;
	
   	      for ( auto x = 0u; x < 3u; ++x )
   	      {
 			  unsigned i = 0u; 
 			  for (const auto& nodemig : mig.vertices()//boost::make_iterator_range (boost::vertices(mig)))
 			  {		
 				  if (children[x].node == nodemig) {
                      unsigned g = children[x].node;
 					  solver.add(sel[x][count] == ctx.bv_val(g,bw));
                      //std::cout << "count" << count <<std::endl;
                      //std::cout << "x" << x<< "=" << i << "e invece children[x].node ="<< children[x].node << std::endl;
 				  }
 				  else 
 					  i++; 
 			  }

 		  }
			  
 			
 	  count++;
   } 
   } 
  

    void sum(unsigned r, unsigned level, mign_graph mig, unsigned count_zero)
    {
       // unsigned n = level * 3 + 1; // utneg e i 3 input di ogni nodo 
	  const unsigned inputs = mig.inputs().size(); 
	 
					  
	  unsigned n = level *3 + 1 - count_zero; 
        for ( auto x = 1u; x <= (n+1); ++x)
        {
            var += ctx.bool_const( str( format( "var%d" ) % x ).c_str() );
		}
		
		auto x = 1u; 
		auto i = 0u; 
			
	 		for ( unsigned node=(input+1); node< mig.num_gates();node++)
	 		{
				auto k = x+1; 
				auto j = x+2;
				const auto children = mig.children(node);
	 			//std::cout << node << std::endl;
				if ((children[0u].node == 0u) && (children[1u].node != 0u) && (children[2u].node != 0u))
				{
            // assertion to say that x1 x2 and x3 sono rispettivamente neg[0][i] neg[1][i] etc 
           solver.add((var[x] == neg[1][i]) && (var[k] == neg[2][i]));
		   x = x+ 2; 
		   ++i; 
	   }
		   
		   else if ((children[0u].node != 0u) && (children[1u].node == 0u) && (children[2u].node != 0u))
				{
            // assertion to say that x1 x2 and x3 sono rispettivamente neg[0][i] neg[1][i] etc 
           solver.add((var[x] == neg[0][i]) && (var[k] == neg[2][i]));
		   x = x+ 2; 
		   ++i; 
	   }
		   
		else if ((children[0u].node != 0u) && (children[1u].node != 0u) && (children[2u].node == 0u))
		{
       // assertion to say that x1 x2 and x3 sono rispettivamente neg[0][i] neg[1][i] etc 
      solver.add((var[x] == neg[0][i]) && (var[k] == neg[1][i]));
   x = x+ 2; 
   ++i; 
  }
	  
	else if ((children[0u].node != 0u) && (children[1u].node != 0u) && (children[2u].node != 0u))
	{
  // assertion to say that x1 x2 and x3 sono rispettivamente neg[0][i] neg[1][i] etc 
 solver.add((var[x] == neg[0][i]) && (var[k] == neg[1][i]) && (var[j] == neg[2][i]));
		   x = x+3; 
		   ++i; 
	   }
   }
	   solver.add(var[x] == out_neg); 
	  // return count_zero; 
	  // std::cout << "fine di sum" << std::endl; 
   }
   
   
	
	void define_s (unsigned r, unsigned n, unsigned kappa)
	{
		// unsigned n = level * 3 + 1;
		// std::cout << "n =" << n << std::endl; 
		//  std::cout << "r =" << r << std::endl; 
	unsigned k = kappa; 
		//std::cout << k << std::endl;
		for (unsigned j = 1u; j <= (n-r); ++j) {
			 //std::cout << "ciao4"<< std::endl; 
			 s[j] += ctx.bool_const( str( format( "s%d_%d" ) % j % k).c_str());
		}	
	}
	
	void define_alls (unsigned r, unsigned n, unsigned kappa) // rand n names are the same as ex 30 and 31 libro 
	{
		// unsigned n = level * 3 + 1;
		unsigned k = kappa; 
		//std::cout << "valore di r cioe p " << r<< std::endl; 
		//std::cout << "valore di n " << n<< std::endl; 
		if (k == 0u) {
				//std::cout << "ciao4 siamo nel k = 0"<< std::endl; 
			for (unsigned j = 1u; j <= (n-r); ++j) {
				 solver.add((!var[j+k]) || (s[j][k+1]) ); 
			}
		}
		
		else if (k == r)
		{
			auto j = 0u; 
			for (j = 1u; j < (n-r); ++j) {
				//std::cout << "valore di j =" << j <<std::endl; 
				//std::cout << "ciao4 siamo nel k = r 1"<< std::endl; 
			 solver.add (!var[j+k] || !s[j][k]); 
			// std::cout << "ciao4 siamo nel k = r2"<< std::endl; 
			 solver.add ( !s[j][k] || s[j+1][k] ); 
			  //solver.add (var[j+k] || s[j][k] || !s[j+1][k]);
		}
		//std::cout<< "valore di k = " << k << std::endl << "valore di j =" << j <<std::endl; 
		 solver.add (!var[j+k] || !s[j][k]);
		// std::cout << "fine k = r"<< std::endl; 
	}
			
		else {	
		for (auto j = 1u; j <= (n-r); ++j) {
			 solver.add (!var[j+k] || !s[j][k] || s[j][k+1]); 
			 // solver.add (s[j][k] || !s[j][k+1]);
		}
		for (auto j = 1u; j < (n-r); ++j) {
		solver.add ( !s[j][k] || s[j+1][k] ); 
		//solver.add (var[j+k] || s[j][k] || !s[j+1][k]); 
	}
	}
}
		

  mign_graph extract_solution( const std::string& model_name, const std::string& output_name, bool verbose )
  {
    mign_graph mig;
    mig.set_name(model_name );

    std::vector<mig_function> inputs, nodes;
    for ( auto i = 0u; i < num_vars; ++i )
    {
      inputs += mig.create_pi(str( format( "x%d" ) % i ) );
    }

    const auto m = solver.get_model();

    for ( auto i = 0u; i < sel[0].size(); ++i )
    {
      mig_function children[3];

      for ( auto x = 0u; x < 3u; ++x )
      {
        auto sel_val = to_bitset( m.eval( sel[x][i] ) ).to_ulong();

        if ( sel_val == 0u )
        {
          children[x] = mig.get_constant(false );
        }
        else if ( sel_val > 0 && sel_val <= num_vars )
        {
          children[x] = inputs[sel_val - 1u];
        }
        else
        {
          children[x] = nodes[sel_val - num_vars - 1u];
        }

        if ( expr_to_bool( m.eval( neg[x][i] ) ) )
        {
          children[x] = !children[x];
        }
      }

      nodes += mig.create_maj( children);

      if ( verbose )
      {
        std::cout << format( "added node (%d,%d) for i = %d" ) % nodes.back().node % nodes.back().complemented % i << std::endl;
        for ( auto x = 0u; x < 3u; ++x )
        {
          std::cout << format( "  - child %d = (%d,%d)" ) % x % children[x].node % children[x].complemented << std::endl;
        }
      }
    }

    auto f = nodes.back();

    if ( expr_to_bool( m.eval( out_neg ) ) )
    {
      f = !f;
      }

    mig.create_po(f, output_name );

    return mig;
  }

  unsigned num_level; 
  unsigned num_vars;
   unsigned nr; 
  z3::context ctx;
  z3::solver solver;

  std::vector<std::vector<z3::expr>> out;
  std::vector<std::vector<std::vector<z3::expr>>> in;

  std::vector<std::vector<z3::expr>> sel;
  std::vector<std::vector<z3::expr>> neg;

    std::vector<std::vector<z3::expr>> s;
    
  z3::expr out_neg;
  
  std::vector<z3::expr> var;
  
  
};

/******************************************************************************
 * Private functions                                                          *
 ******************************************************************************/

mign_graph minimum_ce(mign_graph& mig, const tt& spec, unsigned start, const std::string& model_name,const std::string&  output_name, bool verbose )
{
  unsigned p = start;
  
  //const auto& info = mig_info(mig); 
  unsigned ninput = mig.inputs.size(); // numero di inputs
  unsigned k = mig.num_gates() - ninput - 1u ; // numero di nodi 
  
  auto count_zero = 0u; 
	 
	for ( unsigned node=(ninput+1); node< mig.num_gates();node++)
	{
		 
	    const auto children = mig.children(node); 


		for ( auto x = 0u; x < 3u; ++x )
		 {	
			  if (children[x].node == 0u) {
				  ++count_zero; 
			  }
	  }
  }
   
  while ( p > 0)
  {
    if ( verbose )
    {
      std::cout << boost::format( "[i] check for realization with lte %d ce" ) % p << std::endl;
    }

	auto n = k*3 + 2 - count_zero; 
	exact_mig_instance_min inst3( k, ninput, n+p);
	     
		  for ( auto i = 0u; i < k; ++i )
	      {
			   
	        inst3.add_level();
			
	      }
		 inst3.constrain(spec);
		 inst3.extract_sel(mig);
		  
		 inst3.sum(p,k,mig,count_zero);
		
		 
		  for ( auto kappa = 0u; kappa <= p; ++kappa )
		   {
			     inst3.define_s(p,(n-1),kappa);
				  }
	  		  
				  for ( auto kappa = 0u; kappa <= p; ++kappa )
				   {
					   inst3.define_alls(p,(n-1),kappa);
				   }
				 
	
          if ( inst3.solver.check() == z3::sat )
          {
              mig = inst3.extract_solution( model_name, output_name, verbose );
			  auto ce = compute_ce(mig); 
			  --p; 
			  if ( ce < p) 
				  p = ce; 
          }
		  else
		  { 
			  //std::cout << "it is impossible to realize the MIG with" << p << " ce" << std::endl;
			  return mig; 
		  }
	  }
	  return mig; 
  }
  
mign_graph exact_mig_with_sat_explicit( const tt& spec, unsigned start, const std::string& model_name, const std::string& output_name, bool verbose)
{
    auto p = 0u;
    auto k = start;
    while ( true )
    {
      if ( verbose ) {
          if (p == 0u)
          {
           std::cout << boost::format( "[i] check for realization with %d gates" ) % k << std::endl;
          }
          else
          std::cout << boost::format( "[i] check for realization with %d gates and %d ce" ) % (k-1) %p << std::endl;
          }
     
	 exact_mig_instance_min inst( tt_num_vars( spec ),0,0);
	
	  if ( p == 0u) {
          for ( auto i = 0u; i < k; ++i )
          {
           inst.add_level();
          }

          inst.constrain( spec );

          if ( inst.solver.check() == z3::sat )
              {
              ++p;
              }
 
         ++k; 
         }
      
        else {
  		  const auto g = k-1; 
  		  auto n = g*3 + 2 - p; 
  		  exact_mig_instance_min inst2( tt_num_vars( spec ),n,0);
  	      for ( auto i = 0u; i < (k-1); ++i )
  	         {
  	             inst2.add_level();
  	         }

  	         inst2.constrain( spec );
             inst2.sum(p,g);
		 
          for ( auto kappa = 0u; kappa <= p; ++kappa )
             {
             inst2.define_s(p,g,kappa);
             }
  	 
	  		  
       	for ( auto kappa = 0u; kappa <= p; ++kappa )
  	        {
  	         inst2.define_alls(p,g,kappa);
  	        }
        if ( inst2.solver.check() == z3::sat )
            {
                return inst2.extract_solution( model_name, output_name, verbose );
            }
            ++p;
            }
      }
}

  
#endif



/******************************************************************************
 * Public functions                                                           *
 ******************************************************************************/
  
 mign_graph min_ce_with_sat( mign_graph mig, const tt& spec,
                                const properties::ptr& settings,
                                const properties::ptr& statistics)
  {
	
	auto ce = compute_ce(mig);
    /* settings */
  	const auto start         = get( settings, "start",       ce );
    const auto model_name    = get( settings, "model_name",  std::string( "minimum_ce" ) );
  
    const auto output_name   = get( settings, "output_name", std::string( "f" ) );
    const auto minimum_after = get( settings, "minimum_after", true );  // minimize after or together with exact synthesis
    const auto verbose       = get( settings, "verbose",     false );


    /* timing */
    properties_timer t( statistics );
	
#ifdef ADDON_FORMAL

    if ( minimum_after )
    {
  	   return minimum_ce( mig,spec, start, model_name, output_name,  verbose );
       std::cout << format( "[i] run-time: %.2f seconds" ) % statistics->get<double>( "runtime" ) << std::endl;
    }
    else 
  	  return exact_mig_with_sat_explicit( spec, start, model_name, coutput_name,verbose);
      
#else
  std::cout << "[e] z3 solver requires formal addon enabled" << std::endl;
  return mig; 
#endif
  }

}
