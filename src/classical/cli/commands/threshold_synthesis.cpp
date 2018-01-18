/* CirKit: A circuit toolkit
 * Copyright (C) 2009-2015  University of Bremen
 * Copyright (C) 2015-2016  EPFL
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

#include "threshold_synthesis.hpp"

#include <iostream>

#include <core/utils/program_options.hpp>

#include <classical/cli/stores_mign.hpp>
#include <classical/mign/mign.hpp>
#include <classical/mign/threshold_synthesis.hpp>
#include <classical/mign/mign_multiobj_opt.hpp>
#include <classical/mign/mign_flow_map.hpp>
#include <classical/mign/mign_flow_almost_maj.hpp>
#include <classical/mign/mign_cover.hpp> 
#include <classical/mign/mign_utils.hpp>
#include <classical/mign/mign_rewrite.hpp>
#include <cli/stores.hpp>
#include <classical/mign/mign_io.hpp>

#include <classical/abc/abc_api.hpp>
#include <classical/abc/abc_manager.hpp>

namespace abc
{
int Abc_ExorcismMain( Vec_Wec_t * vEsop, int nIns, int nOuts, char * pFileNameOut, int Quality, int Verbosity );
int Extra_ThreshHeuristic( word * t, int nVars, int * pW );
int Extra_ThreshCheck( word * t, int nVars, int * pW );
}

using boost::program_options::value;
using boost::format;

namespace cirkit
{
	/*
thres_synth_command::thres_synth_command ( const environment::ptr& env )
	: cirkit_command( env, "Threshold Synthesis using Neutzling method (arg tt)" )
{
	
}

bool thres_synth_command::execute()
{
	   
		auto statistics = std::make_shared<properties>();
		auto settings = std::make_shared<properties>();
		
        const auto& tts = env->store<tt>();
		
        if ( tts.current_index() == -1 )
        {
          std::cout << "[w] no truth tables in store" << std::endl;
          return true;
        }

        const auto num_vars = tt_num_vars( tts.current() );

		int weights[num_vars];
        std::vector<unsigned long> truth( tts.current().num_blocks() );
        boost::to_block_range( tts.current(), &truth[0] );
     

        abc_manager::get();
        abc::Abc_FrameGetGlobalFrame();
        auto threshold = abc::Extra_ThreshCheck( &truth[0], num_vars, weights );

        std::cout << "[i] Threshold = " << threshold << std::endl;
		
        if ( threshold != 0 )
        {
          std::cout << "[i] Weights =";
          for ( auto i = 0u; i < num_vars; ++i )
          {
            std::cout << " " << weights[i];
          }
          std::cout << std::endl;
		  
	   	  auto& migns = env->store<mign_graph>();
 	
	   	  auto mign = threshold_synthesis(tts.current(),settings, statistics);  
	   	  migns.extend(); 
	   	  migns.current() = mign;
        }
    	 
		else 
		{
			std::cout << " This is not a threshold" << std::endl; 
		}
   	  
	
	return true;
}*/

thres_majn_command::thres_majn_command ( const environment::ptr& env )
	: cirkit_command( env, "Generate Maj-n using from threshold using cut pruning - depth optimization" )
{
	opts.add_options()
		( "cut_size,c",          value_with_default( &cut_size ),       "cut size (maximum and default value 6)" )
	    ( "priority_cut,p",      value_with_default( &priority_cut ),   "priority cut (value, default = 0)" )
		( "allow_almost,a",      value_with_default( &allow_almost ),   "allow usage of almost majority -- not with multiobjective" )
		( "multi_obj_opt,m",     value_with_default( &multi_obj_opt ),  "multi objective optimization (both delay and energy) (1) ");
}

bool thres_majn_command::execute()
{
	auto statistics = std::make_shared<properties>();
	auto settings = std::make_shared<properties>();
	
	if (multi_obj_opt == false)
	{
		if (allow_almost == 0)
		{
		    clock_t t1,t2,t3,t4;
	        std::ostream& os = std::cout;
		
	        auto& migns = env->store<mign_graph>();
       
	        if ( migns.current_index() == -1 )
	        {
	           std::cout << "[w] no MIGn in store" << std::endl;
	           return true;
	        }
		
	   	    auto& mign = migns.current();
	  
	   	    settings->set( "cut_size", cut_size );
		    settings->set( "priority_cut", priority_cut); 
			settings->set( "almost_maj", allow_almost); 
	   	    t1 = clock();
	   	    mign_flow_map( mign, settings, statistics );
	  

	   	   t3 = clock();
	   	   auto mign_new = mign_cover_write (mign); 
	   	   t4 = clock();
	   	   float diff ((float)t4 - (float)t3);
	   	   float seconds = diff / CLOCKS_PER_SEC;
	   	   os << format( "[i] run-time writing graph MIG-n: %.2f seconds" ) % seconds << std::endl;
	 
	   	   auto mign_n = mign_rewrite_top_down( mign_new, settings,statistics ); 
	   	   t2 = clock();
	 
	   	   float diff_two ((float)t2 - (float)t1);
	   	   seconds = diff_two / CLOCKS_PER_SEC;
	   	   os << format( "[i] TOTAL run-time: %.2f seconds" ) % seconds << std::endl;
	 
	   	   migns.extend();
	       migns.current() = mign_n; 
		}
	    else 
		{
		    clock_t t1,t2,t3,t4;
	        std::ostream& os = std::cout;
		
	        auto& migns = env->store<mign_graph>();
       
	        if ( migns.current_index() == -1 )
	        {
	           std::cout << "[w] no MIGn in store" << std::endl;
	           return true;
	        }
		
	   	    auto& mign = migns.current();
	  
	   	    settings->set( "cut_size", cut_size );
		    settings->set( "priority_cut", priority_cut); 
			settings->set( "allow_almost", allow_almost); 
	   	    t1 = clock();
	   	    mign_flow_almost_maj( mign, settings, statistics );
	  

	   	   t3 = clock();
	   	   auto mign_new = mign_cover_write_maj (mign); 
	   	   t4 = clock();
	   	   float diff ((float)t4 - (float)t3);
	   	   float seconds = diff / CLOCKS_PER_SEC;
	   	   os << format( "[i] run-time writing graph MIG-n: %.2f seconds" ) % seconds << std::endl;
	 
	   	   auto mign_n = mign_rewrite_top_down( mign_new, settings,statistics ); 
	   	   t2 = clock();
	 
	   	   float diff_two ((float)t2 - (float)t1);
	   	   seconds = diff_two / CLOCKS_PER_SEC;
	   	   os << format( "[i] TOTAL run-time: %.2f seconds" ) % seconds << std::endl;
	 
	   	   migns.extend();
	       migns.current() = mign_n; 
		}
		
	}
	
	else 
	{
   	  auto& migns = env->store<mign_graph>();
      
         if ( migns.current_index() == -1 )
         {
           std::cout << "[w] no MIGn in store" << std::endl;
           return true;
         }
		 
   	  auto& mign = migns.current();
	  
   	  settings->set( "cut_size", cut_size );
   	  mign_flow_map_opt( mign, settings, statistics );
	   
   	 auto mign_new = mign_cover_write_multi(mign); 
   	 std::vector<mign_graph> mign_n; 

     for( auto m : mign_new)
   	 {
   		 mign_n.push_back(mign_rewrite_top_down( m, settings,statistics )); 
   	 }
	 
	 
   	 mign_n = simplify_vector_mign(mign_n); 
   	 mign_n.push_back(mign); // maybe also the original one is good ;) 
	 
   	 auto i = 0; 
   	 for( auto m : mign_n)
   	 {
   		 std::cout << "[" << i << "]" << " MIGn: "; 
   		 auto energy = evaluate_energy (m); 
   		 auto f = boost::lexical_cast<std::string>(energy);
   		 std::cout << " energy = " << energy; 
   		 auto depth = evaluate_depth_th (m); 
   		 auto g = boost::lexical_cast<std::string>(depth);
   		 auto output_verilog = "_multi_obj.v"; 
   		 auto filename = g + f + output_verilog;  
   		 std::cout << " depth = " << depth << std::endl; 
   		 write_verilog_compact( m, filename );
   		 migns.extend();
   		 migns.current() = m;
   		 i++; 
   	 }
	}
		
	  return true;
}

}

// Local Variables:
// c-basic-offset: 2
// eval: (c-set-offset 'substatement-open 0)
// eval: (c-set-offset 'innamespace 0)
// End:
