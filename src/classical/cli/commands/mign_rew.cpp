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

#include "mign_rew.hpp"

#include <iostream>

#include <core/utils/program_options.hpp>

#include <classical/cli/stores_mign.hpp>
#include <classical/mign/mign.hpp>
#include <classical/mign/threshold_synthesis.hpp>
#include <classical/mign/mign_multiobj_opt.hpp>
#include <classical/mign/mign_flow_map.hpp>
#include <classical/mign/mign_cover.hpp> 
#include <classical/mign/mign_utils.hpp>
#include <classical/mign/mign_rewrite.hpp>
#include <classical/cli/stores.hpp>
#include <classical/mign/mign_io.hpp>
#include <classical/mign/mign_lut_based_synthesis.hpp> 
#include <classical/mign/mign_rewriting.hpp>
#include <classical/mign/mign_rewriting_majn_to_smaller.hpp>
#include <classical/mign/mign_to_mig3.hpp>
#include <classical/mign/mign_invfree.hpp>
#include <classical/mign/create_homo.hpp>

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

rules_rewriting_command::rules_rewriting_command ( const environment::ptr& env )
	: cirkit_command( env, "Rewrite Mig-n using Associativity and Distributivity Rules" )
{
	
}

bool rules_rewriting_command::execute()
{
	   
    auto statistics = std::make_shared<properties>();
	auto& migns = env->store<mign_graph>();
    auto settings = std::make_shared<properties>();
	
    if ( migns.current_index() == -1 )
    {
      std::cout << "[w] no MIGn in store" << std::endl;
      return true;
    }
 
    auto& mign_old = migns.current();

    auto mign_rew = mign_rewriting(mign_old,settings,statistics);

    migns.extend();
	migns.current() = mign_rew; 
	
	return true;
}

homo_rewriting_command::homo_rewriting_command ( const environment::ptr& env )
	: cirkit_command( env, "Write homogeneous MIGn - not always possible" )
{
	opts.add_options()
		( "n_inputs,n",          value_with_default( &n_inputs ),       "number of inputs (0 by default)" );
}

bool homo_rewriting_command::execute()
{
    auto statistics = std::make_shared<properties>();
	auto settings = std::make_shared<properties>();
	
	auto Nin = 0u; 
	
	auto& migns = env->store<mign_graph>();
    
    if ( migns.current_index() == -1 )
    {
      std::cout << "[w] no MIGn in store" << std::endl;
      return true;
    }
	
	if (n_inputs == 0)
	{
        std::cout << " Insert the number of inputs with option -n " << std::endl;
        return true;
	}
	
	auto& mign_old = migns.current();
   
	auto mign_new = create_homo(mign_old,n_inputs);
	
    migns.extend();
    migns.current() = mign_new; 
		
	  return true;
}

reduce_n_inputs_command::reduce_n_inputs_command ( const environment::ptr& env )
	: cirkit_command( env, "Rewrite to smaller n (number different from 3 are allowed - whatever number of n inputs)" )
{
	opts.add_options()
		( "n_max,n",          value_with_default( &max_inputs ),       "number of max inputs in the circuit (default = 3)" );
}

bool reduce_n_inputs_command::execute()
{
    auto statistics = std::make_shared<properties>();

    auto& migns = env->store<mign_graph>();
    auto settings = std::make_shared<properties>();
       
	if ( migns.current_index() == -1 )
    {
       std::cout << "[w] no MIGn in store" << std::endl;
        return true;
     }
 
    auto& mign_old = migns.current();
 
    auto mign_rew = mign_rewriting_majn_to_smaller(mign_old, max_inputs,settings,statistics ); 
 
     migns.extend();
   	 migns.current() = mign_rew; 
		
    return true;
}

to_mig3_command::to_mig3_command ( const environment::ptr& env )
	: cirkit_command( env, "Write maj-5 and maj-7 into maj-3" )
{

}

bool to_mig3_command::execute()
{
	auto& migns = env->store<mign_graph>();
	
    if ( migns.current_index() == -1 )
    {
      std::cout << "[w] no MIGn in store" << std::endl;
      return true;
    }
	
	auto& mign_old = migns.current();

    auto mign_new = mign_to_mig3(mign_old);
	migns.extend();
 	migns.current() = mign_new; 
		
    return true;
}

luts_mign_command::luts_mign_command ( const environment::ptr& env )
	: cirkit_command( env, "Create an MIg-n network using LUTs approach - input blif file" )
{
	opts.add_options()
		( "filename,f",          value_with_default( &filename ),       "file name (mandatory)" );
}

bool luts_mign_command::execute()
{
	 auto statistics = std::make_shared<properties>();

     auto& migns = env->store<mign_graph>();
     auto settings = std::make_shared<properties>();
 
     auto mign = mign_lut_based_synthesis(filename, settings, statistics ); 
	  
     migns.extend();
     migns.current() = mign; 
		
     return true;
}

mign_inv_free_command::mign_inv_free_command ( const environment::ptr& env )
	: cirkit_command( env, "Create an MIg-n with all inversions on Primary Inputs" )
{
	
}

bool mign_inv_free_command::execute()
{
    auto& migns = env->store<mign_graph>();
   
    if ( migns.current_index() == -1 )
    {
      std::cout << "[w] no MIGn in store" << std::endl;
      return true;
    }

    auto& mign = migns.current();
  
    auto mign_new = mign_invfree(mign);
    migns.extend();
  
    migns.current() = mign_new;
		
     return true;
}

mign_fo_restr_command::mign_fo_restr_command ( const environment::ptr& env )
	: cirkit_command( env, "Create an MIg-n network using LUTs approach - input blif file" )
{
	opts.add_options()
		( "depth_constant,d",          value_with_default( &depth_const ),       "constant depth (1) or not (0). Default= 1" )
		( "max_fanout,f",              value_with_default( &max_fanout ),        "maximum fan-out in the network (default = 3)" );
}

bool mign_fo_restr_command::execute()
{
    auto& migns = env->store<mign_graph>();
   
    if ( migns.current_index() == -1 )
    {
      std::cout << "[w] no MIGn in store" << std::endl;
      return true;
    }

    auto& mign = migns.current();
  
    auto mign_new = mign_fo_restricted_opt( mign,max_fanout,depth_const); 
    migns.extend();
  
    migns.current() = mign_new;
		
    return true;
}

}

// Local Variables:
// c-basic-offset: 2
// eval: (c-set-offset 'substatement-open 0)
// eval: (c-set-offset 'innamespace 0)
// End:
