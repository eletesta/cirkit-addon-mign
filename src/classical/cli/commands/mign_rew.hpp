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

/**
 * @file exact_mign.hpp
 *
 * @brief TODO
 *
 * @author Eleonora
 * @since  2.3
 */

#ifndef CLI_MIGN_REWRITING_COMMAND_HPP
#define CLI_MIGN_REWRITING_COMMAND_HPP

#include <cli/cirkit_command.hpp>
#include <classical/utils/truth_table_utils.hpp>

#include <cudd.h>
#include <cuddInt.h>

namespace cirkit
{

class rules_rewriting_command : public cirkit_command
{
public:
	rules_rewriting_command( const environment::ptr& env );
	
protected:
	bool execute();
};

class homo_rewriting_command : public cirkit_command
{
public:
	homo_rewriting_command( const environment::ptr& env );
	
protected:
	bool execute();

private:
	unsigned n_inputs = 0;
};

class reduce_n_inputs_command : public cirkit_command
{
public:
	reduce_n_inputs_command( const environment::ptr& env );
	
protected:
	bool execute();
private: 
	unsigned max_inputs = 3u; 
};

class to_mig3_command : public cirkit_command
{
public:
	to_mig3_command( const environment::ptr& env );
	
protected:
	bool execute();
};

class luts_mign_command : public cirkit_command
{
public:
	luts_mign_command( const environment::ptr& env );
	
protected:
	bool execute();
private:
	std::string filename;
};

class mign_inv_free_command : public cirkit_command
{
public:
	mign_inv_free_command( const environment::ptr& env );
	
protected:
	bool execute();
};

class mign_fo_restr_command : public cirkit_command
{
public:
	mign_fo_restr_command( const environment::ptr& env );
	
protected:
	bool execute();
private:
	bool depth_const = 1; 
	unsigned max_fanout = 3; 
};

class minim_ce_command : public cirkit_command
{
public:
	minim_ce_command( const environment::ptr& env );
	
protected:
	bool execute();
	
private:
	tt spec; 
};

class bdd_to_mign_command : public cirkit_command
{
public:
	bdd_to_mign_command( const environment::ptr& env );
	
protected:
	bool execute();
private:
	bool ce_on = 0; 
	unsigned order_option = 0; 
	unsigned inputs = 3; 
};

}

#endif

// Local Variables:
// c-basic-offset: 2
// eval: (c-set-offset 'substatement-open 0)
// eval: (c-set-offset 'innamespace 0)
// End:
