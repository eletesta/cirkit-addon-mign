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

#ifndef CLI_THRESHOLD_SYNTHESIS_COMMAND_HPP
#define CLI_THRESHOLD_SYNTHESIS_COMMAND_HPP

#include <cli/cirkit_command.hpp>

namespace cirkit
{

/*class thres_synth_command : public cirkit_command
{
public:
	thres_synth_command( const environment::ptr& env );
	
protected:
	bool execute();
};*/

class thres_majn_command : public cirkit_command
{
public:
	thres_majn_command( const environment::ptr& env );
	
protected:
	bool execute();

private:
	unsigned cut_size = 6u;
	unsigned priority_cut = 0u; 
	unsigned allow_almost = 0u;
	bool multi_obj_opt = false; 
};

}

#endif

// Local Variables:
// c-basic-offset: 2
// eval: (c-set-offset 'substatement-open 0)
// eval: (c-set-offset 'innamespace 0)
// End:
