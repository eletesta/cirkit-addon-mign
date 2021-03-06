#include "mig_strash.hpp"

#include "classical/mig/mig_rewrite.hpp"

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

mig_graph mig_strash( const mig_graph& mig )
{
  return mig_rewrite_top_down( mig, mig_rewrite_default_maj );
}

}