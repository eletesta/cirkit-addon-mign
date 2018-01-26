#ifndef FFR_MIN_CE
#define FFR_MIN_CE

#include <fstream>
#include <iostream>
#include <string>

#include <core/properties.hpp>
#include <classical/mign/mign.hpp>

namespace cirkit
{
   mign_graph ffr_min_ce (mign_graph& mig, const properties::ptr& settings = properties::ptr(),const properties::ptr& statistics = properties::ptr() );    
}

#endif