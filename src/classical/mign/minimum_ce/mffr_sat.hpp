#ifndef MFFR_SAT
#define MFFR_SAT

#include <fstream>
#include <iostream>
#include <string>

#include <core/properties.hpp>
#include <classical/mign/mign.hpp>

namespace cirkit
{
   mign_graph mffrc_sat (mign_graph& mign, const properties::ptr& settings = properties::ptr(),const properties::ptr& statistics = properties::ptr() );    
}

#endif