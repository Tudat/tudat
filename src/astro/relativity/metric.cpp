/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <boost/make_shared.hpp>

#include "tudat/astro/relativity/metric.h"

namespace tudat
{

namespace relativity
{

//! Initialize global PPN parameters
std::shared_ptr< PPNParameterSet > ppnParameterSet = std::make_shared< PPNParameterSet >( 1.0, 1.0 );

double equivalencePrincipleLpiViolationParameter = 0.0;

}

}
