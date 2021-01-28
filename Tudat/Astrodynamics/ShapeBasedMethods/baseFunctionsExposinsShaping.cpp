/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */


#include "baseFunctionsExposinsShaping.h"
#include <math.h>


namespace tudat
{
namespace shape_based_methods
{

std::shared_ptr< BaseFunctionExposinsShaping > createBaseFunctionExposinsShaping(
        const baseFunctionExposinsShapingType baseFunctionType )
{
    // Declare return object.
    std::shared_ptr< BaseFunctionExposinsShaping > baseFunctionExposinsShaping;

    // Check which type of base function is to be created.
    switch( baseFunctionType )
    {
    case constantExposinsShaping:
    {
        baseFunctionExposinsShaping = std::make_shared< ConstantFunctionExposinsShaping >( );
        break;
    }
    case linearExposinsShaping:
    {
        baseFunctionExposinsShaping = std::make_shared< LinearFunctionExposinsShaping >( );
        break;
    }
    case squaredExposinsShaping:
    {
        baseFunctionExposinsShaping = std::make_shared< SquaredFunctionExposinsShaping >( );
        break;
    }
    case cosineExposinsShaping:
    {
        baseFunctionExposinsShaping = std::make_shared< CosineFunctionExposinsShaping >( );
        break;
    }
    case powerCosineExposinsShaping:
    {
        baseFunctionExposinsShaping = std::make_shared< PowerCosineFunctionExposinsShaping >( );
        break;
    }
    case sineExposinsShaping:
    {
        baseFunctionExposinsShaping = std::make_shared< SineFunctionExposinsShaping >( );
        break;
    }
    case powerSineExposinsShaping:
    {
        baseFunctionExposinsShaping = std::make_shared< PowerSineFunctionExposinsShaping >( );
        break;
    }
    default:
    {
        throw std::runtime_error(
                    "Error, did not recognize base function type for exposins shaping" );
    }

    }
    return baseFunctionExposinsShaping;
}

} // namespace shape_based_methods
} // namespace tudat
