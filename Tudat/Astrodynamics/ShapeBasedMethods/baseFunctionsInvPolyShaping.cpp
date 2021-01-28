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


#include "baseFunctionsInvPolyShaping.h"
#include <math.h>


namespace tudat
{
namespace shape_based_methods
{

std::shared_ptr< BaseFunctionInvPolyShaping > createBaseFunctionInvPolyShaping(
        const baseFunctionInvPolyShapingType baseFunctionType )
{
    // Declare return object.
    std::shared_ptr< BaseFunctionInvPolyShaping > baseFunctionInvPolyShaping;

    // Check which type of base function is to be created.
    switch( baseFunctionType )
    {
    case constantInvPolyShaping:
    {
        baseFunctionInvPolyShaping = std::make_shared< ConstantFunctionInvPolyShaping >( );
        break;
    }
    case linearInvPolyShaping:
    {
        baseFunctionInvPolyShaping = std::make_shared< LinearFunctionInvPolyShaping >( );
        break;
    }
    case squaredInvPolyShaping:
    {
        baseFunctionInvPolyShaping = std::make_shared< SquaredFunctionInvPolyShaping >( );
        break;
    }
    case cubedInvPolyShaping:
    {
        baseFunctionInvPolyShaping = std::make_shared< CubedFunctionInvPolyShaping >( );
        break;
    }
    case quarticInvPolyShaping:
    {
        baseFunctionInvPolyShaping = std::make_shared< QuarticFunctionInvPolyShaping >( );
        break;
    }
    case quinticInvPolyShaping:
    {
        baseFunctionInvPolyShaping = std::make_shared< QuinticFunctionInvPolyShaping >( );
        break;
    }
    case sexticInvPolyShaping:
    {
        baseFunctionInvPolyShaping = std::make_shared< SexticFunctionInvPolyShaping >( );
        break;
    }
    case septicInvPolyShaping:
    {
        baseFunctionInvPolyShaping = std::make_shared< SepticFunctionInvPolyShaping >( );
        break;
    }
    default:
    {
        throw std::runtime_error(
                    "Error, did not recognize base function type for InvPoly shaping" );
    }

    }
    return baseFunctionInvPolyShaping;
}

} // namespace shape_based_methods
} // namespace tudat
