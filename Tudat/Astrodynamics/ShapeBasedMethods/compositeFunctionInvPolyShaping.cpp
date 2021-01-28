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


#include "compositeFunctionInvPolyShaping.h"
#include <math.h>
#include <iostream>

namespace tudat
{
namespace shape_based_methods
{

void CompositeRadialFunctionInvPolyShaping::resetCompositeFunctionCoefficients( Eigen::VectorXd compositeFunctionCoefficients )
{
    // Check whether the size is correct.
    if( compositeFunctionCoefficients.size() == compositeFunctionComponents_.size() )
    {
        compositeFunctionCoefficients_ = compositeFunctionCoefficients;
    }
    else
    {
        std::cerr << "The number of composite function coefficients is incorrect!" << std::endl
             << compositeFunctionCoefficients.size() << " coefficients were put in, however, "
             << compositeFunctionComponents_.size() << " are required. The composite function coefficients are not reset!\n";
    }
}


double CompositeRadialFunctionInvPolyShaping::evaluateCompositeFunction( const double independentVariable )
{
    double functionValue = 0.0;
    for( unsigned int i = 0; i < compositeFunctionComponents_.size(); i++ )
    {
        functionValue += compositeFunctionCoefficients_[i] * compositeFunctionComponents_[ i ]->evaluateFunction( independentVariable );
    }

    return 1.0 / functionValue;
}


double CompositeRadialFunctionInvPolyShaping::evaluateCompositeFunctionFirstDerivative( const double independentVariable )
{
    double sumComponentsFirstDerivativeValue = 0.0;
    for( unsigned int i = 0; i < compositeFunctionComponents_.size(); i++ )
    {
        sumComponentsFirstDerivativeValue += compositeFunctionCoefficients_[i] * compositeFunctionComponents_[ i ]->evaluateFirstDerivative( independentVariable );
    }
    return - sumComponentsFirstDerivativeValue * std::pow( evaluateCompositeFunction( independentVariable ), 2.0 );
}


double CompositeRadialFunctionInvPolyShaping::evaluateCompositeFunctionSecondDerivative( const double independentVariable )
{
    double functionValue = 0.0;
    for( unsigned int i = 0; i < compositeFunctionComponents_.size(); i++ )
    {
        functionValue += compositeFunctionCoefficients_[i] * compositeFunctionComponents_[ i ]->evaluateFunction( independentVariable );
    }

    double sumComponentsSecondDerivativeValue = 0.0;
    for( unsigned int i = 0; i < compositeFunctionComponents_.size(); i++ )
    {
        sumComponentsSecondDerivativeValue += compositeFunctionCoefficients_[i] * compositeFunctionComponents_[ i ]->evaluateSecondDerivative( independentVariable );
    }
    return - sumComponentsSecondDerivativeValue * std::pow( evaluateCompositeFunction( independentVariable ), 2.0 )
            + 2.0 * functionValue * std::pow( evaluateCompositeFunctionFirstDerivative( independentVariable ), 2.0 );
}

double CompositeRadialFunctionInvPolyShaping::evaluateCompositeFunctionThirdDerivative( const double independentVariable )
{
    double functionValue = 0.0;
    for( unsigned int i = 0; i < compositeFunctionComponents_.size(); i++ )
    {
        functionValue += compositeFunctionCoefficients_[i] * compositeFunctionComponents_[ i ]->evaluateFunction( independentVariable );
    }

    double sumComponentsFirstDerivativeValue = 0.0;
    for( unsigned int i = 0; i < compositeFunctionComponents_.size(); i++ )
    {
        sumComponentsFirstDerivativeValue += compositeFunctionCoefficients_[i] * compositeFunctionComponents_[ i ]->evaluateFirstDerivative( independentVariable );
    }

    double sumComponentsSecondDerivativeValue = 0.0;
    for( unsigned int i = 0; i < compositeFunctionComponents_.size(); i++ )
    {
        sumComponentsSecondDerivativeValue += compositeFunctionCoefficients_[i] * compositeFunctionComponents_[ i ]->evaluateSecondDerivative( independentVariable );
    }

    double sumComponentsThirdDerivativeValue = 0.0;
    for( unsigned int i = 0; i < compositeFunctionComponents_.size(); i++ )
    {
        sumComponentsThirdDerivativeValue += compositeFunctionCoefficients_[i] * compositeFunctionComponents_[ i ]->evaluateThirdDerivative( independentVariable );
    }

    return - sumComponentsThirdDerivativeValue * std::pow( evaluateCompositeFunction( independentVariable ), 2.0 )
            - 2.0 * evaluateCompositeFunction( independentVariable ) * evaluateCompositeFunctionFirstDerivative( independentVariable )
            * sumComponentsSecondDerivativeValue
            + 2.0 * sumComponentsFirstDerivativeValue * std::pow( evaluateCompositeFunctionFirstDerivative( independentVariable ), 2.0 )
            + 4.0 * functionValue * evaluateCompositeFunctionFirstDerivative( independentVariable )
            * evaluateCompositeFunctionSecondDerivative( independentVariable );
}

double CompositeRadialFunctionInvPolyShaping::getComponentFunctionCurrentValue( int componentIndex, double currentTime )
{
    return compositeFunctionComponents_[ componentIndex ]->evaluateFunction( currentTime );
}


double CompositeRadialFunctionInvPolyShaping::getComponentFunctionFirstDerivative( int componentIndex, double currentTime )
{
    return compositeFunctionComponents_[ componentIndex ]->evaluateFirstDerivative( currentTime );
}

double CompositeRadialFunctionInvPolyShaping::getComponentFunctionSecondDerivative( int componentIndex, double currentTime )
{
    return compositeFunctionComponents_[ componentIndex ]->evaluateSecondDerivative( currentTime );
}

double CompositeRadialFunctionInvPolyShaping::getComponentFunctionThirdDerivative( int componentIndex, double currentTime )
{
    return compositeFunctionComponents_[ componentIndex ]->evaluateThirdDerivative( currentTime );
}


void CompositeZFunctionInvPolyShaping::resetCompositeFunctionCoefficients( Eigen::VectorXd compositeFunctionCoefficients )
{
    // Check whether the size is correct.
    if( compositeFunctionCoefficients.size() == compositeFunctionComponents_.size() )
    {
        compositeFunctionCoefficients_ = compositeFunctionCoefficients;
    }
    else
    {
        std::cerr << "The number of composite function coefficients is incorrect!" << std::endl
             << compositeFunctionCoefficients.size() << " coefficients were put in, however, "
             << compositeFunctionComponents_.size() << " are required. The composite function coefficients are not reset!\n";
    }
}


double CompositeZFunctionInvPolyShaping::evaluateCompositeFunction( const double independentVariable )
{
    double functionValue = 0.0;
    for( unsigned int i = 0; i < compositeFunctionComponents_.size(); i++ )
    {
        functionValue += compositeFunctionCoefficients_[i] * compositeFunctionComponents_[ i ]->evaluateFunction( independentVariable );
    }
    return functionValue;
}


double CompositeZFunctionInvPolyShaping::evaluateCompositeFunctionFirstDerivative( const double independentVariable )
{
    double functionValue = 0.0;
    for( unsigned int i = 0; i < compositeFunctionComponents_.size(); i++ )
    {
        functionValue += compositeFunctionCoefficients_[i] * compositeFunctionComponents_[ i ]->evaluateFirstDerivative( independentVariable );
    }
    return functionValue;
}


double CompositeZFunctionInvPolyShaping::evaluateCompositeFunctionSecondDerivative( const double independentVariable )
{
    double functionValue = 0.0;
    for( unsigned int i = 0; i < compositeFunctionComponents_.size(); i++ )
    {
        functionValue += compositeFunctionCoefficients_[i] * compositeFunctionComponents_[ i ]->evaluateSecondDerivative( independentVariable );
    }
    return functionValue;
}

double CompositeZFunctionInvPolyShaping::evaluateCompositeFunctionThirdDerivative( const double independentVariable )
{
    double functionValue = 0.0;
    for( unsigned int i = 0; i < compositeFunctionComponents_.size(); i++ )
    {
        functionValue += compositeFunctionCoefficients_[i] * compositeFunctionComponents_[ i ]->evaluateThirdDerivative( independentVariable );
    }
    return functionValue;
}

double CompositeZFunctionInvPolyShaping::getComponentFunctionCurrentValue( int componentIndex, double currentTime )
{
    return compositeFunctionComponents_[ componentIndex ]->evaluateFunction( currentTime );
}


double CompositeZFunctionInvPolyShaping::getComponentFunctionFirstDerivative( int componentIndex, double currentTime )
{
    return compositeFunctionComponents_[ componentIndex ]->evaluateFirstDerivative( currentTime );
}

double CompositeZFunctionInvPolyShaping::getComponentFunctionSecondDerivative( int componentIndex, double currentTime )
{
    return compositeFunctionComponents_[ componentIndex ]->evaluateSecondDerivative( currentTime );
}

double CompositeZFunctionInvPolyShaping::getComponentFunctionThirdDerivative( int componentIndex, double currentTime )
{
    return compositeFunctionComponents_[ componentIndex ]->evaluateThirdDerivative( currentTime );
}


} // namespace shape_based_methods
} // namespace tudat
