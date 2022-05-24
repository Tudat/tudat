/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/orbit_determination/estimatable_parameters/sphericalHarmonicSineCoefficients.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/sphericalHarmonicCosineCoefficients.h"

namespace tudat
{

namespace estimatable_parameters
{

//! Function to retrieve the current values of the sine coefficients that are to be estimated.
Eigen::VectorXd SphericalHarmonicsSineCoefficients::getParameterValue( )
{
    Eigen::VectorXd parameterVector = Eigen::VectorXd::Zero( parameterSize_ );

    Eigen::MatrixXd coefficientBlock = getSineCoefficients_( );

    for( unsigned int i = 0; i < blockIndices_.size( ); i++ )
    {
        parameterVector( i ) = coefficientBlock( blockIndices_.at( i ).first, blockIndices_.at( i ).second );
    }
    return parameterVector;

}

//! Function to reset the sine coefficients that are to be estimated.
void SphericalHarmonicsSineCoefficients::setParameterValue( const Eigen::VectorXd parameterValue )
{
    Eigen::MatrixXd coefficients = getSineCoefficients_( );

    for( unsigned int i = 0; i < blockIndices_.size( ); i++ )
    {
        coefficients( blockIndices_.at( i ).first, blockIndices_.at( i ).second ) = parameterValue( i );
    }
    setSineCoefficients_( coefficients );
}

//! Function to get a list of Kaula constraint values for gravity field coefficients for given parameter
Eigen::VectorXd getKaulaConstraintVector(
        const std::shared_ptr< SphericalHarmonicsSineCoefficients > parameter,
        const double constraintMultiplier )
{
    return getKaulaConstraintVector( parameter->getBlockIndices( ), constraintMultiplier );
}

}

}
