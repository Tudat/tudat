/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <vector>
#include "tudat/math/basic/mathematicalConstants.h"
#include "tudat/math/basic/coordinateConversions.h"
#include "tudat/astro/basic_astro/stateVectorIndices.h"
#include "tudat/astro/basic_astro/sphericalStateConversions.h"
#include "tudat/astro/reference_frames/referenceFrameTransformations.h"

namespace tudat
{

namespace orbital_element_conversions
{

//! Calculate current heading angle.
double calculateHeadingAngle( const Eigen::Vector3d& velocityInVerticalFrame )
{
    return std::atan2( velocityInVerticalFrame( 1 ), velocityInVerticalFrame( 0 ) );
}

//! Calculatre current flight path angle.
double calculateFlightPathAngle( const Eigen::Vector3d& velocityInVerticalFrame )
{
    return -std::asin( velocityInVerticalFrame( 2 ) / velocityInVerticalFrame.norm( ) );
}

//! Function to convert a Cartesian to a spherical orbital state.
Eigen::Vector6d convertCartesianToSphericalOrbitalState(
        const Eigen::Vector6d& bodyFixedCartesianState )
{
    Eigen::Vector6d sphericalOrbitalState;

    // Compute and set spherical position
    Eigen::Vector3d sphericalPosition = coordinate_conversions::convertCartesianToSpherical< double >(
                bodyFixedCartesianState.segment( 0, 3 ) );
    sphericalOrbitalState( radiusIndex ) = sphericalPosition( 0 );
    sphericalOrbitalState( latitudeIndex ) = mathematical_constants::PI / 2.0 - sphericalPosition( 1 );
    sphericalOrbitalState( longitudeIndex ) = sphericalPosition( 2 );

    // Convert velocity to vertical frame.
    Eigen::Vector3d verticalFrameVelocity =
            reference_frames::getRotatingPlanetocentricToLocalVerticalFrameTransformationQuaternion(
                sphericalOrbitalState( longitudeIndex ), sphericalOrbitalState( latitudeIndex ) ) *
            bodyFixedCartesianState.segment( 3, 3 );

    // Set velocity in spherical orbital state.
    sphericalOrbitalState( speedIndex ) = verticalFrameVelocity.norm( );
    sphericalOrbitalState( flightPathIndex ) = calculateFlightPathAngle( verticalFrameVelocity );
    sphericalOrbitalState( headingAngleIndex ) = calculateHeadingAngle( verticalFrameVelocity );

    return sphericalOrbitalState;

}

} // namespace tudat

} // namespace orbital_element_conversions

