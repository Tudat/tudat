/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#include "tudat/math/basic/mathematicalConstants.h"

#include "tudat/astro/electromagnetism/solarSailAcceleration.h"
#include "tudat/astro/electromagnetism/solarSailForce.h"
#include "tudat/astro/electromagnetism/radiationPressureInterface.h"

namespace tudat
{
namespace electromagnetism
{

//! Compute solar sail acceleration with a non-ideal reflective model.
Eigen::Vector3d computeSolarSailAcceleration(
        const double frontEmissivityCoefficient,
        const double backEmissivityCoefficient,
        const double frontLambertianCoefficient,
        const double backLambertianCoefficient,
        const double reflectivityCoefficient,
        const double specularReflectionCoefficient,
        const Eigen::Vector3d& normalisedVectorToSource,
        const Eigen::Vector3d& normalisedVelocityVector,
        const double radiationPressure,
        const double area,
        const double coneAngle,
        const double clockAngle,
        const double mass )
{

    return computeSolarSailForce(
                frontEmissivityCoefficient, backEmissivityCoefficient, frontLambertianCoefficient,
                backLambertianCoefficient, reflectivityCoefficient, specularReflectionCoefficient, normalisedVectorToSource,
                normalisedVelocityVector, radiationPressure, area, coneAngle, clockAngle) / mass;
}

} // namespace electromagnetism
} // namespace tudat
