/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <iostream>

#include "tudat/astro/gravitation/secondDegreeGravitationalTorque.h"

namespace tudat
{

namespace gravitation
{

//! Function to calculate the gravitational torque exerted by a point mass on a body with degree two gravity field
Eigen::Vector3d calculateSecondDegreeGravitationalTorque(
        const Eigen::Vector3d& relativePositionOfBodySubjectToTorque,
        const double premultipler,
        const Eigen::Vector3d& inertiaTensorTimesRelativePositionOfBody )
{
    return premultipler *
            relativePositionOfBodySubjectToTorque.cross( inertiaTensorTimesRelativePositionOfBody );
}

//! Function to calculate the gravitational torque exerted by a point mass on a body with degree two gravity field
Eigen::Vector3d calculateSecondDegreeGravitationalTorque(
        const Eigen::Vector3d& relativePositionOfBodySubjectToTorque,
        const double gravitationalParameterOfAttractingBody,
        const Eigen::Matrix3d& inertiaTensorOfRotatingBody )
{
    double distanceBetweenBodies = relativePositionOfBodySubjectToTorque.norm( ) ;
    Eigen::Vector3d multipliedRelativePosition =  inertiaTensorOfRotatingBody * relativePositionOfBodySubjectToTorque;
    return calculateSecondDegreeGravitationalTorque( relativePositionOfBodySubjectToTorque,
             3.0 * gravitationalParameterOfAttractingBody / std::pow( distanceBetweenBodies, 5.0 ),
             multipliedRelativePosition );

}

} // namespace gravitation


} // namespace tudat
