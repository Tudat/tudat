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

#include "Tudat/Astrodynamics/ElectroMagnetism/solarSailForce.h"
#include <iostream>

namespace tudat
{
namespace electro_magnetism
{

//! Compute solar sail force using a non-ideal reflective model.
Eigen::Vector3d computeSolarSailForce(const double frontEmissivityCoefficient, const double backEmissivityCoefficient, const double frontLambertianCoefficient,
        const double backLambertianCoefficient , const double reflectivityCoefficient, const double specularReflectionCoefficient,
        const Eigen::Vector3d& vectorToSource,
        const Eigen::Vector3d& velocityUnitVector,
        const double radiationPressure,
        const double area, const double coneAngle, const double clockAngle)
{

    //see literature study for the definitions of these equations:
    //I define a constant to improve the readability of the formulas
    double A= frontLambertianCoefficient*reflectivityCoefficient*(1-specularReflectionCoefficient) + (1-reflectivityCoefficient)*
              (frontEmissivityCoefficient*frontLambertianCoefficient-backEmissivityCoefficient*backLambertianCoefficient) /
              (frontEmissivityCoefficient+backEmissivityCoefficient);


    double cosConeAngle=cos(coneAngle);

    //double radPressure=4.56e-6; //[N/m^2]
    double Force_magnitude = radiationPressure*area*cosConeAngle* std::sqrt( std::pow( cosConeAngle*(1+reflectivityCoefficient*specularReflectionCoefficient) + A, 2 )
                           + std::pow((1-specularReflectionCoefficient*reflectivityCoefficient)*sin(coneAngle) , 2 ) );

    //I now have to compute the phi and theta angles, useful for retrieving the force direction:
    double phi = atan2(( (1-specularReflectionCoefficient*reflectivityCoefficient)*cosConeAngle*sin(coneAngle) ),
                       ( (1+reflectivityCoefficient*specularReflectionCoefficient)*std::pow( cosConeAngle, 2 ) +
                       frontLambertianCoefficient*(1-specularReflectionCoefficient)*reflectivityCoefficient*cosConeAngle +
                       (1-reflectivityCoefficient)*(frontEmissivityCoefficient*frontLambertianCoefficient -
                       backEmissivityCoefficient*backLambertianCoefficient)/( frontEmissivityCoefficient+backEmissivityCoefficient )*cosConeAngle ));
    double theta = coneAngle-phi;


    Eigen::Vector3d ForceDirectionLocalFrame = { cos(theta), sin(theta)*sin(clockAngle), sin(theta)*cos(clockAngle) }; //this is the unit vector of the force defined in the (r,theta,k)

    // For efficiency, I define the vectors' components
    //I have to put a minus sign in front of the vectorToSource, because this is the vector which goes from the body to the source of
    //radiation, whereas I need the one from the source of radiation to the body
    double sX=-vectorToSource[0];
    double sY=-vectorToSource[1];
    double sZ=-vectorToSource[2];
    double vX=velocityUnitVector[0];
    double vY=velocityUnitVector[1];
    double vZ=velocityUnitVector[2];


//    std::cout << "pos x,y,z" << std::endl;
//    std::cout << -vectorToSource << std::endl;
//    std::cout << "vel x,y,z" << std::endl;
//    std::cout << velocityUnitVector << std::endl;


    Eigen::Matrix3d rotationMatrixFromLocalToInertial;
    rotationMatrixFromLocalToInertial(0,0) = sX;
    rotationMatrixFromLocalToInertial(0,1) = vX*std::pow(sZ,2)-vZ*sX*sZ-vY*sX*sY+vX*std::pow(sY,2);
    rotationMatrixFromLocalToInertial(0,2) = vZ*sY-vY*sZ;
    rotationMatrixFromLocalToInertial(1,0) = sY;
    rotationMatrixFromLocalToInertial(1,1) = vY*std::pow(sZ,2)-vZ*sY*sZ+vY*std::pow(sX,2)-vX*sY*sX;
    rotationMatrixFromLocalToInertial(1,2) = vX*sZ-vZ*sX;
    rotationMatrixFromLocalToInertial(2,0) = sZ;
    rotationMatrixFromLocalToInertial(2,1) = vZ*std::pow(sY,2)-vY*sY*sZ-vX*sZ*sX+vZ*std::pow(sX,2);
    rotationMatrixFromLocalToInertial(2,2) = vY*sX-vX*sY;
    Eigen::Vector3d ForceDirection=rotationMatrixFromLocalToInertial*ForceDirectionLocalFrame;

    return Force_magnitude*ForceDirection;
}

} // namespace electro_magnetism
} // namespace tudat
