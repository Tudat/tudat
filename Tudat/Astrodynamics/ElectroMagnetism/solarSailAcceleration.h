/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      Eigen. Structures having Eigen members,
 *          http://eigen.tuxfamily.org/dox/TopicStructHavingEigenMembers.html, last accessed: 5th
 *          March, 2013.
 *
 */

#ifndef TUDAT_SOLAR_SAIL_ACCELERATION_H
#define TUDAT_SOLAR_SAIL_ACCELERATION_H

#include <functional>
#include <memory>
#include <boost/lambda/lambda.hpp>
#include <boost/shared_ptr.hpp>
#include <iostream>

#include <Eigen/Core>

#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h"

#include "Tudat/Astrodynamics/ElectroMagnetism/solarSailForce.h"

namespace tudat
{
namespace electro_magnetism
{

//! Compute radiation pressure acceleration using an ideal and non-ideal model.
/*!
 * Computes radiation pressure acceleration using a cannon-ball model, i.e. assuming force to be in
 * opposite direction of the vector to the source. This function is essentially a wrapper for the
 * function that computes the force.
 * opposite direction of the vector to the source.
 * \param frontEmissivity Parameter determining the emissivity of the front of the sail         [-]
 * \param backEmissivity  Parameter determining the emissivity of the back of the sail          [-]
 * \param frontLambertianCoefficient  Parameter determining the Lambertian coefficient
 *        of the front of the sail                                                              [-]
 * \param backLambertianCoefficient   Parameter determining the Lambertian coefficient
 *        of the back of the sail
 * \param reflectivityCoefficient   Coefficient determining the front reflectivity of
 *        the sail
 * \param specularReflection Coefficient  Coefficient of specular reflection                    [-]
 * \param vectorToSource Vector pointing from target to source. N.B: this must be a unit
 *          vector! To compute the unit vector based on a given position vector, you can
 *          use the .normalize() or .normalized() member functions of an Eigen::Vector3d
 *          object.
 * \param radiationPressure Radiation pressure at targert                                   [N/m^2]                                                  [-]
 * \param area Area on which radiation pressure is assumed to act.                            [m^2]
 * \param coneAngle Sail cone angle                                                           [rad]
 * \param mass Mass of accelerated body.                                                       [kg]
 * \return Acceleration due to radiation pressure.                                          [m/s^2]
 * \sa computeSolarSailAccelerationForce().
 */
Eigen::Vector3d computeSolarSailAcceleration(
        const double frontEmissivityCoefficient,
        const double backEmissivityCoefficient,
        const double frontLambertianCoefficient,
        const double backLambertianCoefficient,
        const double reflectivityCoefficient,
        const double specularReflectionCoefficient,
        const Eigen::Vector3d& vectorToSource,
        const Eigen::Vector3d& velocityUnitVector,
        const double radiationPressure,
        const double area,
        const double coneAngle,
        const double clockAngle,
        const double mass );

//! Solar sail acceleration model class.
/*!
 * Class that can be used to compute the solar sail acceleration using a non-ideal model
 */
class SolarSailAcceleration: public basic_astrodynamics::AccelerationModel3d
{
private:

    //! Typedef for double-returning function.
    typedef std::function< double( ) > DoubleReturningFunction;

    //! Typedef for Eigen::Vector3d-returning function.
    typedef std::function< Eigen::Vector3d( ) > Vector3dReturningFunction;

    //! Typedef for Eigen::Vector6d-returning function.
    typedef std::function< Eigen::VectorXd( ) > VectorXdReturningFunction;

public:

    // Ensure that correctly aligned pointers are generated (Eigen, 2013).
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    //! Constructor taking function pointers for all variables.
    /*!
     * Constructor taking function pointers for all variables.
     * \param sourcePositionFunction Function returning position of radiation source.
     * \param acceleratedBodyPositionFunction Function returning position of body undergoing
     *          acceleration.
     * \param centralBodyVelocityFunction Function returning velocity of central body
     * \param radiationPressureFunction Function returning current radiation pressure.
     * \param coneAngleFunction Function returning current cone angle
     * \param clockAngleFunction Function returning current clock angle
     * \param frontEmissivityCoefficientFunction Function returning current front emissivity coefficient
     * \param backEmissivityCoefficientFunction Function returning current back emissivity coefficient
     * \param frontLambertianCoefficientFunction Function returning current front Lambertian coefficient
     * \param backLambertianCoefficientFunction Function returning current back Lambertian coefficient
     * \param reflectivityCoefficientFunction Function returning current reflectivity coefficient
     * \param specularReflectionCoefficientFunction Function returning current specular reflection coefficient
     * \param areaFunction Function returning current area assumed to undergo radiation pressure.
     * \param massFunction Function returning current mass of body undergoing acceleration.
     */
    SolarSailAcceleration(
            Vector3dReturningFunction sourcePositionFunction,
            Vector3dReturningFunction acceleratedBodyPositionFunction,
            Vector3dReturningFunction acceleratedBodyVelocityFunction,
            Vector3dReturningFunction centralBodyVelocityFunction,
            DoubleReturningFunction radiationPressureFunction,
            DoubleReturningFunction coneAngleFunction,
            DoubleReturningFunction clockAngleFunction,
            DoubleReturningFunction frontEmissivityCoefficientFunction,
            DoubleReturningFunction backEmissivityCoefficientFunction,
            DoubleReturningFunction frontLambertianCoefficientFunction,
            DoubleReturningFunction backLambertianCoefficientFunction,
            DoubleReturningFunction reflectivityCoefficientFunction,
            DoubleReturningFunction specularReflectionCoefficientFunction,
            DoubleReturningFunction areaFunction,
            DoubleReturningFunction massFunction )
        : sourcePositionFunction_( sourcePositionFunction ),
          acceleratedBodyPositionFunction_( acceleratedBodyPositionFunction ),
          acceleratedBodyVelocityFunction_( acceleratedBodyVelocityFunction ),
          centralBodyVelocityFunction_( centralBodyVelocityFunction ),
          radiationPressureFunction_( radiationPressureFunction ),
          coneAngleFunction_( coneAngleFunction ),
          clockAngleFunction_( clockAngleFunction ),
          frontEmissivityCoefficientFunction_( frontEmissivityCoefficientFunction ),
          backEmissivityCoefficientFunction_( backEmissivityCoefficientFunction ),
          frontLambertianCoefficientFunction_( frontLambertianCoefficientFunction ),
          backLambertianCoefficientFunction_( backLambertianCoefficientFunction ),
          reflectivityCoefficientFunction_( reflectivityCoefficientFunction ),
          specularReflectionCoefficientFunction_( specularReflectionCoefficientFunction ),
          areaFunction_( areaFunction ),
          massFunction_( massFunction )
    {
        this->updateMembers( );
    }

    //! Constructor taking functions pointers and constant values for parameters.
    /*!
     * Constructor taking function pointers for position vectors and radiation pressure and
     * constant values for other parameters.
     * \param sourcePositionFunction Function returning position of radiation source.
     * \param acceleratedBodyPositionFunction Function returning position of body undergoing
     *          acceleration.
     * \param radiationPressureFunction Function returning current radiation pressure.
     * \param coneAngleFunction Function returning current cone angle
     * \param clockAngleFunction Function returning current clock angle
     * \param frontEmissivityCoefficient Constant front emissivity coefficient
     * \param backEmissivityCoefficient Constant back emissivity coefficient
     * \param frontLambertianCoefficient Constant front Lambertian coefficient
     * \param backLambertianCoefficient Constant back Lambertian coefficient
     * \param reflectivityCoefficient Constant reflectivity coefficient
     * \param coefficientSpecularReflection Constant specular reflection coefficient
     * \param area Constant area assumed to undergo radiation pressure.
     * \param mass Constant mass of body undergoing acceleration.
     */
    SolarSailAcceleration(
            Vector3dReturningFunction sourcePositionFunction,
            Vector3dReturningFunction acceleratedBodyPositionFunction,
            Vector3dReturningFunction acceleratedBodyVelocityFunction,
            Vector3dReturningFunction centralBodyVelocityFunction,
            DoubleReturningFunction radiationPressureFunction,
            DoubleReturningFunction coneAngleFunction,
            DoubleReturningFunction clockAngleFunction,
            const double frontEmissivityCoefficient,
            const double backEmissivityCoefficient,
            const double frontLambertianCoefficient,
            const double backLambertianCoefficient,
            const double reflectivityCoefficient,
            const double specularReflectionCoefficient,
            const double area,
            const double mass
            )
        : sourcePositionFunction_( sourcePositionFunction ),
          acceleratedBodyPositionFunction_( acceleratedBodyPositionFunction ),
          acceleratedBodyVelocityFunction_( acceleratedBodyVelocityFunction ),
          centralBodyVelocityFunction_( centralBodyVelocityFunction ),
          radiationPressureFunction_( radiationPressureFunction ),
          coneAngleFunction_( coneAngleFunction ),
          clockAngleFunction_( clockAngleFunction ),
          frontEmissivityCoefficientFunction_( [ = ]( ){ return frontEmissivityCoefficient;} ),
          backEmissivityCoefficientFunction_( [ = ]( ){ return backEmissivityCoefficient;} ),
          frontLambertianCoefficientFunction_( [ = ]( ){ return frontLambertianCoefficient;} ),
          backLambertianCoefficientFunction_( [ = ]( ){ return backLambertianCoefficient;} ),
          reflectivityCoefficientFunction_( [ = ]( ){ return reflectivityCoefficient;} ),
          specularReflectionCoefficientFunction_( [ = ]( ){ return specularReflectionCoefficient;} ),
          areaFunction_( [ = ]( ){ return area;} ),
          massFunction_( [ = ]( ){ return mass;} )
    {
        this->updateMembers( );
    }

    //! Get solar sail acceleration.
    /*!
     * Returns the solar sail acceleration. No arguments are passed to this function.
     * Instead, all data required for computation is to be obtained from pointers to functions/
     * classes/structs, etc., which are to be set in a derived class and evaluated by the
     * updateMembers() function below. This function is essentially a wrapper for the free
     * function that computes the solar sail acceleration.
     * \return Solar sail acceleration.
     * \sa computeSolarSailAcceleration().
     */
    Eigen::Vector3d getAcceleration( )
    {
        return computeSolarSailAcceleration(
                    currentFrontEmissivityCoefficient_, currentBackEmissivityCoefficient_,
                    currentFrontLambertianCoefficient_, currentBackLambertianCoefficient_,
                    currentReflectivityCoefficient_, currentSpecularReflectionCoefficient_,
                    currentVectorToSource_, currentVelocityUnitVector_,
                    currentRadiationPressure_, currentArea_,
                    currentConeAngle_, currentClockAngle_, currentMass_ );
    }

    //! Update member variables used by the radiation pressure acceleration model.
    /*!
     * Updates member variables used by the acceleration model. This function evaluates all
     * dependent variables to the 'current' values of these parameters. Only these current values,
     * not the function-pointers are then used by the getAcceleration( ) function.
     * \param currentTime Time at which acceleration model is to be updated.
     */
    void updateMembers( const double currentTime = TUDAT_NAN )
    {
        if( !( this->currentTime_ == currentTime ) )
        {
            currentFrontEmissivityCoefficient_ = frontEmissivityCoefficientFunction_( );
            currentBackEmissivityCoefficient_ = backEmissivityCoefficientFunction_( );
            currentFrontLambertianCoefficient_ = frontLambertianCoefficientFunction_( );
            currentBackLambertianCoefficient_ = backLambertianCoefficientFunction_( );
            currentReflectivityCoefficient_ = reflectivityCoefficientFunction_( );
            currentSpecularReflectionCoefficient_ = specularReflectionCoefficientFunction_( );
            currentVectorToSource_ = ( sourcePositionFunction_( )
                                       - acceleratedBodyPositionFunction_( ) ).normalized( );
            currentVelocityUnitVector_=(acceleratedBodyVelocityFunction_()-centralBodyVelocityFunction_()).normalized();
            currentRadiationPressure_ = radiationPressureFunction_( );
            currentArea_ = areaFunction_( );
            currentConeAngle_ = coneAngleFunction_( );
            currentClockAngle_ = clockAngleFunction_( );
            currentMass_ = massFunction_( );
        }

    }

    //! Function to retrieve the function pointer returning mass of accelerated body.
    /*!
     * Function to retrieve the function pointer returning mass of accelerated body.
     * \return Function pointer returning mass of accelerated body.
     */
    std::function< double( ) > getMassFunction( )
    {
        return massFunction_;
    }

private:

    //! Function pointer returning position of source.
    /*!
     * Function pointer returning position of source (3D vector).
     */
    const Vector3dReturningFunction sourcePositionFunction_;

    //! Function pointer returning position of accelerated body.
    /*!
     * Function pointer returning position of accelerated body (3D vector).
     */
    const Vector3dReturningFunction acceleratedBodyPositionFunction_;

    //! Function pointer returning velocity of accelerated body.
    /*!
     * Function pointer returning velocity of accelerated body (3D vector).
     */
    const Vector3dReturningFunction acceleratedBodyVelocityFunction_;

    //! Function pointer returning velocity of central body.
    /*!
     * Function pointer returning velocity of central body (3D vector).
     */
    const Vector3dReturningFunction centralBodyVelocityFunction_;

    //! Function pointer returning radiation pressure.
    /*!
     * Function pointer returning radiation pressure [N/m^{2}].
     */
    const DoubleReturningFunction radiationPressureFunction_;

    //! Function pointer returning cone angle.
    /*!
     * Function pointer returning cone angle [rad].
     */
    const DoubleReturningFunction coneAngleFunction_;


    //! Function pointer returning clock angle.
    /*!
     * Function pointer returning clock angle [rad].
     */
    const DoubleReturningFunction clockAngleFunction_;


    //! Function pointer returning front emissivity coefficient.
    /*!
     * Function pointer returning front emissivity coefficient [-].
     */
    const DoubleReturningFunction frontEmissivityCoefficientFunction_;


    //! Function pointer returning back emissivity coefficient.
    /*!
     * Function pointer returning back emissivity coefficient [-].
     */
    const DoubleReturningFunction backEmissivityCoefficientFunction_;


    //! Function pointer returning front Lambertian coefficient.
    /*!
     * Function pointer returning front Lambertian coefficient [-].
     */
    const DoubleReturningFunction frontLambertianCoefficientFunction_;


    //! Function pointer returning back Lambertian coefficient.
    /*!
     * Function pointer returning back Lambertian coefficient [-].
     */
    const DoubleReturningFunction backLambertianCoefficientFunction_;


    //! Function pointer returning reflectivity coefficient.
    /*!
     * Function pointer returning reflectivity coefficient [-].
     */
    const DoubleReturningFunction reflectivityCoefficientFunction_;


    //! Function pointer returning specular reflection coefficient.
    /*!
     * Function pointer returning specular reflection coefficient [-].
     */
    const DoubleReturningFunction specularReflectionCoefficientFunction_;


    //! Function pointer returning area on which radiation pressure is acting.
    /*!
     * Function pointer returning area on which radiation pressure is acting [m^{2}].
     */
    const DoubleReturningFunction areaFunction_;

    //! Function pointer returning mass of accelerated body.
    /*!
     * Function pointer returning mass of accelerated body [kg].
     */
    const DoubleReturningFunction massFunction_;


    //! Current front emissivity coefficient.
    /*!
     * Current front emissivity coefficient.
     */
    double currentFrontEmissivityCoefficient_;

    //! Current back emissivity coefficient.
    /*!
     * Current back emissivity coefficient.
     */
    double currentBackEmissivityCoefficient_;

    //! Current front Lambertian coefficient.
    /*!
     * Current front Lambertian coefficient.
     */
    double currentFrontLambertianCoefficient_;

    //! Current back Lambertian coefficient.
    /*!
     * Current back Lambertian coefficient.
     */
    double currentBackLambertianCoefficient_;

    //! Current reflectivity coefficient.
    /*!
     * Current reflectivity coefficient.
     */
    double currentReflectivityCoefficient_;

    //! Current specular reflection coefficient.
    /*!
     * Current specular reflection coefficient.
     */
    double currentSpecularReflectionCoefficient_;

    //! Current vector from accelerated body to source.
    /*!
     * Current vector from accelerated body to source (3D vector).
     */
    Eigen::Vector3d currentVectorToSource_;

    //! Current velocity unit vector of the propagated body.
    /*!
     * Current velocity unit vector of the propagated body (3D vector).
     */
    Eigen::Vector3d currentVelocityUnitVector_;

    //! Current radiation pressure.
    /*!
     * Current radiation pressure [N/m^{2}].
     */
    double currentRadiationPressure_;

    //! Current area on which radiation pressure is acting.
    /*!
     * Current area on which radiation pressure is acting [m^{2}].
     */
    double currentArea_;

    //! Current cone angle.
    /*!
     * Current cone angle [rad].
     */
    double currentConeAngle_;

    //! Current clock angle
    /*!
     * Current clock angle [rad].
     */
    double currentClockAngle_;

    //! Current mass of accelerated body.
    /*!
     * Current mass of accelerated body [kg].
     */
    double currentMass_;

};

//! Typedef for shared-pointer to SolarSailAcceleration.
typedef std::shared_ptr< SolarSailAcceleration > SolarSailAccelerationPointer;

} // namespace electro_magnetism
} // namespace tudat

#endif // TUDAT_SOLAR_SAIL_ACCELERATION_H
