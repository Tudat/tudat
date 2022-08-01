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

#ifndef TUDAT_CENTRAL_J2_GRAVITY_MODEL_H
#define TUDAT_CENTRAL_J2_GRAVITY_MODEL_H

#include <boost/lambda/lambda.hpp>
#include <memory>

#include "tudat/astro/basic_astro/accelerationModel.h"
#include "tudat/astro/gravitation/centralGravityModel.h"
#include "tudat/astro/gravitation/sphericalHarmonicsGravityModelBase.h"

namespace tudat
{
namespace gravitation
{

//! Compute gravitational acceleration due to J2.
/*!
 * Computes gravitational acceleration experienced by body1, due to its interaction with another
 * body (body2) with an irregular gravity field. The body exerting the acceleration has an
 * irregular gravity field. The acceleration due to the J2-coefficient is given by (Melman, 2012):
 * \f{eqnarray*}{
 *      {a}_{gravity,x} &=& -\mu_{2}*\frac{x_{1}-x_{2}}{r^{3}} * \frac{3}{2}*J_{2}
 *                          * \left(\frac{R}{r}\right)^{2} * (1-5*\hat{z}^2)  \\
 *      {a}_{gravity,y} &=& -\mu_{2}*\frac{y_{1}-y_{2}}{r^{3}} * \frac{3}{2}*J_{2}
 *                          * \left(\frac{R}{r}\right)^{2} * (1-5*\hat{z}^2)  \\
 *      {a}_{gravity,z} &=& \frac{-\mu_{2}}{r^{2}} * \frac{3}{2}*J_{2}
 *                          * \left(\frac{R}{r}\right)^{2} * (3-5*\hat{z}^2) * \hat{z}  \\
 * \f}
 * where \f$\mu_{2}\f$ is the gravitational parameter of the body exerting acceleration,
 * \f$\hat{z} = \frac{z}{r}\f$, \f$x\f$, \f$y\f$ and \f$z\f$ are the Cartesian position components,
 * \f$r\f$ is the radial position, and \f$J_{2}\f$ is the second zonal coefficient of the gravity
 * field of body2. The positions and accelerations are given with respect to an inertial
 * (barycentric) reference frame.
 * \param positionOfBodySubjectToAcceleration Position vector of body subject to J2-acceleration
 *          (body1) [m].
 * \param gravitationalParameterOfBodyExertingAcceleration Gravitational parameter of body exerting
 *          acceleration (body2) [m^3 s^-2].
 * \param equatorialRadiusOfBodyExertingAcceleration Equatorial radius of body exerting
 *          acceleration (body2), in formulation of spherical harmonics expansion [m].
 * \param j2CoefficientOfGravityField J2-coefficient, describing irregularity of the gravity field
 *          of body2 [-].
 * \param positionOfBodyExertingAcceleration Position vector of body exerting acceleration
 *          (body2) [m].
 * \return Gravitational acceleration exerted on body1 due to J2-effect [m s^-2].
 */
Eigen::Vector3d computeGravitationalAccelerationDueToJ2(
        const Eigen::Vector3d& positionOfBodySubjectToAcceleration,
        const double gravitationalParameterOfBodyExertingAcceleration,
        const double equatorialRadiusOfBodyExertingAcceleration,
        const double j2CoefficientOfGravityField,
        const Eigen::Vector3d& positionOfBodyExertingAcceleration );

//! Central + J2 gravitational acceleration model class.
/*!
 * This class implements a gravitational acceleration model that includes the central and J2
 * (unnormalized coefficient of general spherical harmonics expansion) terms.
 */
class CentralJ2GravitationalAccelerationModel
        : public basic_astrodynamics::AccelerationModel< Eigen::Vector3d >,
        public SphericalHarmonicsGravitationalAccelerationModelBase< Eigen::Vector3d >
{
private:

    //! Typedef for base class.
    typedef SphericalHarmonicsGravitationalAccelerationModelBase< Eigen::Vector3d > Base;

public:

    //! Constructor taking position-functions for bodies, and constant parameters of spherical
    //! harmonics expansion.
    /*!
     * Constructor taking a pointer to a function returning the position of the body subject to
     * gravitational acceleration, constant gravitational parameter, J2-coefficient, and equatorial
     * radius of the body exerting the acceleration, and a pointer to a function returning the
     * position of the body exerting the gravitational acceleration (typically the central body).
     * This constructor uses the Boost::lambda library to create a function on-the-fly that returns
     * the constant gravitational parameter, J2-coefficient, and equatorial radius provided. The
     * constructor also updates all the internal members. The position of the body exerting the
     * gravitational acceleration is an optional parameter; the default position is the origin.
     * \param positionOfBodySubjectToAccelerationFunction Pointer to function returning position of
     *          body subject to gravitational acceleration.
     * \param aGravitationalParameter A (constant) gravitational parameter [m^2 s^-3].
     * \param anEquatorialRadius A (constant) equatorial radius [m].
     * \param aJ2GravityCoefficient A (constant) J2-coefficient.
     * \param positionOfBodyExertingAccelerationFunction Pointer to function returning position of
     *          body exerting gravitational acceleration (default = (0,0,0)).
     */
    CentralJ2GravitationalAccelerationModel(
            const Base::StateFunction positionOfBodySubjectToAccelerationFunction,
            const double aGravitationalParameter,
            const double anEquatorialRadius,
            const double aJ2GravityCoefficient,
            const Base::StateFunction positionOfBodyExertingAccelerationFunction =
            [ ]( ){ return Eigen::Vector3d::Zero( ); } )
        : Base( positionOfBodySubjectToAccelerationFunction,
                aGravitationalParameter,
                positionOfBodyExertingAccelerationFunction, false ),
          equatorialRadius( anEquatorialRadius ),
          j2GravityCoefficient( aJ2GravityCoefficient )
    {
        this->updateMembers( );
    }

    //! Get gravitational acceleration.
    /*!
     * Returns the gravitational acceleration computed using the input parameters provided to the
     * class. This function serves as a wrapper for the computeGravitationalAcceleration() and
     * computeGravitationalAccelerationDueToJ2() functions.
     * \return Computed gravitational acceleration vector.
     */
    Eigen::Vector3d getAcceleration( );

    //! Update members.
    /*!
     * Updates class members relevant for computing the central gravitational acceleration. In this
     * case the function simply updates the members in the base class.
     * \sa SphericalHarmonicsGravitationalAccelerationModelBase.
     * \param currentTime Time at which acceleration model is to be updated.
     */
    void updateMembers( const double currentTime = TUDAT_NAN )
    {
        if( !( this->currentTime_ == currentTime ) )
        {
            this->updateBaseMembers( );
        }
    }

protected:

    //! Equatorial radius [m].
    /*!
     * Equatorial radius of unnormalized spherical harmonics gravity field
     * representation [m].
     */
    const double equatorialRadius;

    //! J2 gravity coefficient.
    /*!
     * J2 coefficient of unnormalized spherical harmonics gravity field representation.
     */
    const double j2GravityCoefficient;

private:
};

//! Typedef for shared-pointer to CentralJ2GravitationalAccelerationModel.
typedef std::shared_ptr< CentralJ2GravitationalAccelerationModel >
CentralJ2GravitationalAccelerationModelPointer;

} // namespace gravitation
} // namespace tudat

#endif // TUDAT_CENTRAL_J2_GRAVITY_MODEL_H
