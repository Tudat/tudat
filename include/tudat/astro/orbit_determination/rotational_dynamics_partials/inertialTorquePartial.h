/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_INERTIALTORQUEPARTIALS_H
#define TUDAT_INERTIALTORQUEPARTIALS_H

#include "tudat/astro/gravitation/secondDegreeGravitationalTorque.h"
#include "tudat/astro/orbit_determination/rotational_dynamics_partials/torquePartial.h"
#include "tudat/astro/orbit_determination/rotational_dynamics_partials/inertiaTensorPartial.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/sphericalHarmonicCosineCoefficients.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/sphericalHarmonicSineCoefficients.h"
#include "tudat/math/basic/linearAlgebra.h"

namespace tudat
{

namespace acceleration_partials
{

//! Class to compute the partial derivatives of the inertial torque w.r.t. states and parameters
class InertialTorquePartial: public TorquePartial
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param angularVelocityFunction Function that returns body's body-fixed angular velocity vector
     * \param inertiaTensorFunction Function that returns body's inertia tensor
     * \param inertiaTensorNormalizationFunction Function the inertia tensor normalization factor
     * \param bodyGravitationalParameterFunction Function returning body gravitational parameter
     * \param acceleratedBody Name of body undergoing torque
     */
    InertialTorquePartial(
            const std::function< Eigen::Vector3d( ) > angularVelocityFunction,
            const std::function< Eigen::Matrix3d( ) > inertiaTensorFunction,
            const std::function< double( ) > inertiaTensorNormalizationFunction,
            const std::function< double( ) > bodyGravitationalParameterFunction,
            const std::string acceleratedBody ):
        TorquePartial( acceleratedBody, acceleratedBody, basic_astrodynamics::torque_free ),
        angularVelocityFunction_( angularVelocityFunction ),
        inertiaTensorFunction_( inertiaTensorFunction ),
        getInertiaTensorNormalizationFactor_( inertiaTensorNormalizationFunction ),
        bodyGravitationalParameterFunction_( bodyGravitationalParameterFunction ){ }

    ~InertialTorquePartial( ){ }

    //! Function for setting up and retrieving a function returning a partial w.r.t. a double parameter.
    /*!
     *  Function for setting up and retrieving a function returning a partial w.r.t. a double parameter.
     *  Function returns empty function and zero size indicator for parameters with no dependency for current torque.
     *  \param parameter Parameter w.r.t. which partial is to be taken.
     *  \return Pair of parameter partial function and number of columns in partial (0 for no dependency, 1 otherwise).
     */
    std::pair< std::function< void( Eigen::MatrixXd& ) >, int >
    getParameterPartialFunction( std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter );

    //! Function for setting up and retrieving a function returning a partial w.r.t. a vector parameter.
    /*!
     *  Function for setting up and retrieving a function returning a partial w.r.t. a vector parameter.
     *  Function returns empty function and zero size indicator for parameters with no dependency for current torque.
     *  \param parameter Parameter w.r.t. which partial is to be taken.
     *  \return Pair of parameter partial function and number of columns in partial (0 for no dependency).
     */
    std::pair< std::function< void( Eigen::MatrixXd& ) >, int > getParameterPartialFunction(
            std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter );

    //! Function for calculating the partial of the torque w.r.t. the orientation of the accelerated body.
    /*!
     *  Function for calculating the partial of the torque w.r.t. the orientation of the accelerated body and
     *  adding it to the existing partial block.
     *  \param partialMatrix Block of partial derivatives of torque w.r.t. Orientation of body
     *  undergoing torque where current partial is to be added.
     *  \param addContribution Variable denoting whether to return the partial itself (true) or the negative partial (false).
     *  \param startRow First row in partialMatrix block where the computed partial is to be added.
     *  \param startColumn First column in partialMatrix block where the computed partial is to be added.
     */
    void wrtOrientationOfAcceleratedBody(
            Eigen::Block< Eigen::MatrixXd > partialMatrix,
            const bool addContribution = 1, const int startRow = 0, const int startColumn = 0 )
    { }

    //! Function for calculating the partial of the torque w.r.t. the angular velocity of the accelerated body.
    /*!
     *  Function for calculating the partial of the torque w.r.t. the angular velocity of the accelerated body and
     *  adding it to the existing partial block.
     *  \param partialMatrix Block of partial derivatives of torque w.r.t. angular velocity of body
     *  undergoing torque where current partial is to be added.
     *  \param addContribution Variable denoting whether to return the partial itself (true) or the negative partial (false).
     *  \param startRow First row in partialMatrix block where the computed partial is to be added.
     *  \param startColumn First column in partialMatrix block where the computed partial is to be added.
     */
    void wrtRotationalVelocityOfAcceleratedBody(
            Eigen::Block< Eigen::MatrixXd > partialMatrix,
            const bool addContribution = 1, const int startRow = 0, const int startColumn = 3 )
    {
        if( addContribution )
        {
            partialMatrix.block( startRow, startColumn, 3, 3 ) += currentPartialDerivativeWrtAngularVelocity_;
        }
        else
        {
            partialMatrix.block( startRow, startColumn, 3, 3 ) -= currentPartialDerivativeWrtAngularVelocity_;
        }
    }

    //! Update partial model to current time
    /*!
     * Update partial model to current time, retrieving information from consituent functions
     * \param currentTime Time to which partial is to be updated
     */
    void update( const double currentTime = TUDAT_NAN )
    {
        if( !( currentTime_ == currentTime ) )
        {
            currentAngularVelocityVector_ = angularVelocityFunction_( );
            currentAngularVelocityCrossProductMatrix_ = linear_algebra::getCrossProductMatrix(
                        currentAngularVelocityVector_ );

            currentInertiaTensorNormalizationFactor_ = getInertiaTensorNormalizationFactor_( );
            currentGravitationalParameter_ = bodyGravitationalParameterFunction_( );

            currentInertiaTensor_ = inertiaTensorFunction_( );
            currentInverseInertiaTensor_ = currentInertiaTensor_.inverse( );

            currentPartialDerivativeWrtAngularVelocity_ =
                    -linear_algebra::getCrossProductMatrix( currentAngularVelocityVector_ ) * currentInertiaTensor_ +
                    linear_algebra::getCrossProductMatrix( currentInertiaTensor_ * currentAngularVelocityVector_ );
        }
    }

protected:

    //! Function to compute partial of torque w.r.t. mean moment of inertia
    /*!
     * Function to compute partial of torque w.r.t. mean moment of inertia
     * \param momentOfInertiaPartial Computed partial of torque w.r.t. mean moment of inertia (returned by reference)
     */
    void wrtMeanMomentOfInertia(
            Eigen::MatrixXd& momentOfInertiaPartial );

    //! Function to compute partial of torque w.r.t. gravitational parameter
    /*!
     * Function to compute partial of torque w.r.t. gravitational parameter
     * \param gravitationalParameterPartial Computed partial of torque w.r.t. gravitational parameter (returned by reference)
     */
    void wrtGravitationalParameter(
            Eigen::MatrixXd& gravitationalParameterPartial );

    //! Function to compute partial of torque w.r.t. spherical harmonic cosine coefficients
    /*!
     * Function to compute partial of torque w.r.t. spherical harmonic cosine coefficients
     * \param sphericalHarmonicCoefficientPartial Computed partial of torque w.r.t. spherical harmonic cosine coefficients
     * (returned by reference)
     * \param c20Index for degree=2,order=0 coefficient
     * \param c21Index for degree=2,order=1 coefficient
     * \param c22Index for degree=2,order=2 coefficient
     */
    void wrtCosineSphericalHarmonicCoefficientsOfCentralBody(
            Eigen::MatrixXd& sphericalHarmonicCoefficientPartial,
            const int c20Index, const int c21Index, const int c22Index );

    //! Function to compute partial of torque w.r.t. spherical harmonic sine coefficients
    /*!
     * Function to compute partial of torque w.r.t. spherical harmonic sine coefficients
     * \param sphericalHarmonicCoefficientPartial Computed partial of torque w.r.t. spherical harmonic sine coefficients
     * (returned by refe
     * \param s21Index for degree=2,order=1 coefficient
     * \param s22Index for degree=2,order=2 coefficient)
     */
    void wrtSineSphericalHarmonicCoefficientsOfCentralBody(
            Eigen::MatrixXd& sphericalHarmonicCoefficientPartial,
            const int s21Index, const int s22Index );

    //! Function returning body angular velocity vector in body fixed frame.
    std::function< Eigen::Vector3d( ) > angularVelocityFunction_;

    //! Function returning body inertia tensor
    std::function< Eigen::Matrix3d( ) > inertiaTensorFunction_;

    //! Function the inertia tensor normalization factor
    std::function< double( ) > getInertiaTensorNormalizationFactor_;

    //!  Function returning body gravitational parameter
    std::function< double( ) > bodyGravitationalParameterFunction_;

    //! Current angular velocity vector
    Eigen::Vector3d currentAngularVelocityVector_;

    //! Current cross-product matrix from angular velocity vector
    Eigen::Matrix3d currentAngularVelocityCrossProductMatrix_;

    //! Current inertia tensor
    Eigen::Matrix3d currentInertiaTensor_;

    //! Current inverse inertia tensor
    Eigen::Matrix3d currentInverseInertiaTensor_;

    //! Current partial derivative w.r.t. angular velocity vector
    Eigen::Matrix3d currentPartialDerivativeWrtAngularVelocity_;

    //! Current inertia tensor normalization
    double currentInertiaTensorNormalizationFactor_;

    //! Current gravitational parameter
    double currentGravitationalParameter_;

};

} // namespace acceleration_partials

} // namespace tudat

#endif // TUDAT_INERTIALTORQUEPARTIALS_H
