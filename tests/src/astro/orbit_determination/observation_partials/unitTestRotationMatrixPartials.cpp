/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <limits>
#include <boost/test/unit_test.hpp>
#include <boost/make_shared.hpp>
#include <boost/lambda/lambda.hpp>

#include "tudat/interface/spice/spiceInterface.h"
#include "tudat/simulation/estimation_setup/createCartesianStatePartials.h"
#include "tudat/astro/orbit_determination/observation_partials/rotationMatrixPartial.h"
#include "tudat/astro/reference_frames/referenceFrameTransformations.h"
#include "tudat/io/basicInputOutput.h"

namespace tudat
{
namespace unit_tests
{

using namespace tudat::ephemerides;
using namespace tudat::estimatable_parameters;
using namespace tudat::observation_partials;

BOOST_AUTO_TEST_SUITE( test_rotation_matrix_partaisl )

//! Test whether partial derivatives of rotation matrix computed by SimpleRotationalEphemeris works correctly
BOOST_AUTO_TEST_CASE( testSimpleRotationalEphemerisPartials )
{
    // Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Create rotation model
    double nominalRotationRate = 2.0 * mathematical_constants::PI / 86400.0;
    std::shared_ptr< SimpleRotationalEphemeris > rotationalEphemeris =
            std::make_shared< SimpleRotationalEphemeris >(
                spice_interface::computeRotationQuaternionBetweenFrames( "ECLIPJ2000", "IAU_Earth", 1.0E7 ),
                nominalRotationRate, 1.0E7, "ECLIPJ2000", "IAU_Earth" );

    {
        // Create partial object.
        std::shared_ptr< RotationMatrixPartialWrtConstantRotationRate > rotationMatrixPartialObject =
                std::make_shared< RotationMatrixPartialWrtConstantRotationRate >( rotationalEphemeris );

        // Compute partial analytically
        double testTime = 1.0E6;
        Eigen::Matrix3d rotationMatrixPartial =
                rotationMatrixPartialObject->calculatePartialOfRotationMatrixToBaseFrameWrParameter(
                    testTime ).at( 0 );

        Eigen::Matrix3d rotationMatrixDerivativePartial =
                rotationMatrixPartialObject->calculatePartialOfRotationMatrixDerivativeToBaseFrameWrParameter(
                    testTime ).at( 0 );

        // Compute partial numerically.
        double perturbation = 1.0E-12;
        rotationalEphemeris->resetRotationRate( nominalRotationRate + perturbation );
        Eigen::Matrix3d upperturbedRotationMatrix = rotationalEphemeris->getRotationToBaseFrame(
                    testTime).toRotationMatrix( );
        Eigen::Matrix3d upperturbedRotationMatrixDerivative = rotationalEphemeris->getDerivativeOfRotationToBaseFrame(
                    testTime );

        rotationalEphemeris->resetRotationRate( nominalRotationRate - perturbation );
        Eigen::Matrix3d downperturbedRotationMatrix = rotationalEphemeris->getRotationToBaseFrame(
                    testTime).toRotationMatrix( );
        Eigen::Matrix3d downperturbedRotationMatrixDerivative = rotationalEphemeris->getDerivativeOfRotationToBaseFrame(
                    testTime );

        Eigen::Matrix3d numericalRotationMatrixPartial =
                ( upperturbedRotationMatrix - downperturbedRotationMatrix ) / ( 2.0 * perturbation );
        Eigen::Matrix3d numericalRotationMatrixDerivativePartial =
                ( upperturbedRotationMatrixDerivative - downperturbedRotationMatrixDerivative ) / ( 2.0 * perturbation );

        Eigen::Matrix3d matrixDifference = rotationMatrixPartial - numericalRotationMatrixPartial;

        // Compare analytical and numerical result.
        for( unsigned int i = 0; i < 3; i++ )
        {
            for( unsigned int j = 0; j < 3; j++ )
            {
                BOOST_CHECK_SMALL( std::fabs( matrixDifference( i, j ) ), 0.1 );
            }
        }

        matrixDifference = rotationMatrixDerivativePartial - numericalRotationMatrixDerivativePartial;

        // Compare analytical and numerical result.
        for( unsigned int i = 0; i < 3; i++ )
        {
            for( unsigned int j = 0; j < 3; j++ )
            {
                BOOST_CHECK_SMALL( std::fabs( matrixDifference( i, j ) ), 1.0E-5 );
            }
        }
    }

    {
        // Create partial object.
        std::shared_ptr< RotationMatrixPartialWrtPoleOrientation > rotationMatrixPartialObject =
                std::make_shared< RotationMatrixPartialWrtPoleOrientation >( rotationalEphemeris );

        // Compute partial analytically
        double testTime = 1.0E6;
        std::vector< Eigen::Matrix3d > rotationMatrixPartials =
                rotationMatrixPartialObject->calculatePartialOfRotationMatrixToBaseFrameWrParameter(
                    testTime );

        std::vector< Eigen::Matrix3d > rotationMatrixDerivativePartials =
                rotationMatrixPartialObject->calculatePartialOfRotationMatrixDerivativeToBaseFrameWrParameter(
                    testTime );

        Eigen::Vector3d nominalEulerAngles = rotationalEphemeris->getInitialEulerAngles( );
        double perturbedAngle;

        // Compute partial numerically.
        double perturbation = 1.0E-6;
        {


            // Compute partial for right ascension numerically.
            {
                perturbedAngle = nominalEulerAngles( 0 ) + perturbation;
                rotationalEphemeris->resetInitialPoleRightAscensionAndDeclination( perturbedAngle, nominalEulerAngles( 1 ) );
                Eigen::Matrix3d upperturbedRotationMatrix = rotationalEphemeris->getRotationToBaseFrame(
                            testTime).toRotationMatrix( );
                Eigen::Matrix3d upperturbedRotationMatrixDerivative = rotationalEphemeris->getDerivativeOfRotationToBaseFrame(
                            testTime );

                perturbedAngle = nominalEulerAngles( 0 ) - perturbation;
                rotationalEphemeris->resetInitialPoleRightAscensionAndDeclination( perturbedAngle, nominalEulerAngles( 1 ) );
                Eigen::Matrix3d downperturbedRotationMatrix = rotationalEphemeris->getRotationToBaseFrame(
                            testTime).toRotationMatrix( );
                Eigen::Matrix3d downperturbedRotationMatrixDerivative =
                        rotationalEphemeris->getDerivativeOfRotationToBaseFrame(
                            testTime );

                Eigen::Matrix3d numericalRotationMatrixPartial =
                        ( upperturbedRotationMatrix - downperturbedRotationMatrix ) / ( 2.0 * perturbation );
                Eigen::Matrix3d numericalRotationMatrixDerivativePartial =
                        ( upperturbedRotationMatrixDerivative - downperturbedRotationMatrixDerivative ) /
                        ( 2.0 * perturbation );

                Eigen::Matrix3d matrixDifference = rotationMatrixPartials.at( 0 ) - numericalRotationMatrixPartial;

                // Compare analytical and numerical result.
                for( unsigned int i = 0; i < 3; i++ )
                {
                    for( unsigned int j = 0; j < 3; j++ )
                    {
                        BOOST_CHECK_SMALL( std::fabs( matrixDifference( i, j ) ), 1.0E-8 );
                    }
                }

                matrixDifference = rotationMatrixDerivativePartials.at( 0 ) - numericalRotationMatrixDerivativePartial;

                // Compare analytical and numerical result.
                for( unsigned int i = 0; i < 3; i++ )
                {
                    for( unsigned int j = 0; j < 3; j++ )
                    {
                        BOOST_CHECK_SMALL( std::fabs( matrixDifference( i, j ) ), 1.0E-13 );
                    }
                }
            }

            // Compute partial for declination numerically.
            {
                perturbedAngle = nominalEulerAngles( 1 ) + perturbation;
                rotationalEphemeris->resetInitialPoleRightAscensionAndDeclination( nominalEulerAngles( 0 ), perturbedAngle );
                Eigen::Matrix3d upperturbedRotationMatrix = rotationalEphemeris->getRotationToBaseFrame(
                            testTime).toRotationMatrix( );
                Eigen::Matrix3d upperturbedRotationMatrixDerivative = rotationalEphemeris->getDerivativeOfRotationToBaseFrame(
                            testTime );

                perturbedAngle = nominalEulerAngles( 1 ) - perturbation;
                rotationalEphemeris->resetInitialPoleRightAscensionAndDeclination( nominalEulerAngles( 0 ), perturbedAngle );
                Eigen::Matrix3d downperturbedRotationMatrix = rotationalEphemeris->getRotationToBaseFrame(
                            testTime).toRotationMatrix( );
                Eigen::Matrix3d downperturbedRotationMatrixDerivative =
                        rotationalEphemeris->getDerivativeOfRotationToBaseFrame( testTime );

                Eigen::Matrix3d numericalRotationMatrixPartial =
                        ( upperturbedRotationMatrix - downperturbedRotationMatrix ) / ( 2.0 * perturbation );
                Eigen::Matrix3d numericalRotationMatrixDerivativePartial =
                        ( upperturbedRotationMatrixDerivative -
                          downperturbedRotationMatrixDerivative ) / ( 2.0 * perturbation );

                Eigen::Matrix3d matrixDifference = rotationMatrixPartials.at( 1 ) - numericalRotationMatrixPartial;

                // Compare analytical and numerical result.
                for( unsigned int i = 0; i < 3; i++ )
                {
                    for( unsigned int j = 0; j < 3; j++ )
                    {
                        BOOST_CHECK_SMALL( std::fabs( matrixDifference( i, j ) ), 1.0E-8 );
                    }
                }

                matrixDifference = rotationMatrixDerivativePartials.at( 1 ) - numericalRotationMatrixDerivativePartial;

                // Compare analytical and numerical result.
                for( unsigned int i = 0; i < 3; i++ )
                {
                    for( unsigned int j = 0; j < 3; j++ )
                    {
                        BOOST_CHECK_SMALL( std::fabs( matrixDifference( i, j ) ), 1.0E-13 );
                    }
                }
            }
        }
    }
}

//! Test whether partial derivatives of rotation matrix computed by SynchronousRotationalEphemeris works correctly
BOOST_AUTO_TEST_CASE( testSynchronousRotationPartials )
{
    // Define nominal state
    Eigen::Vector6d nominalState =
            tudat::spice_interface::getBodyCartesianStateAtEpoch(
                                       "Mercury", "SSB", "ECLIPJ2000", "None", 1.0E7 );

    // Define nominal state function
    Eigen::Vector6d currentState = nominalState;
    std::function< Eigen::Vector6d( const double, bool ) > relativeStateFunction =
            [ & ]( const double, bool ){ return currentState; };

    // Create rotation model
    std::shared_ptr< tudat::ephemerides::SynchronousRotationalEphemeris > synchronousRotationModel =
            std::make_shared< ephemerides::SynchronousRotationalEphemeris >(
                relativeStateFunction, "SSB", "Mercury_Fixed", "ECLIPJ2000" );

    // Create rotation partial model
    std::shared_ptr< RotationMatrixPartial > rotationMatrixPartialObject =
            std::make_shared< SynchronousRotationMatrixPartialWrtTranslationalState >( synchronousRotationModel );

    // Define test settings
    double testTime = 1.0E7;
    double positionPerturbation = 10000.0;
    double velocityPerturbation = 0.1;

    // Test partials w.r.t. position and velocity components
    std::vector< Eigen::Matrix3d > rotationMatrixPartials =
            rotationMatrixPartialObject->calculatePartialOfRotationMatrixToBaseFrameWrParameter( testTime );
    for( int i = 0; i < 3; i++ )
    {
        currentState = nominalState;
        currentState( i ) += positionPerturbation;
        Eigen::Matrix3d upPerturbedRotationMatrix = synchronousRotationModel->getRotationToBaseFrame(
                    1.0E7 ).toRotationMatrix( );

        currentState = nominalState;
        currentState( i ) -= positionPerturbation;
        Eigen::Matrix3d downPerturbedRotationMatrix = synchronousRotationModel->getRotationToBaseFrame(
                    1.0E7 ).toRotationMatrix( );

        Eigen::Matrix3d relativePartialError =
                ( ( upPerturbedRotationMatrix - downPerturbedRotationMatrix ) / ( 2.0 * positionPerturbation ) -
                rotationMatrixPartials.at( i ) ) / rotationMatrixPartials.at( i ).norm( );

        for( int j = 0; j < 3; j++ )
        {
            for( int k = 0; k < 3; k++ )
            {
                BOOST_CHECK_SMALL( std::fabs( relativePartialError( j, k ) ), 1.0E-8 );
            }
        }

        currentState = nominalState;
        currentState( i + 3 ) += velocityPerturbation;
        upPerturbedRotationMatrix = synchronousRotationModel->getRotationToBaseFrame(
                    1.0E7 ).toRotationMatrix( );

        currentState = nominalState;
        currentState( i + 3 ) -= velocityPerturbation;
        downPerturbedRotationMatrix = synchronousRotationModel->getRotationToBaseFrame(
                    1.0E7 ).toRotationMatrix( );

        relativePartialError =
                        ( ( upPerturbedRotationMatrix - downPerturbedRotationMatrix ) / ( 2.0 * velocityPerturbation ) -
                        rotationMatrixPartials.at( i + 3 ) ) / rotationMatrixPartials.at( i + 3 ).norm( );

        for( int j = 0; j < 3; j++ )
        {
            for( int k = 0; k < 3; k++ )
            {
                BOOST_CHECK_SMALL( std::fabs( relativePartialError( j, k ) ), 1.0E-8 );
            }
        }
    }
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat





