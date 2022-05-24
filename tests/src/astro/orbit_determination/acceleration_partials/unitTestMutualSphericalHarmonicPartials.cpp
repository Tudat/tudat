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
#include <string>
#include "tudat/basics/testMacros.h"
#include "tudat/astro/basic_astro/orbitalElementConversions.h"
#include "tudat/astro/basic_astro/unitConversions.h"

#include <boost/test/unit_test.hpp>
#include <boost/make_shared.hpp>
#include <boost/lambda/lambda.hpp>

#include "tudat/astro/aerodynamics/exponentialAtmosphere.h"
#include "tudat/astro/gravitation/centralGravityModel.h"
#include "tudat/astro/gravitation/sphericalHarmonicsGravityField.h"
#include "tudat/interface/spice/spiceInterface.h"
#include "tudat/io/basicInputOutput.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/gravitationalParameter.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/constantRotationRate.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/constantRotationalOrientation.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/sphericalHarmonicCosineCoefficients.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/sphericalHarmonicSineCoefficients.h"
#include "tudat/astro/orbit_determination/acceleration_partials/mutualSphericalHarmonicGravityPartial.h"
#include "tudat/astro/orbit_determination/acceleration_partials/numericalAccelerationPartial.h"

#include "tudat/astro/gravitation/mutualSphericalHarmonicGravityModel.h"
#include "tudat/simulation/estimation_setup/createEstimatableParameters.h"
#include "tudat/simulation/simulation.h"


using namespace tudat::ephemerides;
using namespace tudat::simulation_setup;
using namespace tudat::orbit_determination;
using namespace tudat::spice_interface;
using namespace tudat::orbit_determination;
using namespace tudat::estimatable_parameters;
using namespace tudat::gravitation;
using namespace tudat::acceleration_partials;

namespace tudat
{
namespace unit_tests
{



BOOST_AUTO_TEST_SUITE( test_mutual_sh_acceleration_partials )

//! Retrieve Phobos gravity field
std::shared_ptr< GravityFieldSettings > getPhobosGravityFieldSettings( )
{
    Eigen::MatrixXd phobosCosineCoefficients = Eigen::MatrixXd::Zero( 20, 20 );

    phobosCosineCoefficients( 0, 0 ) = 1.0;
    phobosCosineCoefficients( 2, 0 ) = -0.072 / basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 0 );
    phobosCosineCoefficients( 2, 2 ) = -0.048 /basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 2 );


    return std::make_shared< SphericalHarmonicsGravityFieldSettings >(
                7.087546066894452E05, 12.0E3, phobosCosineCoefficients,  Eigen::MatrixXd::Zero( 20, 20 ), "IAU_Phobos" );
}

//! Retrieve Moon gravity field and use as Mars gravity field (makes no difference for test purposes)
std::shared_ptr< GravityFieldSettings > getMarsGravityFieldSettings( )
{
    std::pair< Eigen::MatrixXd, Eigen::MatrixXd > coefficients;
    std::pair< double, double > moonGravityFieldSettings = readGravityFieldFile(
                paths::getGravityModelsPath( ) + "/Moon/lpe200.txt", 50, 50, coefficients, 0, 1 );

    return std::make_shared< SphericalHarmonicsGravityFieldSettings >
            ( moonGravityFieldSettings.first, moonGravityFieldSettings.second,
              coefficients.first, coefficients.second, "IAU_Moon" );
}

//! Test whether partials of mutual spherical harmonic acceleration are computed correctly
BOOST_AUTO_TEST_CASE( testMutualSphericalHarmonicGravityPartials )
{
    // Load spice kernel.
    spice_interface::loadStandardSpiceKernels( { paths::getSpiceKernelPath( ) + "/de430_mar097_small.bsp" } );

    std::vector< std::string > bodyList;
    bodyList.push_back( "Sun" );
    bodyList.push_back( "Mars" );

    // Run test for inertial and Mars-centered acceleration.
    for( int testCase = 0; testCase < 2; testCase++ )
    {

        double initialTime = 1.0E6 - 1.0E5;

        // Retrieve body settings
        BodyListSettings bodySettings = getDefaultBodySettings(
                    bodyList, 1.0E6 - 1.0E5, 1.0E6 + 1.0E5 );
        bodySettings.addSettings( "Phobos" );
        bodySettings.at( "Phobos" )->ephemerisSettings = getDefaultEphemerisSettings(
                    "Phobos" );

        // Update rotation models
        bodySettings.at( "Mars" )->rotationModelSettings = std::make_shared< SimpleRotationModelSettings >(
                    "ECLIPJ2000", "IAU_Mars",
                    spice_interface::computeRotationQuaternionBetweenFrames(
                        "ECLIPJ2000", "IAU_Mars", initialTime ),
                    initialTime, 2.0 * mathematical_constants::PI /
                    ( physical_constants::JULIAN_DAY + 40.0 * 60.0 ) );
        bodySettings.at( "Phobos" )->rotationModelSettings = std::make_shared< SimpleRotationModelSettings >(
                    "ECLIPJ2000", "IAU_Phobos",
                    spice_interface::computeRotationQuaternionBetweenFrames(
                        "ECLIPJ2000", "IAU_Phobos", initialTime ),
                    initialTime, 2.0 * mathematical_constants::PI /
                    ( physical_constants::JULIAN_DAY / 4.0 ) );

        // Update gravity field settings.
        bodySettings.at( "Mars" )->gravityFieldSettings = getMarsGravityFieldSettings( );
        bodySettings.at( "Phobos" )->gravityFieldSettings = getPhobosGravityFieldSettings( );

        // Create body objects
        SystemOfBodies bodies = createSystemOfBodies( bodySettings );
        

        // Create links to set and get state functions of bodies.
        std::shared_ptr< Body > mars = bodies.at( "Mars" );
        std::function< void( Eigen::Vector6d ) > marsStateSetFunction = std::bind( &Body::setState, mars, std::placeholders::_1  );
        std::function< Eigen::Vector6d( ) > marsStateGetFunction = std::bind( &Body::getState, mars );
        mars->setStateFromEphemeris( 1.0E6 );

        std::shared_ptr< Body > phobos = std::dynamic_pointer_cast< Body >( bodies.at( "Phobos" ) );
        std::function< void( Eigen::Vector6d ) > phobosStateSetFunction = std::bind( &Body::setState, phobos, std::placeholders::_1  );
        std::function< Eigen::Vector6d( ) > phobosStateGetFunction = std::bind( &Body::getState, phobos );
        phobos->setStateFromEphemeris( 1.0E6 );

        // Calculate and set mars orientation at epoch.
        mars->setCurrentRotationToLocalFrameFromEphemeris( 1.0E6 );
        mars->setStateFromEphemeris( 1.0E6 );
        phobos->setCurrentRotationToLocalFrameFromEphemeris( 1.0E6 );
        phobos->setStateFromEphemeris( 1.0E6 );

        // Define central body
        std::string centralBody = "";
        std::shared_ptr< Body > centralBodyObject;
        if( testCase == 1 )
        {
            centralBody = "Mars";
            centralBodyObject = mars;
        }

        // Create acceleration model.
        std::shared_ptr< AccelerationSettings > accelerationSettings =
                std::make_shared< MutualSphericalHarmonicAccelerationSettings >( 10, 10, 6, 6, 0, 0 );
        std::shared_ptr< MutualSphericalHarmonicsGravitationalAccelerationModel > accelerationModel =
                std::dynamic_pointer_cast< MutualSphericalHarmonicsGravitationalAccelerationModel >(
                    createAccelerationModel( phobos, mars, accelerationSettings, "Phobos", "Mars",
                                             centralBodyObject, centralBody ) );

        // Retrieve gravity fields
        std::shared_ptr< SphericalHarmonicsGravityField > marsGravityField =
                std::dynamic_pointer_cast< SphericalHarmonicsGravityField  >(
                    bodies.at( "Mars" )->getGravityFieldModel( ) );

        std::shared_ptr< SphericalHarmonicsGravityField > phobosGravityField =
                std::dynamic_pointer_cast< SphericalHarmonicsGravityField  >(
                    bodies.at( "Phobos" )->getGravityFieldModel( ) );

        // Create gravitational parameter estimation settings.
        std::shared_ptr< EstimatableParameter< double > > marsGravitationalParameterObject = std::make_shared<
                GravitationalParameter >( marsGravityField, "Mars" );
        std::shared_ptr< EstimatableParameter< double > > phobosGravitationalParameterObject = std::make_shared<
                GravitationalParameter >( phobosGravityField, "Phobos" );

        // Create rotation parameter estimation settings.
        std::shared_ptr< EstimatableParameter< double > > marsRotationRate = std::make_shared<
                RotationRate >( std::dynamic_pointer_cast< SimpleRotationalEphemeris >(
                                    mars->getRotationalEphemeris( ) ), "Mars" );
        std::shared_ptr< EstimatableParameter< Eigen::VectorXd > > marsPolePosition = std::make_shared<
                ConstantRotationalOrientation >( std::dynamic_pointer_cast< SimpleRotationalEphemeris >(
                                                     mars->getRotationalEphemeris( ) ), "Mars" );
        std::shared_ptr< EstimatableParameter< double > > phobosRotationRate = std::make_shared<
                RotationRate >( std::dynamic_pointer_cast< SimpleRotationalEphemeris >(
                                    phobos->getRotationalEphemeris( ) ), "Phobos" );
        std::shared_ptr< EstimatableParameter< Eigen::VectorXd > > phobosPolePosition = std::make_shared<
                ConstantRotationalOrientation >( std::dynamic_pointer_cast< SimpleRotationalEphemeris >(
                                                     phobos->getRotationalEphemeris( ) ), "Phobos" );

        // Create Mars gravity field coefficient estimation settings.
        std::function< Eigen::MatrixXd( ) > getSineCoefficientsFunction =
                std::bind( &SphericalHarmonicsGravityField::getSineCoefficients, marsGravityField );
        std::function< void( Eigen::MatrixXd ) > setSineCoefficientsFunction =
                std::bind( &SphericalHarmonicsGravityField::setSineCoefficients, marsGravityField, std::placeholders::_1 );
        std::function< Eigen::MatrixXd( ) > getCosineCoefficientsFunction =
                std::bind( &SphericalHarmonicsGravityField::getCosineCoefficients, marsGravityField );
        std::function< void( Eigen::MatrixXd ) > setCosineCoefficientsFunction =
                std::bind( &SphericalHarmonicsGravityField::setCosineCoefficients, marsGravityField, std::placeholders::_1 );

        std::shared_ptr< EstimatableParameter< Eigen::VectorXd > > marsSineCoefficients =
                std::make_shared< SphericalHarmonicsSineCoefficients >(
                    getSineCoefficientsFunction, setSineCoefficientsFunction,
                    getSphericalHarmonicBlockIndices( 3, 1, 5, 5 ), "Mars" );
        std::shared_ptr< EstimatableParameter< Eigen::VectorXd > > marsCosineCoefficients =
                std::make_shared< SphericalHarmonicsCosineCoefficients >(
                    getCosineCoefficientsFunction, setCosineCoefficientsFunction,
                    getSphericalHarmonicBlockIndices( 3, 0, 5, 5 ), "Mars" );

        // Create Phobos gravity field coefficient estimation settings.
        getSineCoefficientsFunction =
                std::bind( &SphericalHarmonicsGravityField::getSineCoefficients, phobosGravityField );
        setSineCoefficientsFunction =
                std::bind( &SphericalHarmonicsGravityField::setSineCoefficients, phobosGravityField, std::placeholders::_1 );
        getCosineCoefficientsFunction =
                std::bind( &SphericalHarmonicsGravityField::getCosineCoefficients, phobosGravityField );
        setCosineCoefficientsFunction =
                std::bind( &SphericalHarmonicsGravityField::setCosineCoefficients, phobosGravityField, std::placeholders::_1 );

        std::shared_ptr< EstimatableParameter< Eigen::VectorXd > > phobosSineCoefficients =
                std::make_shared< SphericalHarmonicsSineCoefficients >(
                    getSineCoefficientsFunction, setSineCoefficientsFunction,
                    getSphericalHarmonicBlockIndices( 2, 1, 2, 2 ), "Phobos" );
        std::shared_ptr< EstimatableParameter< Eigen::VectorXd > > phobosCosineCoefficients =
                std::make_shared< SphericalHarmonicsCosineCoefficients >(
                    getCosineCoefficientsFunction, setCosineCoefficientsFunction,
                    getSphericalHarmonicBlockIndices( 2, 0, 2, 2 ), "Phobos" );

        // Create parameter objects
        std::vector< std::shared_ptr< EstimatableParameter< double > > > doubleParameters;
        doubleParameters.push_back( marsGravitationalParameterObject );
        doubleParameters.push_back( phobosGravitationalParameterObject );
        doubleParameters.push_back( marsRotationRate );
        doubleParameters.push_back( phobosRotationRate );

        std::vector< std::shared_ptr< EstimatableParameter< Eigen::VectorXd > > > vectorParameters;
        vectorParameters.push_back( marsSineCoefficients );
        vectorParameters.push_back( marsCosineCoefficients );
        vectorParameters.push_back( phobosSineCoefficients );
        vectorParameters.push_back( phobosCosineCoefficients );
        vectorParameters.push_back( marsPolePosition );
        vectorParameters.push_back( phobosPolePosition );

        std::shared_ptr< EstimatableParameterSet< double > > parameterSet =
                std::make_shared< EstimatableParameterSet< double > >( doubleParameters, vectorParameters );

        // Create acceleration partial object.
        std::shared_ptr< MutualSphericalHarmonicsGravityPartial > accelerationPartial =
                std::dynamic_pointer_cast< MutualSphericalHarmonicsGravityPartial > (
                    createAnalyticalAccelerationPartial(
                        accelerationModel, std::make_pair( "Phobos", phobos ),
                        std::make_pair( "Mars", mars ), bodies, parameterSet ) );


        // Calculate analytical partials.
        accelerationModel->updateMembers( 1.0E6 );
        accelerationPartial->update( 1.0E6 );

        // Get analytical partials w.r.t. positions and velocities
        Eigen::MatrixXd partialWrtMarsPosition = Eigen::Matrix3d::Zero( );
        accelerationPartial->wrtPositionOfAcceleratingBody( partialWrtMarsPosition.block( 0, 0, 3, 3 ) );
        Eigen::MatrixXd partialWrtMarsVelocity = Eigen::Matrix3d::Zero( );
        accelerationPartial->wrtVelocityOfAcceleratingBody( partialWrtMarsVelocity.block( 0, 0, 3, 3 ) );
        Eigen::MatrixXd partialWrtPhobosPosition = Eigen::Matrix3d::Zero( );
        accelerationPartial->wrtPositionOfAcceleratedBody( partialWrtPhobosPosition.block( 0, 0, 3, 3 ) );
        Eigen::MatrixXd partialWrtPhobosVelocity = Eigen::Matrix3d::Zero( );
        accelerationPartial->wrtVelocityOfAcceleratedBody( partialWrtPhobosVelocity.block( 0, 0, 3, 3 ) );

        // Get analytical partials w.r.t. parameters
        Eigen::Vector3d partialWrtMarsGravitationalParameter = accelerationPartial->wrtParameter(
                    marsGravitationalParameterObject );
        Eigen::Vector3d partialWrtPhobosGravitationalParameter = accelerationPartial->wrtParameter(
                    phobosGravitationalParameterObject );
        Eigen::Vector3d partialWrtMarsRotationRate = accelerationPartial->wrtParameter(
                    marsRotationRate );
        Eigen::MatrixXd partialWrtMarsPolePosition = accelerationPartial->wrtParameter(
                    marsPolePosition );
        Eigen::Vector3d partialWrtPhobosRotationRate = accelerationPartial->wrtParameter(
                    phobosRotationRate );
        Eigen::MatrixXd partialWrtPhobosPolePosition = accelerationPartial->wrtParameter(
                    phobosPolePosition );
        Eigen::MatrixXd partialWrtMarsCosineCoefficients = accelerationPartial->wrtParameter(
                    marsCosineCoefficients );
        Eigen::MatrixXd partialWrtMarsSineCoefficients = accelerationPartial->wrtParameter(
                    marsSineCoefficients );
        Eigen::MatrixXd partialWrtPhobosCosineCoefficients = accelerationPartial->wrtParameter(
                    phobosCosineCoefficients );
        Eigen::MatrixXd partialWrtPhobosSineCoefficients = accelerationPartial->wrtParameter(
                    phobosSineCoefficients );

        // Set perturbations in position and velocity for numerical partial
        Eigen::Vector3d positionPerturbation;
        positionPerturbation << 100.0, 100.0, 100.0;
        Eigen::Vector3d velocityPerturbation;
        velocityPerturbation << 1.0, 1.0, 1.0;

        // Calculate numerical partials.
        Eigen::Matrix3d testPartialWrtPhobosPosition = calculateAccelerationWrtStatePartials(
                    phobosStateSetFunction, accelerationModel, phobos->getState( ), positionPerturbation, 0 );
        Eigen::Matrix3d testPartialWrtPhobosVelocity = calculateAccelerationWrtStatePartials(
                    phobosStateSetFunction, accelerationModel, phobos->getState( ), velocityPerturbation, 3 );
        Eigen::Matrix3d testPartialWrtMarsPosition = calculateAccelerationWrtStatePartials(
                    marsStateSetFunction, accelerationModel, mars->getState( ), positionPerturbation, 0 );
        Eigen::Matrix3d testPartialWrtMarsVelocity = calculateAccelerationWrtStatePartials(
                    marsStateSetFunction, accelerationModel, mars->getState( ), velocityPerturbation, 3 );

        Eigen::Vector3d testPartialWrtMarsGravitationalParameter = calculateAccelerationWrtParameterPartials(
                    marsGravitationalParameterObject, accelerationModel, 1.0E14 );
        Eigen::Vector3d testPartialWrtPhobosGravitationalParameter = calculateAccelerationWrtParameterPartials(
                    phobosGravitationalParameterObject, accelerationModel, 1.0E14 );

        Eigen::Vector3d testPartialWrtMarsRotationRate = calculateAccelerationWrtParameterPartials(
                    marsRotationRate, accelerationModel, 1.0E-9, &emptyFunction, 1.0E6, std::bind(
                        &Body::setCurrentRotationToLocalFrameFromEphemeris, mars, std::placeholders::_1 ) );
        Eigen::MatrixXd testPartialWrtPolePosition = calculateAccelerationWrtParameterPartials(
                    marsPolePosition, accelerationModel, ( Eigen::VectorXd( 2 ) << 1.0E-4, 1.0E-4 ).finished( ),
                    &emptyFunction, 1.0E6, std::bind(
                        &Body::setCurrentRotationToLocalFrameFromEphemeris, mars, std::placeholders::_1 ) );

        Eigen::Vector3d testPartialWrtPhobosRotationRate = calculateAccelerationWrtParameterPartials(
                    phobosRotationRate, accelerationModel, 1.0E-8, &emptyFunction, 1.0E6, std::bind(
                        &Body::setCurrentRotationToLocalFrameFromEphemeris, phobos, std::placeholders::_1 ) );
        Eigen::MatrixXd testPartialWrtPhobosPolePosition = calculateAccelerationWrtParameterPartials(
                    phobosPolePosition, accelerationModel, ( Eigen::VectorXd( 2 ) << 1.0E-3, 1.0E-3 ).finished( ),
                    &emptyFunction, 1.0E6, std::bind(
                        &Body::setCurrentRotationToLocalFrameFromEphemeris, phobos, std::placeholders::_1 ) );

        Eigen::MatrixXd testPartialWrtMarsCosineCoefficients = calculateAccelerationWrtParameterPartials(
                    marsCosineCoefficients, accelerationModel, Eigen::VectorXd::Constant(
                        marsCosineCoefficients->getParameterValue( ).size( ), 1, 1.0 ) );
        Eigen::MatrixXd testPartialWrtMarsSineCoefficients = calculateAccelerationWrtParameterPartials(
                    marsSineCoefficients, accelerationModel, Eigen::VectorXd::Constant(
                        marsSineCoefficients->getParameterValue( ).size( ), 1, 1.0 ) );
        Eigen::MatrixXd testPartialWrtPhobosCosineCoefficients = calculateAccelerationWrtParameterPartials(
                    phobosCosineCoefficients, accelerationModel, Eigen::VectorXd::Constant(
                        phobosCosineCoefficients->getParameterValue( ).size( ), 1, 1.0 ) );
        Eigen::MatrixXd testPartialWrtPhobosSineCoefficients = calculateAccelerationWrtParameterPartials(
                    phobosSineCoefficients, accelerationModel, Eigen::VectorXd::Constant(
                        phobosSineCoefficients->getParameterValue( ).size( ), 1, 1.0 ) );

        // Compare numerical and analytical partials of position and velocity partials.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtMarsPosition, partialWrtMarsPosition, 1.0e-9 );
        for( unsigned int i = 0; i < 3; i++ )
        {
            for( unsigned int j = 0; j < 3; j++ )
            {
                BOOST_CHECK_SMALL( std::fabs( testPartialWrtMarsVelocity( i, j ) - partialWrtMarsVelocity( i, j ) ),
                                   std::numeric_limits< double >::epsilon( ) );
                BOOST_CHECK_SMALL( std::fabs( partialWrtMarsVelocity( i, j ) ),
                                   std::numeric_limits< double >::epsilon( ) );
            }

        }

        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtPhobosPosition, partialWrtPhobosPosition, 1.0e-9 );
        for( unsigned int i = 0; i < 3; i++ )
        {
            for( unsigned int j = 0; j < 3; j++ )
            {
                BOOST_CHECK_SMALL( std::fabs( testPartialWrtPhobosVelocity( i, j ) - partialWrtPhobosVelocity( i, j ) ),
                                   std::numeric_limits< double >::epsilon( ) );
                BOOST_CHECK_SMALL( std::fabs( partialWrtPhobosVelocity( i, j ) ),
                                   std::numeric_limits< double >::epsilon( ) );
            }
        }

        // Compare numerical and analytical partials of gravitational parameters
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtMarsGravitationalParameter,
                                           partialWrtMarsGravitationalParameter, 1.0e-14 );

        if( testCase == 0 )
        {
            for( unsigned int j = 0; j < 3; j++ )
            {
                BOOST_CHECK_SMALL( std::fabs( testPartialWrtPhobosGravitationalParameter( j ) -
                                              partialWrtPhobosGravitationalParameter( j ) ),
                                   std::numeric_limits< double >::epsilon( ) );
                BOOST_CHECK_SMALL( std::fabs( partialWrtPhobosGravitationalParameter( j ) ),
                                   std::numeric_limits< double >::epsilon( ) );
            }
        }
        else
        {
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtPhobosGravitationalParameter,
                                               partialWrtPhobosGravitationalParameter, 1.0e-14 );

        }

        // Compare numerical and analytical partials of rotation parameters
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtMarsRotationRate, partialWrtMarsRotationRate, 1.0e-5 );
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtPolePosition, partialWrtMarsPolePosition, 1.0e-5 );

        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtPhobosRotationRate, partialWrtPhobosRotationRate, 1.0e-5 );
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtPhobosPolePosition, partialWrtPhobosPolePosition, 1.0e-5 );

        // Compare numerical and analytical partials of gravity field coefficients
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtMarsSineCoefficients, partialWrtMarsSineCoefficients, 1.0E-10 );
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtMarsCosineCoefficients, partialWrtMarsCosineCoefficients, 1.0E-10 );

        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtPhobosSineCoefficients, partialWrtPhobosSineCoefficients, 1.0E-8 );
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtPhobosCosineCoefficients, partialWrtPhobosCosineCoefficients, 1.0E-10 );
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}

}

