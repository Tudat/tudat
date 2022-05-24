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

#include <Eigen/Core>

#include "tudat/astro/aerodynamics/exponentialAtmosphere.h"

#if TUDAT_BUILD_WITH_NRLMSISE
#include "tudat/astro/aerodynamics/nrlmsise00Atmosphere.h"
#include "tudat/astro/aerodynamics/nrlmsise00InputFunctions.h"
#endif
#include "tudat/astro/aerodynamics/tabulatedAtmosphere.h"
#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/astro/basic_astro/unitConversions.h"
#include "tudat/astro/basic_astro/geodeticCoordinateConversions.h"
#include "tudat/astro/ephemerides/approximatePlanetPositions.h"
#include "tudat/astro/ephemerides/tabulatedEphemeris.h"
#include "tudat/astro/ephemerides/simpleRotationalEphemeris.h"
#include "tudat/astro/ephemerides/itrsToGcrsRotationModel.h"
#include "tudat/astro/gravitation/centralGravityModel.h"
#include "tudat/astro/gravitation/timeDependentSphericalHarmonicsGravityField.h"
#include "tudat/astro/gravitation/basicSolidBodyTideGravityFieldVariations.h"
#include "tudat/basics/testMacros.h"
#include "tudat/astro/ephemerides/synchronousRotationalEphemeris.h"

#include "tudat/interface/spice/spiceEphemeris.h"
#include "tudat/astro/gravitation/triAxialEllipsoidGravity.h"
#include "tudat/io/basicInputOutput.h"
#include "tudat/io/matrixTextFileReader.h"
#include "tudat/io/solarActivityData.h"
#include "tudat/io/parseSolarActivityData.h"
#include "tudat/math/basic/coordinateConversions.h"
#include "tudat/math/interpolators/lagrangeInterpolator.h"
#include "tudat/simulation/environment_setup/createAtmosphereModel.h"
#include "tudat/simulation/environment_setup/createEphemeris.h"
#include "tudat/simulation/environment_setup/createGravityField.h"
#include "tudat/simulation/environment_setup/createRotationModel.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"

namespace tudat
{

namespace unit_tests
{

using namespace simulation_setup;
using namespace input_output;
using namespace reference_frames;
using namespace basic_astrodynamics;

BOOST_AUTO_TEST_SUITE( test_environment_model_setup )

//! Test set up of atmosphere environment models.
BOOST_AUTO_TEST_CASE( test_atmosphereModelSetup )
{
    // Create settings for tabulated atmosphere.
    std::string atmosphereTableFile =
            tudat::paths::getAtmosphereTablesPath( ) + "/USSA1976Until100kmPer100mUntil1000kmPer1000m.dat";
    std::shared_ptr< TabulatedAtmosphereSettings > tabulatedAtmosphereSettings =
            std::make_shared< TabulatedAtmosphereSettings >( atmosphereTableFile );

    // Create settings for exponential atmosphere
    double densityScaleHeight = 8.0E3;
    double constantTemperature = 270.0;
    double densityAtZeroAltitude = 1.225;
    double specificGasConstant = 287.1;
    std::shared_ptr< ExponentialAtmosphereSettings > exponentialAtmosphereSettings =
            std::make_shared< ExponentialAtmosphereSettings >(
                densityScaleHeight, constantTemperature,
                densityAtZeroAltitude, specificGasConstant );


    // Create atmpshere models using setup function
    std::shared_ptr< aerodynamics::AtmosphereModel > exponentialAtmosphere =
            createAtmosphereModel( exponentialAtmosphereSettings, "Earth" );
    std::shared_ptr< aerodynamics::AtmosphereModel > tabulatedAtmosphere =
            createAtmosphereModel( tabulatedAtmosphereSettings, "Earth" );

    // Create atmosphere models manually.
    aerodynamics::TabulatedAtmosphere manualTabulatedAtmosphere( atmosphereTableFile );
    aerodynamics::ExponentialAtmosphere manualExponentialAtmosphere(
                densityScaleHeight, constantTemperature, densityAtZeroAltitude,
                specificGasConstant );

    // Verify equivalence of automatically set up and manual models.
    BOOST_CHECK_EQUAL( manualTabulatedAtmosphere.getDensity( 32.0, 0.0, 0.0, 0.0 ),
                       tabulatedAtmosphere->getDensity( 32.0, 0.0, 0.0, 0.0 ) );
    BOOST_CHECK_EQUAL( manualTabulatedAtmosphere.getPressure( 32.0, 0.0, 0.0, 0.0 ),
                       tabulatedAtmosphere->getPressure( 32.0, 0.0, 0.0, 0.0 ) );
    BOOST_CHECK_EQUAL( manualTabulatedAtmosphere.getTemperature( 32.0, 0.0, 0.0, 0.0 ),
                       tabulatedAtmosphere->getTemperature( 32.0, 0.0, 0.0, 0.0 ) );

    BOOST_CHECK_EQUAL( manualExponentialAtmosphere.getDensity( 32.0, 0.0, 0.0, 0.0 ),
                       exponentialAtmosphere->getDensity( 32.0, 0.0, 0.0, 0.0 ) );
    BOOST_CHECK_EQUAL( manualExponentialAtmosphere.getPressure( 32.0, 0.0, 0.0, 0.0 ),
                       exponentialAtmosphere->getPressure( 32.0, 0.0, 0.0, 0.0 ) );
    BOOST_CHECK_EQUAL( manualExponentialAtmosphere.getTemperature( 32.0, 0.0, 0.0, 0.0 ),
                       exponentialAtmosphere->getTemperature( 32.0, 0.0, 0.0, 0.0 ) );

#if TUDAT_BUILD_WITH_NRLMSISE
    std::shared_ptr< AtmosphereSettings > nrlmsise00AtmosphereSettings;
    for( int atmosphereTest = 0; atmosphereTest < 2; atmosphereTest++ )
    {
        if( atmosphereTest == 0 )
        {
            nrlmsise00AtmosphereSettings =
                    std::make_shared< AtmosphereSettings >( nrlmsise00 );
        }
        else
        {
            nrlmsise00AtmosphereSettings =
                    std::make_shared< NRLMSISE00AtmosphereSettings >(
                        paths::getSpaceWeatherDataPath( ) + "/sw19571001.txt" );
        }
        std::shared_ptr< aerodynamics::AtmosphereModel > nrlmsiseAtmosphere =
                createAtmosphereModel( nrlmsise00AtmosphereSettings, "Earth" );

        // Compute properties using NRLMSISE00
        double julianDaysSinceJ2000 = convertCalendarDateToJulianDay( 2005, 5, 3, 12, 32, 32.3 ) -
                basic_astrodynamics::JULIAN_DAY_ON_J2000;
        nrlmsiseAtmosphere->getDensity( 150.0E3, 1.0, 0.1, julianDaysSinceJ2000 * physical_constants::JULIAN_DAY );

        // Check if input to NRLMSISE00 is correctly computed (actual density computations tested in dedicated test).
        aerodynamics::NRLMSISE00Input nrlMSISE00Input = std::dynamic_pointer_cast< aerodynamics::NRLMSISE00Atmosphere >(
                    nrlmsiseAtmosphere )->getNRLMSISE00Input( );
        BOOST_CHECK_EQUAL( nrlMSISE00Input.year, 2005 );
        BOOST_CHECK_EQUAL( nrlMSISE00Input.dayOfTheYear, 31 + 28 + 31 + 30 + 3 );
        BOOST_CHECK_SMALL( nrlMSISE00Input.secondOfTheDay - ( 12.0 * 3600.0 + 32.0 * 60.0 + 32.3 ), 1.0E-3 );

        BOOST_CHECK_SMALL( nrlMSISE00Input.f107 - 112.3, 1.0E-14 );
        BOOST_CHECK_SMALL( nrlMSISE00Input.f107a - 93.3, 1.0E-14 );
        BOOST_CHECK_SMALL( nrlMSISE00Input.apDaily - 9.0, 1.0E-14 );
    }
#endif

}

Eigen::Vector6d computeCustomState(
        const double time, const double angularVelocity, const double radius, const double speed )
{
    Eigen::Vector6d currentState = Eigen::Vector6d::Zero( );
    currentState( 0 ) = radius * cos( angularVelocity * time );
    currentState( 1 ) = radius * sin( angularVelocity * time );

    currentState( 3 ) = -speed * sin( angularVelocity * time );
    currentState( 4 ) = speed * cos( angularVelocity * time );

    return currentState;

}

//! Test set up of ephemeris environment models.
BOOST_AUTO_TEST_CASE( test_ephemerisSetup )
{
    spice_interface::loadStandardSpiceKernels( );

    {
        // Create settings for approximate planet positions.
        std::string bodyIdentifier = "Mars";
        bool useCircularCoplanarApproximation = 0;
        std::shared_ptr< ApproximateJplEphemerisSettings > approximateEphemerisSettings =
                std::make_shared< ApproximateJplEphemerisSettings >(
                    bodyIdentifier, useCircularCoplanarApproximation );

        // Create ephemeris using setup function.
        std::shared_ptr< ephemerides::Ephemeris > approximateEphemeris =
                createBodyEphemeris( approximateEphemerisSettings, "Earth" );

        // Create manual ephemeris.
        ephemerides::ApproximateJplEphemeris manualApproximateEphemeris(
                    bodyIdentifier );

        // Verify equivalence of automatically set up and manual models.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                    ( manualApproximateEphemeris.getCartesianState( 1.0E7 ) ),
                    ( approximateEphemeris->getCartesianState( 1.0E7 ) ),
                    std::numeric_limits< double >::epsilon( ) );
    }

    {
        // Create spice ephemeris.
        std::shared_ptr< EphemerisSettings > spiceEphemerisSettings =
                std::make_shared< DirectSpiceEphemerisSettings >( "Earth", "J2000" );
        std::shared_ptr< ephemerides::Ephemeris > spiceEphemeris =
                createBodyEphemeris( spiceEphemerisSettings, "Moon" );

        // Compare spice ephemeris against direct spice state.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                    ( spice_interface::getBodyCartesianStateAtEpoch(
                          "Moon", "Earth", "J2000", "None", 1.0E7 ) ),
                    ( spiceEphemeris->getCartesianState( 1.0E7 ) ),
                    std::numeric_limits< double >::epsilon( ) );
    }

    {
        // Create scaled spice ephemeris.
        std::shared_ptr< EphemerisSettings > spiceEphemerisSettings =
                std::make_shared< DirectSpiceEphemerisSettings >( "Earth", "J2000" );

        Eigen::Vector6d absoluteScaling = ( Eigen::Vector6d( ) << 1.0, 2.0, 3.0, 4.0, 5.0, 6.0 ).finished( );
        std::shared_ptr< EphemerisSettings > scaledEphemerisSettings =
                std::make_shared< ScaledEphemerisSettings >( spiceEphemerisSettings, absoluteScaling, true );
        std::shared_ptr< ephemerides::Ephemeris > scaledSpiceEphemeris =
                createBodyEphemeris( scaledEphemerisSettings, "Moon" );

        // Compare spice ephemeris against direct spice state.
        Eigen::Vector6d spiceState = spice_interface::getBodyCartesianStateAtEpoch(
                    "Moon", "Earth", "J2000", "None", 1.0E7 );
        Eigen::Vector6d scaledState = scaledSpiceEphemeris->getCartesianState( 1.0E7 );

        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                    ( spiceState + absoluteScaling ),
                    ( scaledState ),
                    std::numeric_limits< double >::epsilon( ) );
    }

    {
        // Create scaled spice ephemeris.
        std::shared_ptr< EphemerisSettings > spiceEphemerisSettings =
                std::make_shared< DirectSpiceEphemerisSettings >( "Earth", "J2000" );

        Eigen::Vector6d absoluteScaling = ( Eigen::Vector6d( ) << 1.0, 2.0, 3.0, 4.0, 5.0, 6.0 ).finished( );
        std::shared_ptr< EphemerisSettings > scaledEphemerisSettings =
                std::make_shared< ScaledEphemerisSettings >( spiceEphemerisSettings, absoluteScaling, false );
        std::shared_ptr< ephemerides::Ephemeris > scaledSpiceEphemeris =
                createBodyEphemeris( scaledEphemerisSettings, "Moon" );

        // Compare spice ephemeris against direct spice state.
        Eigen::Vector6d spiceState = spice_interface::getBodyCartesianStateAtEpoch(
                    "Moon", "Earth", "J2000", "None", 1.0E7 );
        Eigen::Vector6d scaledState = scaledSpiceEphemeris->getCartesianState( 1.0E7 );


        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                    ( scaledState.cwiseQuotient( spiceState ) ),
                    ( absoluteScaling ),
                    std::numeric_limits< double >::epsilon( ) );
    }



    {
        // Create custom ephemeris
        double angularVelocity = 2.0 * tudat::mathematical_constants::PI / ( 2.0 * 3600.0 );
        double radius = 8000.0E3;
        double speed = 5000.0;

        std::shared_ptr< EphemerisSettings > customEphemerisSettings =
                std::make_shared< CustomEphemerisSettings >(
                    std::bind( &computeCustomState, std::placeholders::_1,  angularVelocity, radius, speed ),
                    "Earth", "J2000" );
        std::shared_ptr< ephemerides::Ephemeris > customEphemeris =
                createBodyEphemeris( customEphemerisSettings, "Satellite" );

        double testTime = 4.0E8;

        Eigen::Vector6d currentState = customEphemeris->getCartesianState(
                    testTime );
        Eigen::Vector6d currentTestState = computeCustomState(
                    testTime, angularVelocity, radius, speed );

        BOOST_CHECK_CLOSE_FRACTION( currentState( 0 ), currentTestState( 0 ), 2.0 * std::numeric_limits< double >::epsilon( ) );
        BOOST_CHECK_CLOSE_FRACTION( currentState( 1 ), currentTestState( 1 ), 2.0 * std::numeric_limits< double >::epsilon( ) );
        BOOST_CHECK_SMALL( currentState( 2 ),  2.0 * std::numeric_limits< double >::epsilon( ) );
        BOOST_CHECK_CLOSE_FRACTION( currentState( 3 ), currentTestState( 3 ), 2.0 * std::numeric_limits< double >::epsilon( ) );
        BOOST_CHECK_CLOSE_FRACTION( currentState( 4 ), currentTestState( 4 ), 2.0 * std::numeric_limits< double >::epsilon( ) );
        BOOST_CHECK_SMALL( currentState( 5 ), 2.0 * std::numeric_limits< double >::epsilon( ) );

    }

    {
        // Create tabulated spice ephemeris
        std::shared_ptr< EphemerisSettings > spiceEphemerisSettings =
                std::make_shared< InterpolatedSpiceEphemerisSettings >(
                    1.0E7 - 50.0 * 600.0, 1.0E7 + 50.0 * 600.0, 600.0, "Earth", "J2000" );
        std::shared_ptr< ephemerides::Ephemeris > spiceEphemeris =
                createBodyEphemeris( spiceEphemerisSettings, "Moon" );

        // Compare tabulated spice ephemeris against direct spice state on node point.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                    ( spice_interface::getBodyCartesianStateAtEpoch(
                          "Moon", "Earth", "J2000", "None", 1.0E7 ) ),
                    ( spiceEphemeris->getCartesianState( 1.0E7 ) ),
                    std::numeric_limits< double >::epsilon( ) );

        // Manually create table of states from spice
        std::map< double, Eigen::Vector6d > tabulatedStates;
        double currentTime = 1.0E7 - 50.0 * 600.0;
        while( currentTime <= 1.0E7 + 50.0 * 600.0 )
        {
            tabulatedStates[ currentTime ] = spice_interface::getBodyCartesianStateAtEpoch(
                        "Moon", "Earth", "J2000", "None", currentTime );
            currentTime += 600.0;
        }

        // Create tabulated ephemeris.
        std::shared_ptr< EphemerisSettings > tabulatedEphemerisSettings =
                std::make_shared< TabulatedEphemerisSettings >(
                    tabulatedStates, "Earth", "J2000" );
        std::shared_ptr< ephemerides::Ephemeris > tabulatedEphemeris =
                createBodyEphemeris( tabulatedEphemerisSettings, "Moon" );

        // Manually create tabulated ephemeris.
        std::shared_ptr< ephemerides::Ephemeris > manualTabulatedEphemeris =
                std::make_shared< ephemerides::TabulatedCartesianEphemeris< > >(
                    std::make_shared< interpolators::LagrangeInterpolator
                    < double, Eigen::Vector6d > >( tabulatedStates, 6 ),
                    "Earth", "J2000" );

        // Compare ephemerides away from node point.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                    ( spiceEphemeris->getCartesianState( 1.0E7 + 110.0 ) ),
                    ( tabulatedEphemeris->getCartesianState( 1.0E7 + 110.0 ) ),
                    std::numeric_limits< double >::epsilon( ) );
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                    ( spiceEphemeris->getCartesianState( 1.0E7 + 110.0 ) ),
                    ( manualTabulatedEphemeris->getCartesianState( 1.0E7 + 110.0 ) ),
                    std::numeric_limits< double >::epsilon( ) );
    }


}


//! Test set up of gravity field model environment models.
BOOST_AUTO_TEST_CASE( test_gravityFieldSetup )
{
    // Load Spice kernel
    spice_interface::loadStandardSpiceKernels( );

    // Create settings for spice central gravity field model.
    std::shared_ptr< GravityFieldSettings > spiceCentralGravityFieldSettings =
            std::make_shared< GravityFieldSettings >( central_spice );

    // Create spice central gravity field model from setup function.
    std::shared_ptr< gravitation::GravityFieldModel > spiceCentralGravityField =
            createGravityFieldModel( spiceCentralGravityFieldSettings, "Venus" );

    // Check correct creation of gravity field.
    BOOST_CHECK_EQUAL(
                ( spice_interface::getBodyGravitationalParameter( "Venus" ) ),
                ( spiceCentralGravityField->getGravitationalParameter( ) ) );

    // Settings for spherical harmonic acceleration.
    double gravitationalParameter = 398600.4418E9;
    Eigen::Vector3d testPosition( 7.0e6, 8.0e6, 9.0e6 );
    Eigen::MatrixXd cosineCoefficients =
            ( Eigen::MatrixXd( 6, 6 ) <<
              1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
              0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
              -4.841651437908150e-4, -2.066155090741760e-10, 2.439383573283130e-6, 0.0, 0.0, 0.0,
              9.571612070934730e-7, 2.030462010478640e-6, 9.047878948095281e-7,
              7.213217571215680e-7, 0.0, 0.0, 5.399658666389910e-7, -5.361573893888670e-7,
              3.505016239626490e-7, 9.908567666723210e-7, -1.885196330230330e-7, 0.0,
              6.867029137366810e-8, -6.292119230425290e-8, 6.520780431761640e-7,
              -4.518471523288430e-7, -2.953287611756290e-7, 1.748117954960020e-7
              ).finished( );
    Eigen::MatrixXd sineCoefficients =
            ( Eigen::MatrixXd( 6, 6 ) <<
              0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
              0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
              0.0, 1.384413891379790e-9, -1.400273703859340e-6, 0.0, 0.0, 0.0,
              0.0, 2.482004158568720e-7, -6.190054751776180e-7, 1.414349261929410e-6, 0.0, 0.0,
              0.0, -4.735673465180860e-7, 6.624800262758290e-7, -2.009567235674520e-7,
              3.088038821491940e-7, 0.0, 0.0, -9.436980733957690e-8, -3.233531925405220e-7,
              -2.149554083060460e-7, 4.980705501023510e-8, -6.693799351801650e-7
              ).finished( );

    // Create settings for central gravity field.
    std::shared_ptr< CentralGravityFieldSettings > centralGravityFieldSettings =
            std::make_shared< CentralGravityFieldSettings >( gravitationalParameter );

    // Create central gravity field with setup function.
    std::shared_ptr< gravitation::GravityFieldModel > centralGravityField =
            createGravityFieldModel( centralGravityFieldSettings, "Earth" );

    // Create central gravity field manually.
    gravitation::GravityFieldModel manualCentralGravityField(
                gravitationalParameter );

    // Verify equivalence of automatically set up and manual models.
    BOOST_CHECK_EQUAL(
                ( manualCentralGravityField.getGravitationalParameter( ) ),
                ( centralGravityField->getGravitationalParameter( ) ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                ( manualCentralGravityField.getGradientOfPotential( testPosition ) ),
                ( centralGravityField->getGradientOfPotential( testPosition ) ),
                std::numeric_limits< double >::epsilon( ) );

    // Create settings for sh gravity field.
    std::shared_ptr< SphericalHarmonicsGravityFieldSettings > shGravityFieldSettings =
            std::make_shared< SphericalHarmonicsGravityFieldSettings >(
                gravitationalParameter, 6378.0E3, cosineCoefficients, sineCoefficients,
                "Earth_fixed" );

    // Create sh gravity field with setup function.
    std::shared_ptr< gravitation::GravityFieldModel > shGravityField =
            createGravityFieldModel( shGravityFieldSettings, "Earth" );

    // Create sh gravity field manually.
    gravitation::SphericalHarmonicsGravityField manualShGravityField(
                gravitationalParameter, 6378.0E3, cosineCoefficients, sineCoefficients );

    // Verify equivalence of automatically set up and manual models.
    BOOST_CHECK_EQUAL(
                ( manualShGravityField.getGravitationalParameter( ) ),
                ( shGravityField->getGravitationalParameter( ) ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                ( manualShGravityField.getGradientOfPotential( testPosition ) ),
                ( shGravityField->getGradientOfPotential( testPosition ) ),
                std::numeric_limits< double >::epsilon( ) );

    std::shared_ptr< gravitation::SphericalHarmonicsGravityField > defaultEarthField =
            std::dynamic_pointer_cast< gravitation::SphericalHarmonicsGravityField >(
                createGravityFieldModel( getDefaultGravityFieldSettings(
                                             "Earth", TUDAT_NAN, TUDAT_NAN ), "Earth" ) );
    BOOST_CHECK_EQUAL(
                ( defaultEarthField->getGravitationalParameter( ) ), ( 0.3986004418E15 ) );
    BOOST_CHECK_EQUAL(
                ( defaultEarthField->getReferenceRadius( ) ), ( 6378137.0 ) );
    BOOST_CHECK_EQUAL(
                ( defaultEarthField->getCosineCoefficients( ).rows( ) ), 361 );
    BOOST_CHECK_EQUAL(
                ( defaultEarthField->getCosineCoefficients( ).cols( ) ), 361 );
    BOOST_CHECK_EQUAL(
                ( defaultEarthField->getSineCoefficients( ).rows( ) ), 361 );
    BOOST_CHECK_EQUAL(
                ( defaultEarthField->getSineCoefficients( ).cols( ) ), 361 );
    BOOST_CHECK_EQUAL(
                ( defaultEarthField->getCosineCoefficients( )( 2, 0 ) ), -0.484165371736E-03 );
    BOOST_CHECK_EQUAL(
                ( defaultEarthField->getCosineCoefficients( )( 5, 3 ) ), -0.451955406071E-06 );
    BOOST_CHECK_EQUAL(
                ( defaultEarthField->getSineCoefficients( )( 7, 1 ) ), 0.954336911867E-07 );

    std::shared_ptr< gravitation::SphericalHarmonicsGravityField > defaultMoonField =
            std::dynamic_pointer_cast< gravitation::SphericalHarmonicsGravityField >(
                createGravityFieldModel( getDefaultGravityFieldSettings(
                                             "Moon", TUDAT_NAN, TUDAT_NAN ), "Moon" ) );
    BOOST_CHECK_EQUAL(
                ( defaultMoonField->getGravitationalParameter( ) ), ( 0.4902800238000000E+13 ) );
    BOOST_CHECK_EQUAL(
                ( defaultMoonField->getReferenceRadius( ) ), ( 0.17380E+07 ) );
    BOOST_CHECK_EQUAL(
                ( defaultMoonField->getCosineCoefficients( ).rows( ) ), 201 );
    BOOST_CHECK_EQUAL(
                ( defaultMoonField->getCosineCoefficients( ).cols( ) ), 201 );
    BOOST_CHECK_EQUAL(
                ( defaultMoonField->getSineCoefficients( ).rows( ) ), 201 );
    BOOST_CHECK_EQUAL(
                ( defaultMoonField->getSineCoefficients( ).cols( ) ), 201 );
    BOOST_CHECK_EQUAL(
                ( defaultMoonField->getCosineCoefficients( )( 5, 3 ) ), 0.5493176535439800E-06 );
    BOOST_CHECK_EQUAL(
                ( defaultMoonField->getSineCoefficients( )( 7, 1 ) ), -0.1744763377093700E-06 );

}

//! Test set up of triaxial ellipsoid gravity field model settings
BOOST_AUTO_TEST_CASE( test_triaxialEllipsoidGravityFieldSetup )
{
    double axisA = 26.8E3;
    double axisB = 22.4E3;
    double axisC = 18.4E3;

    double density = 2.7E3;

    std::shared_ptr< SphericalHarmonicsGravityFieldSettings > gravityFieldSettings =
            createHomogeneousTriAxialEllipsoidGravitySettings(
                axisA, axisB, axisC, density, 11, 11, "TestFrame" );

    // Check size of coefficient blocks
    BOOST_CHECK_EQUAL( gravityFieldSettings->getCosineCoefficients( ).rows( ), 12 );
    BOOST_CHECK_EQUAL( gravityFieldSettings->getCosineCoefficients( ).cols( ), 12 );
    BOOST_CHECK_EQUAL( gravityFieldSettings->getSineCoefficients( ).rows( ), 12 );
    BOOST_CHECK_EQUAL( gravityFieldSettings->getSineCoefficients( ).cols( ), 12 );

    // Check frame ID
    BOOST_CHECK_EQUAL( gravityFieldSettings->getAssociatedReferenceFrame( ), "TestFrame"  );

    // Compute expected low-degree unnormalized coefficients
    double referenceRadius = gravitation::calculateTriAxialEllipsoidReferenceRadius( axisA, axisB, axisC );

    double expectedC20 = 1.0 / ( 5.0 * referenceRadius * referenceRadius ) *
            ( axisC * axisC - ( axisA * axisA + axisB * axisB ) / 2.0 );
    double expectedC22 = 1.0 / ( 20.0 * referenceRadius * referenceRadius ) *
            ( axisA * axisA - axisB * axisB );
    double expectedC40 = 15.0 / 7.0 * ( expectedC20 * expectedC20 +
                                        2.0 * expectedC22 * expectedC22 );
    double expectedC42 = 5.0 / 7.0 * ( expectedC20 * expectedC22 );
    double expectedC44 = 5.0 / 28.0 * ( expectedC22 * expectedC22 );

    // Test low-degree coefficients (normalizing them first)
    BOOST_CHECK_CLOSE_FRACTION(
                expectedC20, gravityFieldSettings->getCosineCoefficients( )( 2, 0 ) *
                basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 0 ),
                2.0 * std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION(
                expectedC22, gravityFieldSettings->getCosineCoefficients( )( 2, 2 ) *
                basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 2 ),
                2.0 * std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION(
                expectedC40, gravityFieldSettings->getCosineCoefficients( )( 4, 0 ) *
                basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 4, 0 ),
                2.0 * std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION(
                expectedC42, gravityFieldSettings->getCosineCoefficients( )( 4, 2 ) *
                basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 4, 2 ),
                2.0 * std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION(
                expectedC44, gravityFieldSettings->getCosineCoefficients( )( 4, 4 ) *
                basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 4, 4 ),
                2.0 * std::numeric_limits< double >::epsilon( ) );

    // Check whether C-coefficients with odd degree or order, and all S-coefficients, are zero.
    for( unsigned int i = 0; i < 12; i++ )
    {
        for( unsigned int j = 0; j < 12; j++ )
        {
            if( ( i % 2 != 0 ) || ( j % 2 ) != 0 )
            {
                BOOST_CHECK_EQUAL( gravityFieldSettings->getCosineCoefficients( )( i, j ), 0.0 );
            }
            BOOST_CHECK_EQUAL( gravityFieldSettings->getSineCoefficients( )( i, j ), 0.0 );
        }
    }

    double gravitationalParameter = 4.0 / 3.0 * mathematical_constants::PI * axisA * axisB * axisC *
            density * physical_constants::GRAVITATIONAL_CONSTANT;
    BOOST_CHECK_CLOSE_FRACTION(
                gravitationalParameter, gravityFieldSettings->getGravitationalParameter( ),
                2.0 * std::numeric_limits< double >::epsilon( ) );
}




//! Test set up of gravity field model variations environment models.
BOOST_AUTO_TEST_CASE( test_gravityFieldVariationSetup )
{
    // Load Spice kernels
    spice_interface::loadStandardSpiceKernels( );

    // Settings for spherical harmonic acceleration.
    double gravitationalParameter = 398600.4418E9;
    double referenceRadius = 6378.0E3;
    Eigen::MatrixXd cosineCoefficients =
            ( Eigen::MatrixXd( 6, 6 ) <<
              1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
              0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
              -4.841651437908150e-4, -2.066155090741760e-10, 2.439383573283130e-6, 0.0, 0.0, 0.0,
              9.571612070934730e-7, 2.030462010478640e-6, 9.047878948095281e-7,
              7.213217571215680e-7, 0.0, 0.0, 5.399658666389910e-7, -5.361573893888670e-7,
              3.505016239626490e-7, 9.908567666723210e-7, -1.885196330230330e-7, 0.0,
              6.867029137366810e-8, -6.292119230425290e-8, 6.520780431761640e-7,
              -4.518471523288430e-7, -2.953287611756290e-7, 1.748117954960020e-7
              ).finished( );
    Eigen::MatrixXd sineCoefficients =
            ( Eigen::MatrixXd( 6, 6 ) <<
              0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
              0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
              0.0, 1.384413891379790e-9, -1.400273703859340e-6, 0.0, 0.0, 0.0,
              0.0, 2.482004158568720e-7, -6.190054751776180e-7, 1.414349261929410e-6, 0.0, 0.0,
              0.0, -4.735673465180860e-7, 6.624800262758290e-7, -2.009567235674520e-7,
              3.088038821491940e-7, 0.0, 0.0, -9.436980733957690e-8, -3.233531925405220e-7,
              -2.149554083060460e-7, 4.980705501023510e-8, -6.693799351801650e-7
              ).finished( );


    // Define body settings.
    BodyListSettings bodySettings;
    bodySettings.addSettings( getDefaultSingleBodySettings( "Moon", 0.0, 1.0E7 ), "Moon" );
    bodySettings.addSettings( getDefaultSingleBodySettings( "Sun", 0.0, 1.0E7 ), "Sun" );
    bodySettings.addSettings( getDefaultSingleBodySettings( "Earth", 0.0, 1.0E7 ), "Earth" );

    bodySettings.at( "Earth" )->gravityFieldSettings = std::make_shared< SphericalHarmonicsGravityFieldSettings >(
                gravitationalParameter, referenceRadius, cosineCoefficients, sineCoefficients, "IAU_Earth" );


    // Define settings for tidal variations.
    double loveNumber = 0.25;
    double loveNumberDegreeThree = 100.0;

    std::map< int, std::vector< std::complex< double > > > fullLoveNumberVector =
            getFullLoveNumbersVector( loveNumber, 3, 2 );
    double testTime = 0.5E7;


    Eigen::MatrixXd cosineCorrections1, cosineCorrections2, cosineCorrections3;
    Eigen::MatrixXd sineCorrections1, sineCorrections2, sineCorrections3;

    // Calculate non-interpolated corrections from two bodies in single object.
    {
        // Define correction settings
        std::vector< std::string > deformingBodies;
        deformingBodies.push_back( "Moon" );
        deformingBodies.push_back( "Sun" );

        bodySettings.at( "Earth" )->gravityFieldVariationSettings.push_back(
                    std::make_shared< BasicSolidBodyGravityFieldVariationSettings >(
                        deformingBodies, fullLoveNumberVector ) );

        // Create bodies
        SystemOfBodies bodies = createSystemOfBodies( bodySettings );

        // Update states.
        bodies.at( "Earth" )->setStateFromEphemeris( testTime );
        bodies.at( "Earth" )->setCurrentRotationalStateToLocalFrameFromEphemeris( testTime );
        bodies.at( "Sun" )->setStateFromEphemeris( testTime );
        bodies.at( "Moon" )->setStateFromEphemeris( testTime );

        // Update gravity field
        std::shared_ptr< gravitation::TimeDependentSphericalHarmonicsGravityField > earthGravityField =
                std::dynamic_pointer_cast< gravitation::TimeDependentSphericalHarmonicsGravityField >(
                    bodies.at( "Earth" )->getGravityFieldModel( ) );
        earthGravityField->update( testTime );

        // Retrieve corrections.
        cosineCorrections1 = earthGravityField->getCosineCoefficients( ) - cosineCoefficients;
        sineCorrections1 = earthGravityField->getSineCoefficients( ) - sineCoefficients;


    }

    // Calculate non-interpolated corrections from two bodies in separate objects.
    {
        bodySettings.at( "Earth" )->gravityFieldVariationSettings.clear( );

        // Define correction settings
        std::vector< std::string > deformingBodies;
        deformingBodies.push_back( "Moon" );
        bodySettings.at( "Earth" )->gravityFieldVariationSettings.push_back(
                    std::make_shared< BasicSolidBodyGravityFieldVariationSettings >(
                        deformingBodies, fullLoveNumberVector ) );

        deformingBodies.clear( );
        deformingBodies.push_back( "Sun" );
        bodySettings.at( "Earth" )->gravityFieldVariationSettings.push_back(
                    std::make_shared< BasicSolidBodyGravityFieldVariationSettings >(
                        deformingBodies, fullLoveNumberVector ) );

        // Create bodies
        SystemOfBodies bodies = createSystemOfBodies( bodySettings );
        

        // Update states.
        bodies.at( "Earth" )->setStateFromEphemeris( testTime );
        bodies.at( "Earth" )->setCurrentRotationalStateToLocalFrameFromEphemeris( testTime );
        bodies.at( "Sun" )->setStateFromEphemeris( testTime );
        bodies.at( "Moon" )->setStateFromEphemeris( testTime );

        // Update gravity field
        std::shared_ptr< gravitation::TimeDependentSphericalHarmonicsGravityField > earthGravityField =
                std::dynamic_pointer_cast< gravitation::TimeDependentSphericalHarmonicsGravityField >(
                    bodies.at( "Earth" )->getGravityFieldModel( ) );
        earthGravityField->update( testTime );

        // Retrieve corrections.
        cosineCorrections2 = earthGravityField->getCosineCoefficients( ) - cosineCoefficients;
        sineCorrections2 = earthGravityField->getSineCoefficients( ) - sineCoefficients;
    }

    {
        bodySettings.at( "Earth" )->gravityFieldVariationSettings.clear( );

        // Define correction settings
        std::vector< std::string > deformingBodies;
        deformingBodies.push_back( "Moon" );
        deformingBodies.push_back( "Sun" );
        bodySettings.at( "Earth" )->gravityFieldVariationSettings.push_back(
                    std::make_shared< BasicSolidBodyGravityFieldVariationSettings >(
                        deformingBodies, fullLoveNumberVector,
                        std::make_shared< ModelInterpolationSettings >(
                            0.25E7, 0.75E7, 600.0,
                            std::make_shared< interpolators::LagrangeInterpolatorSettings >( 8 ) ) ) );

        // Create bodies
        SystemOfBodies bodies = createSystemOfBodies( bodySettings );



        // Update gravity field
        std::shared_ptr< gravitation::TimeDependentSphericalHarmonicsGravityField > earthGravityField =
                std::dynamic_pointer_cast< gravitation::TimeDependentSphericalHarmonicsGravityField >(
                    bodies.at( "Earth" )->getGravityFieldModel( ) );

        earthGravityField->update( testTime );

        // Retrieve corrections.
        cosineCorrections3 = earthGravityField->getCosineCoefficients( ) - cosineCoefficients;
        sineCorrections3= earthGravityField->getSineCoefficients( ) - sineCoefficients;

    }

    // Mutually compare results of three methods.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( cosineCorrections1.block( 2, 0, 2, 3 ),  cosineCorrections2.block( 2, 0, 2, 3 ), 1.0E-10 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( cosineCorrections1.block( 2, 0, 2, 3 ),  cosineCorrections3.block( 2, 0, 2, 3 ), 1.0E-10 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( sineCorrections1.block( 2, 1, 2, 3 ),  sineCorrections2.block( 2, 1, 2, 3 ), 1.0E-10 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( sineCorrections1.block( 2, 1, 2, 3 ),  sineCorrections2.block( 2, 1, 2, 3 ), 1.0E-10 );

    // Check whether 0 coefficients are actually 0.
    for( unsigned int i = 0; i < 6; i++ )
    {
        for( unsigned int j = 0; j < 6; j++ )
        {
            if( i < 2 || i > 3 || j > 2 )
            {
                BOOST_CHECK_EQUAL( cosineCorrections1( i, j ), 0.0 );
                BOOST_CHECK_EQUAL( cosineCorrections2( i, j ), 0.0 );
                BOOST_CHECK_EQUAL( cosineCorrections3( i, j ), 0.0 );

                BOOST_CHECK_EQUAL( sineCorrections1( i, j ), 0.0 );
                BOOST_CHECK_EQUAL( sineCorrections2( i, j ), 0.0 );
                BOOST_CHECK_EQUAL( sineCorrections3( i, j ), 0.0 );
            }
            else if( j == 0 )
            {
                BOOST_CHECK_EQUAL( sineCorrections1( i, j ), 0.0 );
                BOOST_CHECK_EQUAL( sineCorrections2( i, j ), 0.0 );
                BOOST_CHECK_EQUAL( sineCorrections3( i, j ), 0.0 );
            }

        }
    }

    // Calculate corrections manually and compare against created results.
    std::pair< Eigen::MatrixXd, Eigen::MatrixXd > directMoonTide =
            gravitation::calculateSolidBodyTideSingleCoefficientSetCorrectionFromAmplitude(
                fullLoveNumberVector, 0.4902800238000000E+13 / gravitationalParameter,
                referenceRadius,  spice_interface::getBodyCartesianPositionAtEpoch(
                    "Moon", "Earth", "IAU_Earth", "None", testTime ), 3, 2 );
    std::pair< Eigen::MatrixXd, Eigen::MatrixXd > directSunTide =
            gravitation::calculateSolidBodyTideSingleCoefficientSetCorrectionFromAmplitude(
                fullLoveNumberVector, spice_interface::getBodyGravitationalParameter( "Sun" ) / gravitationalParameter,
                referenceRadius,  spice_interface::getBodyCartesianPositionAtEpoch(
                    "Sun", "Earth", "IAU_Earth", "None", testTime ), 3, 2 );

    for( unsigned int n = 2; n <= 3; n++ )
    {
        for( unsigned m = 0; m <=2; m++ )
        {
            BOOST_CHECK_SMALL( directMoonTide.first( n, m ) + directSunTide.first( n, m ) - cosineCorrections1( n, m ), 1.0E-18 );
            BOOST_CHECK_SMALL( directMoonTide.second( n, m ) + directSunTide.second( n, m ) - sineCorrections1( n, m ), 1.0E-18 );
        }
    }


    //! Test combinations of factory functions, check if results are identical
    std::vector< Eigen::MatrixXd > cosineCoefficientsFunctionTestDegreeTwo;
    std::vector< Eigen::MatrixXd > sineCoefficientsFunctionTestDegreeTwo;

    // Test for different factory functions
    for( int functionTest = 0; functionTest < 4; functionTest++ )
    {
        std::vector< Eigen::MatrixXd > cosineCoefficientsBodyTest;
        std::vector< Eigen::MatrixXd > sineCoefficientsBodyTest;

        // Test for separate or joint deforming bodies
        for( int bodyTest = 0; bodyTest < 3; bodyTest++ )
        {

            // Clear for current test
            bodySettings.at( "Earth" )->gravityFieldVariationSettings.clear( );

            // Define deforming bodies (one or two)
            std::vector< std::string > deformingBodies;
            if( bodyTest == 0 || bodyTest == 2 )
            { deformingBodies.push_back( "Moon" ); }
            if( bodyTest == 1 || bodyTest == 2  )
            { deformingBodies.push_back( "Sun" ); }

            // Create settings for each deforming body
            for( unsigned int i = 0; i < deformingBodies.size( ); i++ )
            {
                if( functionTest == 0 )
                {
                    bodySettings.at( "Earth" )->gravityFieldVariationSettings.push_back(
                                fixedSingleDegreeLoveNumberGravityFieldVariationSettings(
                                    deformingBodies.at( i ), loveNumber, 2 ) );
                }
                else if( functionTest == 1 )
                {
                    bodySettings.at( "Earth" )->gravityFieldVariationSettings.push_back(
                                fixedSingleDegreeLoveNumberGravityFieldVariationSettings(
                                    deformingBodies.at( i ), std::complex< double >( loveNumber, 0.0 ), 2 ) );
                }
                else if( functionTest == 2 )
                {
                    std::map< int, double > loveNumberMap;
                    loveNumberMap[ 2 ] = 0.25;
                    bodySettings.at( "Earth" )->gravityFieldVariationSettings.push_back(
                                fixedSingleDegreeLoveNumberGravityFieldVariationSettings(
                                    deformingBodies.at( i ), loveNumberMap ) );
                }
                else if( functionTest == 3 )
                {
                    std::map< int, std::complex< double > > loveNumberMap;
                    loveNumberMap[ 2 ] = std::complex< double >( loveNumber, 0.0 );
                    bodySettings.at( "Earth" )->gravityFieldVariationSettings.push_back(
                                fixedSingleDegreeLoveNumberGravityFieldVariationSettings(
                                    deformingBodies.at( i ), loveNumberMap ) );
                }
            }

            // Create bodies
            SystemOfBodies bodies = createSystemOfBodies( bodySettings );

            // Update states.
            bodies.at( "Earth" )->setStateFromEphemeris( testTime );
            bodies.at( "Earth" )->setCurrentRotationalStateToLocalFrameFromEphemeris( testTime );
            bodies.at( "Sun" )->setStateFromEphemeris( testTime );
            bodies.at( "Moon" )->setStateFromEphemeris( testTime );

            // Update gravity field
            std::shared_ptr< gravitation::TimeDependentSphericalHarmonicsGravityField > earthGravityField =
                    std::dynamic_pointer_cast< gravitation::TimeDependentSphericalHarmonicsGravityField >(
                        bodies.at( "Earth" )->getGravityFieldModel( ) );
            earthGravityField->update( testTime );

            // Retrieve corrections.
            cosineCoefficientsBodyTest.push_back(
                        earthGravityField->getCosineCoefficients( ) - cosineCoefficients );
            sineCoefficientsBodyTest.push_back(
                        earthGravityField->getSineCoefficients( ) - sineCoefficients );
        }


        // Test whether corrections from separate bodies, or two bodies together, give identical results
        Eigen::MatrixXd cosineDegreeTwoBodyDifference =
                cosineCoefficientsBodyTest.at( 2 ) -
                cosineCoefficientsBodyTest.at( 1 ) - cosineCoefficientsBodyTest.at( 0 );
        Eigen::MatrixXd sineDegreeTwoBodyDifference =
                sineCoefficientsBodyTest.at( 2 ) -
                sineCoefficientsBodyTest.at( 1 ) - sineCoefficientsBodyTest.at( 0 );

        for( unsigned int n = 2; n <= 2; n++ )
        {
            for( unsigned m = 0; m <=2; m++ )
            {
                BOOST_CHECK_SMALL( cosineDegreeTwoBodyDifference( n, m ), 1.0E-23 );
                BOOST_CHECK_SMALL( sineDegreeTwoBodyDifference( n, m ), 1.0E-23 );
            }
        }

        // Check different function types for two-body case (most complex)
        cosineCoefficientsFunctionTestDegreeTwo.push_back( cosineCoefficientsBodyTest.at( 2 ) );
        sineCoefficientsFunctionTestDegreeTwo.push_back( sineCoefficientsBodyTest.at( 2 ) );
        if( functionTest != 0 )
        {
            // Check if the four constant love number interfaces are the same
            Eigen::MatrixXd cosineDegreeTwoBodyDifference =
                    cosineCoefficientsFunctionTestDegreeTwo.at( functionTest ) - cosineCoefficientsFunctionTestDegreeTwo.at( 0 );
            Eigen::MatrixXd sineDegreeTwoBodyDifference =
                    sineCoefficientsFunctionTestDegreeTwo.at( functionTest ) - sineCoefficientsFunctionTestDegreeTwo.at( 0 );

            for( unsigned int n = 2; n <= 2; n++ )
            {
                for( unsigned m = 0; m <=2; m++ )
                {
                    BOOST_CHECK_SMALL( cosineDegreeTwoBodyDifference( n, m ), 1.0E-23 );
                    BOOST_CHECK_SMALL( sineDegreeTwoBodyDifference( n, m ), 1.0E-23 );
                }
            }
        }
    }

    std::vector< Eigen::MatrixXd > cosineCoefficientsDegreeTest;
    std::vector< Eigen::MatrixXd > sineCoefficientsDegreeTest;

    for( int degreeTest = 0; degreeTest < 4; degreeTest++ )
    {
        bodySettings.at( "Earth" )->gravityFieldVariationSettings.clear( );
        std::string deformingBody = "Moon";

        if( degreeTest == 0 )
        {
            bodySettings.at( "Earth" )->gravityFieldVariationSettings.push_back(
                        fixedSingleDegreeLoveNumberGravityFieldVariationSettings(
                            deformingBody, loveNumber, 2 ) );
        }
        else if( degreeTest == 1 )
        {
            bodySettings.at( "Earth" )->gravityFieldVariationSettings.push_back(
                        fixedSingleDegreeLoveNumberGravityFieldVariationSettings(
                            deformingBody, loveNumberDegreeThree, 3 ) );
        }
        else if( degreeTest == 2 )
        {
            std::map< int, double > loveNumberMap;
            loveNumberMap[ 2 ] = 0.25;
            loveNumberMap[ 3 ] = loveNumberDegreeThree;

            bodySettings.at( "Earth" )->gravityFieldVariationSettings.push_back(
                        fixedSingleDegreeLoveNumberGravityFieldVariationSettings(
                            deformingBody, loveNumberMap ) );
        }
        else if( degreeTest == 3 )
        {
            std::map< int, double > loveNumberMap;
            loveNumberMap[ 2 ] = 0.25;
            bodySettings.at( "Earth" )->gravityFieldVariationSettings.push_back(
                        fixedSingleDegreeLoveNumberGravityFieldVariationSettings(
                            deformingBody, loveNumberMap ) );

            loveNumberMap.clear( );
            loveNumberMap[ 3 ] = loveNumberDegreeThree;
            bodySettings.at( "Earth" )->gravityFieldVariationSettings.push_back(
                        fixedSingleDegreeLoveNumberGravityFieldVariationSettings(
                            deformingBody, loveNumberMap ) );

        }

        // Create bodies
        SystemOfBodies bodies = createSystemOfBodies( bodySettings );

        // Update states.
        bodies.at( "Earth" )->setStateFromEphemeris( testTime );
        bodies.at( "Earth" )->setCurrentRotationalStateToLocalFrameFromEphemeris( testTime );
        bodies.at( "Sun" )->setStateFromEphemeris( testTime );
        bodies.at( "Moon" )->setStateFromEphemeris( testTime );

        // Update gravity field
        std::shared_ptr< gravitation::TimeDependentSphericalHarmonicsGravityField > earthGravityField =
                std::dynamic_pointer_cast< gravitation::TimeDependentSphericalHarmonicsGravityField >(
                    bodies.at( "Earth" )->getGravityFieldModel( ) );
        earthGravityField->update( testTime );

        // Retrieve corrections.
        cosineCoefficientsDegreeTest.push_back(
                    earthGravityField->getCosineCoefficients( ) - cosineCoefficients );
        sineCoefficientsDegreeTest.push_back(
                    earthGravityField->getSineCoefficients( ) - sineCoefficients );
    }

    Eigen::MatrixXd firstCosineDifference =
            cosineCoefficientsDegreeTest.at( 0 ) + cosineCoefficientsDegreeTest.at( 1 ) -
            cosineCoefficientsDegreeTest.at( 2 );
    Eigen::MatrixXd secondCosineDifference =
            cosineCoefficientsDegreeTest.at( 0 ) + cosineCoefficientsDegreeTest.at( 1 ) -
            cosineCoefficientsDegreeTest.at( 2 );

    Eigen::MatrixXd firstSineDifference =
            sineCoefficientsDegreeTest.at( 0 ) + sineCoefficientsDegreeTest.at( 1 ) -
            sineCoefficientsDegreeTest.at( 2 );
    Eigen::MatrixXd secondSineDifference =
            sineCoefficientsDegreeTest.at( 0 ) + sineCoefficientsDegreeTest.at( 1 ) -
            sineCoefficientsDegreeTest.at( 2 );

    for( unsigned int n = 2; n <= 3; n++ )
    {
        for( unsigned m = 0; m <=3; m++ )
        {
            BOOST_CHECK_SMALL( firstCosineDifference( n, m ), 1.0E-23 );
            BOOST_CHECK_SMALL( firstSineDifference( n, m ), 1.0E-23 );
            BOOST_CHECK_SMALL( secondCosineDifference( n, m ), 1.0E-23 );
            BOOST_CHECK_SMALL( secondSineDifference( n, m ), 1.0E-23 );
        }
    }

    {
        std::vector< Eigen::MatrixXd > cosineCoefficientsOrderTest;
        std::vector< Eigen::MatrixXd > sineCoefficientsOrderTest;

        for( int orderTest = 0; orderTest < 5; orderTest++ )
        {
            // Clear for current test
            bodySettings.at( "Earth" )->gravityFieldVariationSettings.clear( );


            if( orderTest == 0 )
            {
                bodySettings.at( "Earth" )->gravityFieldVariationSettings.push_back(
                            fixedSingleDegreeLoveNumberGravityFieldVariationSettings(
                                "Moon", loveNumber, 2 ) );
            }
            else if( orderTest < 4 )
            {
                std::vector< double > loveNumbers = { 0.0, 0.0, 0.0 };
                loveNumbers[ orderTest - 1 ] = loveNumber;
                bodySettings.at( "Earth" )->gravityFieldVariationSettings.push_back(
                            orderVariableSingleDegreeLoveNumberGravityFieldVariationSettings(
                                "Moon", loveNumbers, 2 ) );
            }
            else
            {
                std::vector< double > loveNumbers = { 0.0, 0.0, 0.0 };
                for( int order = 0; order < 3; order++ )
                {
                    loveNumbers[ order ] = loveNumber;
                    bodySettings.at( "Earth" )->gravityFieldVariationSettings.push_back(
                                orderVariableSingleDegreeLoveNumberGravityFieldVariationSettings(
                                    "Moon", loveNumbers, 2 ) );
                }
            }

            // Create bodies
            SystemOfBodies bodies = createSystemOfBodies( bodySettings );

            // Update states.
            bodies.at( "Earth" )->setStateFromEphemeris( testTime );
            bodies.at( "Earth" )->setCurrentRotationalStateToLocalFrameFromEphemeris( testTime );
            bodies.at( "Sun" )->setStateFromEphemeris( testTime );
            bodies.at( "Moon" )->setStateFromEphemeris( testTime );

            // Update gravity field
            std::shared_ptr< gravitation::TimeDependentSphericalHarmonicsGravityField > earthGravityField =
                    std::dynamic_pointer_cast< gravitation::TimeDependentSphericalHarmonicsGravityField >(
                        bodies.at( "Earth" )->getGravityFieldModel( ) );
            earthGravityField->update( testTime );

            // Retrieve corrections.
            cosineCoefficientsOrderTest.push_back(
                        earthGravityField->getCosineCoefficients( ) - cosineCoefficients );
            sineCoefficientsOrderTest.push_back(
                        earthGravityField->getSineCoefficients( ) - sineCoefficients );

        }

        Eigen::MatrixXd firstCosineDifference =
                cosineCoefficientsOrderTest.at( 0 ) - cosineCoefficientsOrderTest.at( 1 ) -
                cosineCoefficientsOrderTest.at( 2 ) - cosineCoefficientsOrderTest.at( 3 );
        Eigen::MatrixXd firstSineDifference =
                sineCoefficientsOrderTest.at( 0 ) - sineCoefficientsOrderTest.at( 1 ) -
                sineCoefficientsOrderTest.at( 2 ) - sineCoefficientsOrderTest.at( 3 );

        Eigen::MatrixXd secondCosineDifference =
                cosineCoefficientsOrderTest.at( 4 ) - cosineCoefficientsOrderTest.at( 1 ) -
                cosineCoefficientsOrderTest.at( 2 ) - cosineCoefficientsOrderTest.at( 3 );
        Eigen::MatrixXd secondSineDifference =
                sineCoefficientsOrderTest.at( 4 ) - sineCoefficientsOrderTest.at( 1 ) -
                sineCoefficientsOrderTest.at( 2 ) - sineCoefficientsOrderTest.at( 3 );


        for( unsigned int n = 2; n <= 2; n++ )
        {
            for( unsigned m = 0; m <=2; m++ )
            {
                BOOST_CHECK_SMALL( firstCosineDifference( n, m ), 1.0E-23 );
                BOOST_CHECK_SMALL( firstSineDifference( n, m ), 1.0E-23 );
            }
        }
    }
}




//! Test set up of rotation models.
BOOST_AUTO_TEST_CASE( test_rotationModelSetup )
{

    // Create settings for simple rotation model.
    Eigen::Matrix3d spiceInitialRotationToTargetFrameMatrix;
    spiceInitialRotationToTargetFrameMatrix
            << -0.9548214974296336, 0.2665104385944917, 0.1314841974018291,
            -0.296591573568662, -0.882413772579987, -0.3652114078848295,
            0.01869081416890206, -0.3877088083617987, 0.9215923900425707;
    double venusRotationRate = unit_conversions::convertDegreesToRadians( -1.4813688 ) /
            physical_constants::JULIAN_DAY;
    std::shared_ptr< SimpleRotationModelSettings > simpleRotationSettings =
            std::make_shared< SimpleRotationModelSettings >
            ( "IAU_VENUS", "J2000", Eigen::Quaterniond( spiceInitialRotationToTargetFrameMatrix ),
              1.0E7, venusRotationRate );

    // Create rotation model using setup function
    std::shared_ptr< ephemerides::RotationalEphemeris > approximateEphemeris =
            createRotationModel( simpleRotationSettings, "Earth" );

    // Create rotation model manually.
    ephemerides::SimpleRotationalEphemeris manualApproximateEphemeris(
                Eigen::Quaterniond( spiceInitialRotationToTargetFrameMatrix ),
                venusRotationRate, 1.0E7,
                "J2000", "IAU_VENUS" );

    // Verify equivalence of automatically set up and manual models.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                ( Eigen::Matrix3d( manualApproximateEphemeris.getRotationToBaseFrame(
                                       4.0E7) ) ),
                ( Eigen::Matrix3d( approximateEphemeris->getRotationToBaseFrame(
                                       4.0E7) ) ),
                std::numeric_limits< double >::epsilon( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                ( Eigen::Matrix3d( manualApproximateEphemeris.getRotationToTargetFrame(
                                       4.0E7) ) ),
                ( Eigen::Matrix3d( approximateEphemeris->getRotationToTargetFrame(
                                       4.0E7) ) ),
                std::numeric_limits< double >::epsilon( ) );

}



//! Test set up of synchronous rotation model.
BOOST_AUTO_TEST_CASE( test_synchronousRotationModelSetup )
{

    // Load Spice kernels
    spice_interface::loadStandardSpiceKernels( );

    // Create settings for synchronous rotation model.
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Europa" );
    bodiesToCreate.push_back( "Jupiter" );

    // Create bodies needed in simulation
    BodyListSettings bodySettings =
            getDefaultBodySettings( bodiesToCreate );

    bodySettings.at( "Europa" )->rotationModelSettings = std::make_shared< SynchronousRotationModelSettings >( "Jupiter", "ECLIPJ2000", "IAU_Europa" );
    bodySettings.at( "Europa" )->ephemerisSettings->resetFrameOrigin( "Jupiter" );

    SystemOfBodies bodies = createSystemOfBodies( bodySettings );
    

    std::shared_ptr< SynchronousRotationModelSettings > synchronousRotationSettings
            = std::dynamic_pointer_cast< SynchronousRotationModelSettings >( bodySettings.at( "Europa" )->rotationModelSettings );


    // Create rotation model using setup function
    std::shared_ptr< ephemerides::RotationalEphemeris > rotationalEphemerisFromRotationModelSettings =
            createRotationModel( synchronousRotationSettings, "Europa", bodies );


    // Create synchronous rotation model directly.
    std::function< Eigen::Vector6d( const double, bool ) > relativeStateFunction = createRelativeStateFunction( bodies, "Europa", "Jupiter");

    ephemerides::SynchronousRotationalEphemeris synchronousRotationalEphemeris( relativeStateFunction, "Jupiter", "ECLIPJ2000", "IAU_Europa" );


    // Verify equivalence of automatically set up and manual models.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                ( Eigen::Matrix3d( synchronousRotationalEphemeris.getRotationToBaseFrame(
                                       4.0E7) ) ),
                ( Eigen::Matrix3d( rotationalEphemerisFromRotationModelSettings->getRotationToBaseFrame(
                                       4.0E7) ) ),
                std::numeric_limits< double >::epsilon( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                ( Eigen::Matrix3d( synchronousRotationalEphemeris.getRotationToTargetFrame(
                                       4.0E7) ) ),
                ( Eigen::Matrix3d( rotationalEphemerisFromRotationModelSettings->getRotationToTargetFrame(
                                       4.0E7) ) ),
                std::numeric_limits< double >::epsilon( ) );

}



#if TUDAT_BUILD_WITH_SOFA_INTERFACE
//! Test set up of GCRS<->ITRS rotation model
BOOST_AUTO_TEST_CASE( test_earthRotationModelSetup )
{
    std::shared_ptr< GcrsToItrsRotationModelSettings > rotationSettings =
            std::make_shared< GcrsToItrsRotationModelSettings >( );

    // Create rotation model using setup function
    std::shared_ptr< tudat::ephemerides::RotationalEphemeris > earthRotationModel =
            createRotationModel( rotationSettings, "Earth" );

    double testTime = 5.0E7;

    std::shared_ptr< GcrsToItrsRotationModelSettings > gcrsToItrsRotationSettings2000a =
            std::make_shared< GcrsToItrsRotationModelSettings >( basic_astrodynamics::iau_2000_a );
    std::shared_ptr< GcrsToItrsRotationModelSettings > gcrsToItrsRotationSettings2000b =
            std::make_shared< GcrsToItrsRotationModelSettings >( basic_astrodynamics::iau_2000_b );

    std::shared_ptr< tudat::ephemerides::RotationalEphemeris > earthRotationModel2000a =
            createRotationModel( gcrsToItrsRotationSettings2000a, "Earth" );
    std::shared_ptr< tudat::ephemerides::RotationalEphemeris > earthRotationModel2000b =
            createRotationModel( gcrsToItrsRotationSettings2000b, "Earth" );

    Eigen::Quaterniond rotationMatrix = earthRotationModel->getRotationToBaseFrame( testTime );

    Eigen::Quaterniond rotationMatrix2000a = earthRotationModel2000a->getRotationToBaseFrame( testTime );
    Eigen::Matrix3d matrixDeviation2000a = rotationMatrix2000a.toRotationMatrix( ) - rotationMatrix.toRotationMatrix( );

    Eigen::Quaterniond rotationMatrix2000b = earthRotationModel2000b->getRotationToBaseFrame( testTime );
    Eigen::Matrix3d matrixDeviation2000b = rotationMatrix2000b.toRotationMatrix( ) - rotationMatrix.toRotationMatrix( );

    std::shared_ptr< tudat::ephemerides::GcrsToItrsRotationModel >  defaultEarthModel =
            std::make_shared< ephemerides::GcrsToItrsRotationModel >(
                tudat::earth_orientation::createStandardEarthOrientationCalculator( ) );

    Eigen::Quaterniond defaultRotationMatrix = defaultEarthModel->getRotationToBaseFrame( testTime );
    Eigen::Matrix3d matrixDeviation = defaultRotationMatrix.toRotationMatrix( ) - rotationMatrix.toRotationMatrix( );

    for( unsigned int i = 0; i < 3; i++ )
    {
        for( unsigned int j = 0; j < 3; j++ )
        {
            if( i < 2 && j < 2 )
            {
                BOOST_CHECK_SMALL( std::fabs( matrixDeviation2000a( i, j ) ), 1.0E-13 );
                BOOST_CHECK_SMALL( std::fabs( matrixDeviation2000b( i, j ) ), 1.0E-12 );
                BOOST_CHECK_SMALL( std::fabs( matrixDeviation( i, j ) ), std::numeric_limits< double >::epsilon( ) );
            }
            else
            {
                BOOST_CHECK_SMALL( std::fabs( matrixDeviation2000a( i, j ) ), 1.0E-10 );
                BOOST_CHECK_SMALL( std::fabs( matrixDeviation2000b( i, j ) ), 1.0E-8 );
                BOOST_CHECK_SMALL( std::fabs( matrixDeviation( i, j ) ), std::numeric_limits< double >::epsilon( ) );
            }
        }
    }
}
#endif


//! Test set up of radiation pressure interfacel environment models.
BOOST_AUTO_TEST_CASE( test_radiationPressureInterfaceSetup )
{
    // Load Spice kernels
    spice_interface::loadStandardSpiceKernels( );

    // Define body settings.
    BodyListSettings bodySettings;
    bodySettings.addSettings( getDefaultSingleBodySettings( "Earth", 0.0, 1.0E7 ), "Earth" );
    bodySettings.addSettings( getDefaultSingleBodySettings( "Sun", 0.0, 1.0E7 ), "Sun" );

    // Get settings for vehicle
    double area = 2.34;
    double coefficient = 1.2;
    Eigen::Vector6d initialKeplerElements =
            ( Eigen::Vector6d( ) << 12000.0E3, 0.13, 0.3, 0.0, 0.0, 0.0 ).finished( );
    bodySettings.addSettings( "Vehicle" );
    bodySettings.at( "Vehicle" )->radiationPressureSettings[ "Sun" ] =
            std::make_shared< CannonBallRadiationPressureInterfaceSettings >( "Sun", area, coefficient );
    bodySettings.at( "Vehicle" )->ephemerisSettings =
            std::make_shared< KeplerEphemerisSettings >(
                initialKeplerElements,
                0.0, spice_interface::getBodyGravitationalParameter( "Earth" ), "Earth", "ECLIPJ2000" );

    // Create bodies
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );
    

    BOOST_CHECK_EQUAL( bodies.at( "Vehicle" )->getRadiationPressureInterfaces( ).size( ), 1 );
    BOOST_CHECK_EQUAL( bodies.at( "Vehicle" )->getRadiationPressureInterfaces( ).count( "Sun" ), 1 );

    double testTime = 0.5E7;

    // Update environment to current time.
    bodies.at( "Sun" )->setStateFromEphemeris< double, double >( testTime );
    bodies.at( "Earth" )->setStateFromEphemeris< double, double >( testTime );
    bodies.at( "Vehicle" )->setStateFromEphemeris< double, double >( testTime );

    std::shared_ptr< electromagnetism::RadiationPressureInterface > vehicleRadiationPressureInterface =
            bodies.at( "Vehicle" )->getRadiationPressureInterfaces( ).at( "Sun" );

    vehicleRadiationPressureInterface->updateInterface( testTime );
    double sourceDistance = ( ( bodies.at( "Vehicle" )->getState( ) -  bodies.at( "Sun" )->getState( ) ).
                              segment( 0, 3 ) ).norm( );
    double expectedRadiationPressure = electromagnetism::calculateRadiationPressure(
                defaultRadiatedPowerValues.at( "Sun" ), sourceDistance );

    BOOST_CHECK_CLOSE_FRACTION( expectedRadiationPressure,
                                vehicleRadiationPressureInterface->getCurrentRadiationPressure( ),
                                std::numeric_limits< double >::epsilon( ) );

}



//! Test set up of body shape environment models (see testShapeModels).
BOOST_AUTO_TEST_CASE( test_shapeModelSetup )
{
    // Define test values.
    const Eigen::Vector3d testCartesianPosition( 1917032.190, 6029782.349, -801376.113 );
    const double flattening = 1.0 / 298.257223563;
    const double equatorialRadius = 6378137.0;

    // Test spherical body setup
    {
        std::shared_ptr< BodyShapeSettings > shapeSettings = std::make_shared< SphericalBodyShapeSettings >(
                    equatorialRadius );
        std::shared_ptr< BodyShapeModel > shapeModel = createBodyShapeModel(
                    shapeSettings, "Earth" );


        double calculatedAltitude = shapeModel->getAltitude(
                    testCartesianPosition );
        BOOST_CHECK_CLOSE_FRACTION( calculatedAltitude, testCartesianPosition.norm( ) - equatorialRadius,
                                    std::numeric_limits< double >::epsilon( ) );
        BOOST_CHECK_CLOSE_FRACTION( shapeModel->getAverageRadius( ), equatorialRadius,
                                    std::numeric_limits< double >::epsilon( ) );
    }
    // Test oblate spheroid setup
    {
        std::shared_ptr< BodyShapeSettings > shapeSettings = std::make_shared< OblateSphericalBodyShapeSettings >(
                    equatorialRadius, flattening );
        std::shared_ptr< BodyShapeModel > shapeModel = createBodyShapeModel(
                    shapeSettings, "Earth" );


        double calculatedAltitude = shapeModel->getAltitude(
                    testCartesianPosition );
        double manualAltitude = coordinate_conversions::calculateAltitudeOverOblateSpheroid(
                    testCartesianPosition, equatorialRadius, flattening, 1.0E-4 );
        BOOST_CHECK_CLOSE_FRACTION( calculatedAltitude, manualAltitude,
                                    std::numeric_limits< double >::epsilon( ) );
    }

}



//! Test set up of flight conditions object.
BOOST_AUTO_TEST_CASE( test_flightConditionsSetup )
{
    // Load Spice kernels
    spice_interface::loadStandardSpiceKernels( );

    // Define body settings/
    BodyListSettings bodySettings;
    bodySettings.addSettings( getDefaultSingleBodySettings(
                                  "Earth", 0.0, 1.0E7 ), "Earth" );
    bodySettings.addSettings( "Vehicle" );
    bodySettings.at( "Vehicle" ) ->aerodynamicCoefficientSettings =
            std::make_shared< ConstantAerodynamicCoefficientSettings >(
                1.0, 2.0, 3.0, Eigen::Vector3d::Zero( ),
                ( Eigen::Vector3d( ) << -1.1, 0.1, 2.3 ).finished( ),
                Eigen::Vector3d::Zero( ), 1, 1 );

    // Create bodies
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );
    

    // Define expected aerodynamic angles (see testAerodynamicAngleCalculator)
    double testHeadingAngle = 1.229357188236127;
    double testFlightPathAngle = -0.024894033070522;
    double testLatitude = -0.385027359562548;
    double testLongitude = -1.849449608688977;

    double angleOfAttack = 1.232;
    double angleOfSideslip = -0.00322;
    double bankAngle = 2.323432;

    // Create flight conditions object.
    std::shared_ptr< aerodynamics::FlightConditions > vehicleFlightConditions =
            createAtmosphericFlightConditions( bodies.at( "Vehicle" ), bodies.at( "Earth" ),
                                               "Vehicle", "Earth",
                                               [ & ]( ){ return angleOfAttack; },
    [ & ]( ){ return angleOfSideslip; },
    [ & ]( ){ return bankAngle; } );

    // Set vehicle body-fixed state (see testAerodynamicAngleCalculator)
    Eigen::Vector6d vehicleBodyFixedState =
            ( Eigen::Vector6d( ) << -1656517.23153109, -5790058.28764025, -2440584.88186829,
              6526.30784888051, -2661.34558272018, 2377.09572383163 ).finished( );

    // Set states in environment.
    double testTime = 0.5E7;
    Eigen::Vector6d vehicleInertialState =
            ephemerides::transformStateToFrameFromRotations(
                vehicleBodyFixedState,
                bodies.at( "Earth" )->getRotationalEphemeris( )->getRotationToBaseFrame( testTime ),
                bodies.at( "Earth" )->getRotationalEphemeris( )->getDerivativeOfRotationToBaseFrame( testTime ) );
    bodies.at( "Earth" )->setState( Eigen::Vector6d::Zero( ) );
    bodies.at( "Vehicle" )->setState( vehicleInertialState );
    bodies.at( "Earth" )->setCurrentRotationalStateToLocalFrameFromEphemeris( testTime );

    // Update flight conditions
    vehicleFlightConditions->updateConditions( testTime );

    // Check whether calulcate angles and body-fixed state correspond to expected results
    BOOST_CHECK_SMALL(
                std::fabs( vehicleFlightConditions->getAerodynamicAngleCalculator( )
                           ->getAerodynamicAngle( latitude_angle ) - testLatitude),
                10.0 * std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_SMALL(
                std::fabs( vehicleFlightConditions->getAerodynamicAngleCalculator( )
                           ->getAerodynamicAngle( longitude_angle ) - testLongitude),
                10.0 * std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_SMALL(
                std::fabs( vehicleFlightConditions->getAerodynamicAngleCalculator( )
                           ->getAerodynamicAngle( heading_angle ) - testHeadingAngle),
                10.0 * std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_SMALL(
                std::fabs( vehicleFlightConditions->getAerodynamicAngleCalculator( )
                           ->getAerodynamicAngle( flight_path_angle) - testFlightPathAngle),
                10.0 * std::numeric_limits< double >::epsilon( ) );

    BOOST_CHECK_SMALL(
                std::fabs( vehicleFlightConditions->getAerodynamicAngleCalculator( )
                           ->getAerodynamicAngle( angle_of_attack ) - angleOfAttack ),
                10.0 * std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_SMALL(
                std::fabs( vehicleFlightConditions->getAerodynamicAngleCalculator( )
                           ->getAerodynamicAngle( angle_of_sideslip) - angleOfSideslip),
                10.0 * std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_SMALL(
                std::fabs( vehicleFlightConditions->getAerodynamicAngleCalculator( )
                           ->getAerodynamicAngle( bank_angle ) - bankAngle),
                10.0 * std::numeric_limits< double >::epsilon( ) );

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                vehicleFlightConditions->getCurrentBodyCenteredBodyFixedState( ),vehicleBodyFixedState,
                ( 2.0 * std::numeric_limits< double >::epsilon( ) ) );


}


//! Test set up of solar sailing radiation pressure interface environment models.
BOOST_AUTO_TEST_CASE( test_solarSailingRadiationPressureInterfaceSetup )
{

    // Load Spice kernels
    spice_interface::loadStandardSpiceKernels( );

    // Define body settings.
    BodyListSettings bodySettings;
    bodySettings.addSettings( getDefaultSingleBodySettings( "Earth", 0.0, 1.0E7 ), "Earth" );
    bodySettings.addSettings( getDefaultSingleBodySettings( "Sun", 0.0, 1.0E7 ), "Sun" );

    // Get settings for vehicle
    Eigen::Vector6d initialKeplerElements =
            ( Eigen::Vector6d( ) << 12000.0E3, 0.13, 0.3, 0.0, 0.0, 0.0 ).finished( );
    bodySettings.addSettings( "Vehicle" );
    bodySettings.at( "Vehicle" )->ephemerisSettings = std::make_shared< KeplerEphemerisSettings >(
                initialKeplerElements, 0.0, spice_interface::getBodyGravitationalParameter( "Earth" ), "Earth", "ECLIPJ2000" );


    // Define area subject to radiation pressure [m^2].
    double area = 2.0;

    // Set cone angle of the solar sail [rad].
    double coneAngle = 0.25;

    // Get cone angle of the solar sail [rad].
    std::function< double( const double ) > coneAngleFunction
            = [ = ]( const double ){ return coneAngle; };

    // Define clock angle of the solar sail [rad].
    double clockAngle = 0.2;

    // Get clock angle of the solar sail [rad].
    std::function< double( const double ) > clockAngleFunction
            = [ = ]( const double ){ return clockAngle; };

    // Define front emissivity coefficient of the solar sail [-].
    double frontEmissivityCoefficient = 0.4;

    // Define back emissivity coefficient of the solar sail [-].
    double backEmissivityCoefficient = 0.4;

    // Define front Lambertian coefficient of the solar sail [-].
    double frontLambertianCoefficient = 0.4;

    // Define back Lambertian coefficient of the solar sail [-].
    double backLambertianCoefficient = 0.4;

    // Define reflectivity coefficient of the solar sail [-].
    double reflectivityCoefficient = 0.3;

    // Define specular reflection coefficient of the solar sail [-].
    double specularReflectionCoefficient = 1.0;

    // Define central body.
    std::string centralBody = "Earth";

    // Create radiation pressure interface settings.
    std::shared_ptr< SolarSailRadiationInterfaceSettings > radiationPressureInterfaceSettings =
            std::make_shared< SolarSailRadiationInterfaceSettings >(
                "Sun", area, coneAngleFunction, clockAngleFunction, frontEmissivityCoefficient, backEmissivityCoefficient, frontLambertianCoefficient,
                backLambertianCoefficient, reflectivityCoefficient, specularReflectionCoefficient, std::vector< std::string >( ),
                centralBody );

    bodySettings.at( "Vehicle" )->radiationPressureSettings[ "Sun" ] = radiationPressureInterfaceSettings;

    // Create bodies
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );
    

    BOOST_CHECK_EQUAL( bodies.at( "Vehicle" )->getRadiationPressureInterfaces( ).size( ), 1 );
    BOOST_CHECK_EQUAL( bodies.at( "Vehicle" )->getRadiationPressureInterfaces( ).count( "Sun" ), 1 );

    double testTime = 0.5E7;

    // Update environment to current time.
    bodies.at( "Sun" )->setStateFromEphemeris< double, double >( testTime );
    bodies.at( "Earth" )->setStateFromEphemeris< double, double >( testTime );
    bodies.at( "Vehicle" )->setStateFromEphemeris< double, double >( testTime );


    std::shared_ptr< electromagnetism::RadiationPressureInterface > vehicleRadiationPressureInterface =
            bodies.at( "Vehicle" )->getRadiationPressureInterfaces( ).at( "Sun" );

    // Compute expected radiation pressure.
    vehicleRadiationPressureInterface->updateInterface( testTime );
    double sourceDistance = ( ( bodies.at( "Vehicle" )->getState( ) -  bodies.at( "Sun" )->getState( ) ).
                              segment( 0, 3 ) ).norm( );

    double expectedRadiationPressure = electromagnetism::calculateRadiationPressure(
                defaultRadiatedPowerValues.at( "Sun" ), sourceDistance );

    BOOST_CHECK_CLOSE_FRACTION( expectedRadiationPressure,
                                vehicleRadiationPressureInterface->getCurrentRadiationPressure( ),
                                std::numeric_limits< double >::epsilon( ) );

    std::shared_ptr< electromagnetism::SolarSailingRadiationPressureInterface > vehicleSolarSailingRadiationPressureInterface
            = std::dynamic_pointer_cast< electromagnetism::SolarSailingRadiationPressureInterface > ( vehicleRadiationPressureInterface );

    for ( int i = 0 ; i < 3 ; i++ )
    {
        BOOST_CHECK_SMALL( std::fabs( vehicleSolarSailingRadiationPressureInterface->getCentralBodyVelocity( )( )[ i ] -
                                      bodies.at( "Earth" )->getState()[ i + 3 ] ), 1.0E-15 );

        BOOST_CHECK_SMALL( std::fabs( vehicleSolarSailingRadiationPressureInterface->getCurrentVelocityVector( )[ i ] -
                                      ( bodies.at( "Vehicle" )->getState( ) - bodies.at( "Earth" )->getState( ) ).segment(3,3).normalized()[ i ] ), 1.0E-15 );

        BOOST_CHECK_SMALL( std::fabs( vehicleSolarSailingRadiationPressureInterface->getCurrentSolarVector( )[ i ] -
                                      ( - bodies.at( "Vehicle" )->getState( ) + bodies.at( "Sun" )->getState( ) ).segment(0,3)[ i ] ), 1.0E-15 );

    }

}


BOOST_AUTO_TEST_CASE( test_groundStationCreation )
{
    using namespace unit_conversions;
    using namespace coordinate_conversions;

    const double flattening = 1.0 / 298.257223563;
    const double equatorialRadius = 6378137.0;

    BodyListSettings bodySettings;
    bodySettings.addSettings( "Earth" );
    bodySettings.at( "Earth" )->shapeModelSettings = std::make_shared< OblateSphericalBodyShapeSettings >(
                equatorialRadius, flattening );

    Eigen::Vector3d testCartesianPosition( 1917032.190, 6029782.349, -801376.113 );
    Eigen::Vector3d testGeodeticPosition(
                -63.667,  convertDegreesToRadians( -7.26654999 ), convertDegreesToRadians( 72.36312094 ) );
    Eigen::Vector3d testSphericalPosition = coordinate_conversions::convertCartesianToSpherical(
                testCartesianPosition );
    testSphericalPosition( 1 ) = mathematical_constants::PI / 2.0 - testSphericalPosition( 1 );

    bodySettings.at( "Earth" )->groundStationSettings.push_back(
                std::make_shared< GroundStationSettings >( "Station1", testCartesianPosition, cartesian_position  ) );
    bodySettings.at( "Earth" )->groundStationSettings.push_back(
                std::make_shared< GroundStationSettings >( "Station2", testSphericalPosition, spherical_position ) );
    bodySettings.at( "Earth" )->groundStationSettings.push_back(
                std::make_shared< GroundStationSettings >( "Station3", testGeodeticPosition, geodetic_position ) );

    // Create bodies
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );

    

    BOOST_CHECK_EQUAL( bodies.at( "Earth" )->getGroundStationMap( ).count( "Station1" ), 1 );
    BOOST_CHECK_EQUAL( bodies.at( "Earth" )->getGroundStationMap( ).count( "Station2" ), 1 );
    BOOST_CHECK_EQUAL( bodies.at( "Earth" )->getGroundStationMap( ).count( "Station3" ), 1 );
    BOOST_CHECK_EQUAL( bodies.at( "Earth" )->getGroundStationMap( ).size( ), 3 );

    Eigen::Vector3d testPosition1 = bodies.at( "Earth" )->getGroundStationMap( ).at( "Station1" )->
            getNominalStationState( )->getNominalCartesianPosition( );
    Eigen::Vector3d testPosition2= bodies.at( "Earth" )->getGroundStationMap( ).at( "Station2" )->
            getNominalStationState( )->getNominalCartesianPosition( );
    Eigen::Vector3d testPosition3 = bodies.at( "Earth" )->getGroundStationMap( ).at( "Station3" )->
            getNominalStationState( )->getNominalCartesianPosition( );

    for( unsigned int i = 0; i < 3; i++ )
    {
        BOOST_CHECK_SMALL( std::fabs( testPosition1( i ) - testCartesianPosition( i ) ), 1.0E-9 );
        BOOST_CHECK_SMALL( std::fabs( testPosition2( i ) - testCartesianPosition( i ) ), 1.0E-8 );
        BOOST_CHECK_SMALL( std::fabs( testPosition3( i ) - testCartesianPosition( i ) ), 1.0E-3 );
    }
}



//! Test set up of panelled radiation pressure interface environment models.
BOOST_AUTO_TEST_CASE( test_panelledRadiationPressureInterfaceSetup )
{

    // Load Spice kernels
    spice_interface::loadStandardSpiceKernels( );

    // Define body settings.
    BodyListSettings bodySettings;
    bodySettings.addSettings( getDefaultSingleBodySettings( "Earth", 0.0, 1.0E7 ), "Earth" );
    bodySettings.addSettings( getDefaultSingleBodySettings( "Sun", 0.0, 1.0E7 ), "Sun" );

    // Get settings for vehicle
    Eigen::Vector6d initialKeplerElements =
            ( Eigen::Vector6d( ) << 12000.0E3, 0.13, 0.3, 0.0, 0.0, 0.0 ).finished( );
    bodySettings.addSettings( "Vehicle" );
    bodySettings.at( "Vehicle" )->ephemerisSettings = std::make_shared< KeplerEphemerisSettings >(
                initialKeplerElements, 0.0, spice_interface::getBodyGravitationalParameter( "Earth" ), "Earth", "ECLIPJ2000" );




    // Create radiation pressure properties
    std::vector< double > areas;
    areas.push_back( 4.0 );
    areas.push_back( 6.0 );
    areas.push_back( 2.3 );
    areas.push_back( 2.3 );
    areas.push_back( 5.3 );
    areas.push_back( 2.7 );
    areas.push_back( 4.1 );
    areas.push_back( 2.7 );

    std::vector< double > emissivities;
    emissivities.push_back( 0.1 );
    emissivities.push_back( 0.0 );
    emissivities.push_back( 0.1 );
    emissivities.push_back( 0.1 );
    emissivities.push_back( 0.94 );
    emissivities.push_back( 0.1 );
    emissivities.push_back( 0.94 );
    emissivities.push_back( 0.1 );

    std::vector< double > diffuseReflectionCoefficients;
    diffuseReflectionCoefficients.push_back( 0.46 );
    diffuseReflectionCoefficients.push_back( 0.06 );
    diffuseReflectionCoefficients.push_back( 0.46 );
    diffuseReflectionCoefficients.push_back( 0.46 );
    diffuseReflectionCoefficients.push_back( 0.06 );
    diffuseReflectionCoefficients.push_back( 0.46 );
    diffuseReflectionCoefficients.push_back( 0.06 );
    diffuseReflectionCoefficients.push_back( 0.46 );

    std::vector< Eigen::Vector3d > panelSurfaceNormals;
    panelSurfaceNormals.push_back( Eigen::Vector3d::UnitZ( ) );
    panelSurfaceNormals.push_back( - Eigen::Vector3d::UnitZ( ) );
    panelSurfaceNormals.push_back( Eigen::Vector3d::UnitX( ) );
    panelSurfaceNormals.push_back( - Eigen::Vector3d::UnitX( ) );
    panelSurfaceNormals.push_back( Eigen::Vector3d::UnitY( ) );
    panelSurfaceNormals.push_back( Eigen::Vector3d::UnitY( ) );
    panelSurfaceNormals.push_back( - Eigen::Vector3d::UnitY( ) );
    panelSurfaceNormals.push_back( - Eigen::Vector3d::UnitY( ) );



    std::shared_ptr< PanelledRadiationPressureInterfaceSettings > radiationPressureInterfaceSettings =
            std::make_shared< PanelledRadiationPressureInterfaceSettings >(
                "Sun", emissivities, areas, diffuseReflectionCoefficients, panelSurfaceNormals );

    bodySettings.at( "Vehicle" )->radiationPressureSettings[ "Sun" ] = radiationPressureInterfaceSettings;


    // Create bodies
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );
    

    Eigen::Vector7d unitRotationalState = Eigen::Vector7d::Zero( );
    unitRotationalState.segment( 0, 4 ) = linear_algebra::convertQuaternionToVectorFormat(
                Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) ) );
    bodies.at( "Vehicle" )->setCurrentRotationalStateToLocalFrame( unitRotationalState );

    BOOST_CHECK_EQUAL( bodies.at( "Vehicle" )->getRadiationPressureInterfaces( ).size( ), 1 );
    BOOST_CHECK_EQUAL( bodies.at( "Vehicle" )->getRadiationPressureInterfaces( ).count( "Sun" ), 1 );

    double testTime = 0.5E7;

    // Update environment to current time.
    bodies.at( "Sun" )->setStateFromEphemeris< double, double >( testTime );
    bodies.at( "Earth" )->setStateFromEphemeris< double, double >( testTime );
    bodies.at( "Vehicle" )->setStateFromEphemeris< double, double >( testTime );


    std::shared_ptr< electromagnetism::RadiationPressureInterface > vehicleRadiationPressureInterface =
            bodies.at( "Vehicle" )->getRadiationPressureInterfaces( ).at( "Sun" );

    // Compute expected radiation pressure.
    vehicleRadiationPressureInterface->updateInterface( testTime );
    double sourceDistance = ( ( bodies.at( "Vehicle" )->getState( ) -  bodies.at( "Sun" )->getState( ) ).
                              segment( 0, 3 ) ).norm( );

    double expectedRadiationPressure = electromagnetism::calculateRadiationPressure(
                defaultRadiatedPowerValues.at( "Sun" ), sourceDistance );

    BOOST_CHECK_CLOSE_FRACTION( expectedRadiationPressure,
                                vehicleRadiationPressureInterface->getCurrentRadiationPressure( ),
                                std::numeric_limits< double >::epsilon( ) );

}


BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat

