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

#include <string>
#include <thread>

#include <limits>

#include <boost/make_shared.hpp>
#include <boost/test/unit_test.hpp>

#include "tudat/basics/testMacros.h"
#include "tudat/simulation/estimation.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/directTidalTimeLag.h"


namespace tudat
{
namespace unit_tests
{
BOOST_AUTO_TEST_SUITE( test_tidal_propery_estimation )

//Using declarations.
using namespace tudat::observation_models;
using namespace tudat::orbit_determination;
using namespace tudat::estimatable_parameters;
using namespace tudat::interpolators;
using namespace tudat::numerical_integrators;
using namespace tudat::spice_interface;
using namespace tudat::simulation_setup;
using namespace tudat::orbital_element_conversions;
using namespace tudat::ephemerides;
using namespace tudat::propagators;
using namespace tudat::basic_astrodynamics;

//! Unit test to check if tidal time lag parameters are estimated correctly
BOOST_AUTO_TEST_CASE( test_DissipationParameterEstimation )
{
    //Load spice kernels
    spice_interface::loadStandardSpiceKernels( );

    //Define list of bodies
    std::vector< std::string > bodyNames;
    bodyNames.push_back( "Sun" );
    bodyNames.push_back( "Jupiter" );
    bodyNames.push_back( "Io" );
    bodyNames.push_back( "Europa" );

    std::vector< std::string > satelliteNames;
    satelliteNames.push_back( "Io" );
    satelliteNames.push_back( "Europa" );

    // Specify initial time
    double initialEphemerisTime = 0.0;
    double finalEphemerisTime = 0.25 * physical_constants::JULIAN_YEAR;
    double fixedStepSize = 450.0;
    double buffer = 10.0 * fixedStepSize;

    // Get default settings
    BodyListSettings bodySettings =
            getDefaultBodySettings( bodyNames, initialEphemerisTime - buffer, finalEphemerisTime + buffer );
    Eigen::MatrixXd cosineCoefficients = Eigen::MatrixXd::Zero( 3, 3 );

    // Update gravity field settings with sh field, required by tidal models (reference radius)
    cosineCoefficients( 0, 0 ) = 1.0;
    Eigen::MatrixXd sineCoefficients = Eigen::MatrixXd::Zero( 3, 3 );
    for( unsigned int i = 0; i < bodyNames.size( ); i++ )
    {
        if( bodyNames.at( i ) != "Sun" )
        {
            bodySettings.at( bodyNames.at( i ) )->gravityFieldSettings = std::make_shared< SphericalHarmonicsGravityFieldSettings >
                    ( getBodyGravitationalParameter( bodyNames.at( i ) ), getAverageRadius( bodyNames.at( i ) ),
                      cosineCoefficients, sineCoefficients, "IAU_" + bodyNames.at( i ) );
        }
    }

    SystemOfBodies bodies = createSystemOfBodies( bodySettings );

    

    // Test estimation of Jupiter dissipation without and with frequency dependence
    for( unsigned int test = 0; test < 2; test++ )
    {
        // Define accelerations acting on Io and Europa
        double satelliteLoveNumber = 1.0, satelliteTimeLag = 1000.0;
        double jupiterLoveNumber = 1.0, jupiterTimeLag = 100.0;
        SelectedAccelerationMap accelerationMap;
        for( unsigned int i = 0; i < satelliteNames.size( ); i++ )
        {
            std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfSatellite;
            accelerationsOfSatellite[ "Jupiter" ].push_back( std::make_shared< AccelerationSettings >(
                                                                 point_mass_gravity ) );
            accelerationsOfSatellite[ "Sun" ].push_back( std::make_shared< AccelerationSettings >(
                                                             point_mass_gravity ) );
            accelerationsOfSatellite[ "Jupiter" ].push_back( std::make_shared< DirectTidalDissipationAccelerationSettings >(
                                                                 satelliteLoveNumber, satelliteTimeLag, false, false ) );
            accelerationsOfSatellite[ "Jupiter" ].push_back( std::make_shared< DirectTidalDissipationAccelerationSettings >(
                                                                 jupiterLoveNumber, jupiterTimeLag, false, true ) );
            for( unsigned int j = 0; j < satelliteNames.size( ); j++ )
            {
                if( i != j )
                {
                    accelerationsOfSatellite[ satelliteNames.at( j ) ].push_back(
                                std::make_shared< AccelerationSettings >( point_mass_gravity ) );
                }
            }
            accelerationMap[ satelliteNames.at( i ) ] = accelerationsOfSatellite;
        }

        // Set bodies for which initial state is to be estimated and integrated.
        std::vector< std::string > bodiesToEstimate;
        bodiesToEstimate.push_back( "Io" );
        bodiesToEstimate.push_back( "Europa" );
        std::vector< std::string > centralBodies;
        centralBodies.push_back( "Jupiter" );
        centralBodies.push_back( "Jupiter" );
        AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                    bodies, accelerationMap, bodiesToEstimate, centralBodies );


        std::shared_ptr< PropagatorSettings< double > > propagatorSettings =
                std::make_shared< TranslationalStatePropagatorSettings< double > >
                ( centralBodies, accelerationModelMap, bodiesToEstimate,
                  getInitialStatesOfBodies( bodiesToEstimate, centralBodies, bodies, initialEphemerisTime ),
                  finalEphemerisTime );

        // Set parameters that are to be estimated.
        std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames =
                getInitialStateParameterSettings< double >( propagatorSettings, bodies );
        parameterNames.push_back(
                    std::make_shared< tudat::estimatable_parameters::DirectTidalTimeLagEstimatableParameterSettings >(
                        "Io", "" ) );
        parameterNames.push_back(
                    std::make_shared< tudat::estimatable_parameters::DirectTidalTimeLagEstimatableParameterSettings >(
                        "Europa", "" ) );
        if( test == 0 )
        {
            parameterNames.push_back(
                        std::make_shared< tudat::estimatable_parameters::DirectTidalTimeLagEstimatableParameterSettings >(
                            "Jupiter", "" ) );
        }
        else
        {
            parameterNames.push_back(
                        std::make_shared< tudat::estimatable_parameters::DirectTidalTimeLagEstimatableParameterSettings >(
                            "Jupiter", "Io" ) );
            parameterNames.push_back(
                        std::make_shared< tudat::estimatable_parameters::DirectTidalTimeLagEstimatableParameterSettings >(
                            "Jupiter", "Europa" ) );
        }
        std::shared_ptr< tudat::estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate =
                createParametersToEstimate< double >( parameterNames, bodies, propagatorSettings );


        // Define links in simulation.
        std::vector< LinkEnds > linkEnds;
        linkEnds.resize( 2 );
        linkEnds[ 0 ][ observed_body ] = std::make_pair( "Io", "" );
        linkEnds[ 1 ][ observed_body ] = std::make_pair( "Europa", "" );
        std::vector< std::shared_ptr< ObservationModelSettings > > observationSettingsList;
        observationSettingsList.push_back( std::make_shared< ObservationModelSettings >(
                                           position_observable, linkEnds[ 0 ] ) );
        observationSettingsList.push_back( std::make_shared< ObservationModelSettings >(
                                           position_observable, linkEnds[ 1 ] ) );

        // Define integrator and propagator settings.
        std::shared_ptr< IntegratorSettings< > > integratorSettings =
                std::make_shared< RungeKuttaVariableStepSizeSettings< > >
                ( 0.0, fixedStepSize,
                  RungeKuttaCoefficients::rungeKuttaFehlberg78, fixedStepSize, fixedStepSize, 1.0, 1.0 );

        // Create orbit determination object.
        OrbitDeterminationManager< double, double > orbitDeterminationManager =
                OrbitDeterminationManager< double, double >(
                    bodies, parametersToEstimate,
                    observationSettingsList, integratorSettings, propagatorSettings );

        Eigen::VectorXd initialParameterEstimate =
                parametersToEstimate->template getFullParameterValues< double >( );

        // Test if tidal time lag parameters are correctly linked to acceleration models
        std::vector< std::shared_ptr< EstimatableParameter< double > > > dissipationEstimatedParameters =
                parametersToEstimate->getEstimatedDoubleParameters( );
        if( test == 0 )
        {
            BOOST_CHECK_EQUAL( dissipationEstimatedParameters.size( ), 3 );
        }
        else
        {
            BOOST_CHECK_EQUAL( dissipationEstimatedParameters.size( ), 4 );
        }
        for( unsigned int i = 0; i < dissipationEstimatedParameters.size( ); i++ )
        {
            std::shared_ptr< DirectTidalTimeLag > currentTidalTimeLagParameter =
                    std::dynamic_pointer_cast< DirectTidalTimeLag >( dissipationEstimatedParameters.at( i ) );
            BOOST_CHECK_EQUAL( ( currentTidalTimeLagParameter == nullptr ), false );
            int numberOfAccelerationModels = currentTidalTimeLagParameter->getTidalAccelerationModels( ).size( );

            if( ( i == 2 ) && ( test == 0 ) )
            {
                BOOST_CHECK_EQUAL( numberOfAccelerationModels, 2 );
            }
            else
            {
                BOOST_CHECK_EQUAL( numberOfAccelerationModels, 1 );
            }


        }

        // Define observatoion simulation times
        std::vector< double > observationTimes;
        double observationTime = initialEphemerisTime + 10.0 * fixedStepSize;
        while( observationTime < finalEphemerisTime - 10.0 * fixedStepSize  )
        {
            observationTimes.push_back( observationTime );
            observationTime += 6.0 * 3600.0;;
        }


        std::vector< std::shared_ptr< ObservationSimulationSettings< > > > measurementSimulationInput;
        measurementSimulationInput.push_back(
                    std::make_shared< TabulatedObservationSimulationSettings< > >(
                        position_observable, linkEnds[ 0 ], observationTimes, observed_body ) );
        measurementSimulationInput.push_back(
                    std::make_shared< TabulatedObservationSimulationSettings< > >(
                        position_observable, linkEnds[ 1 ], observationTimes, observed_body ) );

        // Simulate observations
        std::shared_ptr< ObservationCollection< > >  observationsAndTimes = simulateObservations< double, double >(
                    measurementSimulationInput, orbitDeterminationManager.getObservationSimulators( ), bodies );


        // Perturb parameter values
        Eigen::VectorXd truthParameters = initialParameterEstimate;
        for( unsigned int i = 0; i < bodiesToEstimate.size( ); i++ )
        {
            initialParameterEstimate[ 0 + 6 * i ] += 1.0E0;
            initialParameterEstimate[ 1 + 6 * i ] += 1.0E0;
            initialParameterEstimate[ 2 + 6 * i ] += 1.0E0;
            initialParameterEstimate[ 3 + 6 * i ] += 1.0E-5;
            initialParameterEstimate[ 4 + 6 * i ] += 1.0E-5;
            initialParameterEstimate[ 5 + 6 * i ] += 1.0E-5;
        }
        for( unsigned int i = 6 * bodiesToEstimate.size( );
             i < static_cast< unsigned int >( initialParameterEstimate.rows( ) ); i++ )
        {
            initialParameterEstimate[ i ] *= 10.0;
        }
        parametersToEstimate->resetParameterValues( initialParameterEstimate );

        // Estimate initial states and tidal parameters
        std::shared_ptr< PodInput< double, double > > podInput =
                std::make_shared< PodInput< double, double > >(
                    observationsAndTimes, ( initialParameterEstimate ).rows( ) );
        std::shared_ptr< PodOutput< double > > podOutput = orbitDeterminationManager.estimateParameters(
                    podInput, std::make_shared< EstimationConvergenceChecker >( 3 ) );

        // Check if parameters are correctly estimated
        Eigen::VectorXd estimatedParametervalues = podOutput->parameterEstimate_;
        for( unsigned int i = 0; i < 3; i++ )
        {
            BOOST_CHECK_SMALL( std::fabs( truthParameters( i ) - podOutput->parameterEstimate_( i ) ), 0.1 );
            BOOST_CHECK_SMALL( std::fabs( truthParameters( i + 6 ) - podOutput->parameterEstimate_( i + 6 ) ), 0.1 );
            BOOST_CHECK_SMALL( std::fabs( truthParameters( i + 3 ) - podOutput->parameterEstimate_( i + 3 ) ), 1.0E-6 );
            BOOST_CHECK_SMALL( std::fabs( truthParameters( i + 9 ) - podOutput->parameterEstimate_( i + 9 ) ), 1.0E-6 );
        }
        BOOST_CHECK_SMALL( std::fabs( truthParameters( 12 ) - podOutput->parameterEstimate_( 12 ) ), 0.1 );
        BOOST_CHECK_SMALL( std::fabs( truthParameters( 13 ) - podOutput->parameterEstimate_( 13 ) ), 10.0 );
        BOOST_CHECK_SMALL( std::fabs( truthParameters( 14 ) - podOutput->parameterEstimate_( 14 ) ), 1.0E-3 );
        if( test == 1 )
        {
            BOOST_CHECK_SMALL( std::fabs( truthParameters( 15 ) - podOutput->parameterEstimate_( 15 ) ), 1.0E-2 );
        }

        std::cout << "Parameter error: " << ( truthParameters -  podOutput->parameterEstimate_ ).transpose( ) << std::endl;
    }
}

//! Function to get tidal deformation model for Earth
std::vector< std::shared_ptr< GravityFieldVariationSettings > > getEarthGravityFieldVariationSettings( )
{
    std::vector< std::shared_ptr< GravityFieldVariationSettings > > gravityFieldVariations;

    std::vector< std::string > deformingBodies;
    deformingBodies.push_back( "Moon" );

    std::map< int, std::vector< std::complex< double > > > loveNumbers;

    std::vector< std::complex< double > > degreeTwoLoveNumbers_;
    degreeTwoLoveNumbers_.push_back( std::complex< double >( 0.29525, -0.00087 ) );
    degreeTwoLoveNumbers_.push_back( std::complex< double >( 0.29525, -0.00087 ) );
    degreeTwoLoveNumbers_.push_back( std::complex< double >( 0.29525, -0.00087 ) );
    std::vector< std::complex< double > > degreeThreeLoveNumbers_;
    degreeThreeLoveNumbers_.push_back( std::complex< double >( 0.093, 0.0 ) );
    degreeThreeLoveNumbers_.push_back( std::complex< double >( 0.093, 0.0 ) );
    degreeThreeLoveNumbers_.push_back( std::complex< double >( 0.093, 0.0 ) );
    degreeThreeLoveNumbers_.push_back( std::complex< double >( 0.093, 0.0 ) );
    loveNumbers[ 2 ] = degreeTwoLoveNumbers_;
    loveNumbers[ 3 ] = degreeThreeLoveNumbers_;


    std::shared_ptr< GravityFieldVariationSettings > moonGravityFieldVariation =
            std::make_shared< BasicSolidBodyGravityFieldVariationSettings >(
                deformingBodies, loveNumbers );
    gravityFieldVariations.push_back( moonGravityFieldVariation );

    deformingBodies[ 0 ] = "Sun";
    std::shared_ptr< GravityFieldVariationSettings > sunSingleGravityFieldVariation =
            std::make_shared< BasicSolidBodyGravityFieldVariationSettings >(
                deformingBodies, loveNumbers );
    gravityFieldVariations.push_back( sunSingleGravityFieldVariation );

    return gravityFieldVariations;
}

//! Unit test to check if real/complex Love numbers at different degrees/orders/forcing frequencies are being correctly estimated.
BOOST_AUTO_TEST_CASE( test_LoveNumberEstimationFromOrbiterData )
{

    //Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Define bodies in simulation
    std::vector< std::string > bodyNames;
    bodyNames.push_back( "Earth" );
    bodyNames.push_back( "Sun" );
    bodyNames.push_back( "Moon" );
    bodyNames.push_back( "Mars" );

    // Specify initial time
    int numberOfDaysOfData = 5;
    double initialEphemerisTime = 1.0E7;
    double finalEphemerisTime = initialEphemerisTime + numberOfDaysOfData * 86400.0;

    // Create bodies needed in simulation
    BodyListSettings bodySettings =
            getDefaultBodySettings( bodyNames, "Earth", "ECLIPJ2000"  );
    bodySettings.at( "Earth" )->gravityFieldVariationSettings = getEarthGravityFieldVariationSettings( );
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );
    bodies.createEmptyBody( "Vehicle" );

    // Set accelerations on Vehicle that are to be taken into account.
    SelectedAccelerationMap accelerationMap;
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfVehicle;
    accelerationsOfVehicle[ "Earth" ].push_back( std::make_shared< SphericalHarmonicAccelerationSettings >( 3, 3 ) );
    accelerationsOfVehicle[ "Sun" ].push_back( std::make_shared< AccelerationSettings >(
                                                   basic_astrodynamics::point_mass_gravity ) );
    accelerationsOfVehicle[ "Moon" ].push_back( std::make_shared< AccelerationSettings >(
                                                    basic_astrodynamics::point_mass_gravity ) );
    accelerationMap[ "Vehicle" ] = accelerationsOfVehicle;

    // Set bodies for which initial state is to be estimated and integrated.
    std::vector< std::string > bodiesToIntegrate;
    std::vector< std::string > centralBodies;
    bodiesToIntegrate.push_back( "Vehicle" );
    centralBodies.push_back( "Earth" );

    // Create acceleration models
    AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodies, accelerationMap, bodiesToIntegrate, centralBodies );

    // Set Keplerian and Cartesian elements for spacecraft.
    Eigen::Vector6d asterixInitialStateInKeplerianElements;
    asterixInitialStateInKeplerianElements( semiMajorAxisIndex ) = 6600.0E3;
    asterixInitialStateInKeplerianElements( eccentricityIndex ) = 0.05;
    asterixInitialStateInKeplerianElements( inclinationIndex ) = unit_conversions::convertDegreesToRadians( 85.3 );
    asterixInitialStateInKeplerianElements( argumentOfPeriapsisIndex )
            = unit_conversions::convertDegreesToRadians( 235.7 );
    asterixInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex )
            = unit_conversions::convertDegreesToRadians( 23.4 );
    asterixInitialStateInKeplerianElements( trueAnomalyIndex ) = unit_conversions::convertDegreesToRadians( 139.87 );
    double earthGravitationalParameter = bodies.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( );
    Eigen::Vector6d systemInitialState = convertKeplerianToCartesianElements(
                asterixInitialStateInKeplerianElements, earthGravitationalParameter );

    // Create propagator settings
    std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings< double > >
            ( centralBodies, accelerationModelMap, bodiesToIntegrate, systemInitialState,
              double( finalEphemerisTime ), cowell );

    // Create integrator settings
    std::shared_ptr< IntegratorSettings< double > > integratorSettings =
            std::make_shared< RungeKuttaVariableStepSizeSettings< double > >
            ( initialEphemerisTime, 60.0,
              RungeKuttaCoefficients::CoefficientSets::rungeKuttaFehlberg78,
              60.0, 60.0, 1.0, 1.0 );

    // Define link ends to use
    LinkEnds linkEnds;
    linkEnds[ observed_body ] = std::make_pair( "Vehicle", "" );
    std::map< ObservableType, std::vector< LinkEnds > > linkEndsPerObservable;
    linkEndsPerObservable[ position_observable ].push_back( linkEnds );

    // Define parameters to estimate (vehicle initial state and various Love number combinations)
    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames =
            getInitialStateParameterSettings< double >( propagatorSettings, bodies );
    parameterNames.push_back( std::make_shared< SingleDegreeVariableTidalLoveNumberEstimatableParameterSettings >(
                                  "Earth", 2, std::vector< int >{ 2 }, "Moon", true ) );
    parameterNames.push_back( std::make_shared< SingleDegreeVariableTidalLoveNumberEstimatableParameterSettings >(
                                  "Earth", 2, std::vector< int >{ 2 }, "Sun", true ) );
    parameterNames.push_back( std::make_shared< SingleDegreeVariableTidalLoveNumberEstimatableParameterSettings >(
                                  "Earth", 2, std::vector< int >{ 0, 1 }, "Moon", false ) );
    parameterNames.push_back( std::make_shared< FullDegreeTidalLoveNumberEstimatableParameterSettings >(
                                  "Earth", 3, "Moon", true ) );

    // Create parameters
    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate =
            createParametersToEstimate( parameterNames, bodies );

    // Define observation settings
    std::vector< std::shared_ptr< ObservationModelSettings > > observationSettingsList;
    observationSettingsList.push_back( std::make_shared< ObservationModelSettings >(
                                       position_observable, linkEnds ) );
    // Create orbit determination object.
    OrbitDeterminationManager< double, double > orbitDeterminationManager =
            OrbitDeterminationManager< double, double >(
                bodies, parametersToEstimate, observationSettingsList,
                integratorSettings, propagatorSettings );

    // define observation simulation times
    std::vector< double > baseTimeList;
    double currentObservationTime = initialEphemerisTime + 1000.0;
    double observationInterval = 60.0;
    while( currentObservationTime < finalEphemerisTime - 1000.0 )
    {
        baseTimeList.push_back( currentObservationTime );
        currentObservationTime += observationInterval;
    }

    std::vector< std::shared_ptr< ObservationSimulationSettings< > > > measurementSimulationInput;
    measurementSimulationInput.push_back(
                std::make_shared< TabulatedObservationSimulationSettings< > >(
                    position_observable, linkEnds, baseTimeList, observed_body ) );

    // Simulate observations
    std::shared_ptr< ObservationCollection< > >  observationsAndTimes = simulateObservations< double, double >(
                measurementSimulationInput, orbitDeterminationManager.getObservationSimulators( ), bodies );

    // Perturb parameter estimate
    Eigen::VectorXd initialParameterEstimate =
            parametersToEstimate->template getFullParameterValues< double >( );
    Eigen::VectorXd truthParameters = initialParameterEstimate;
    Eigen::VectorXd parameterPerturbation =
            Eigen::VectorXd::Zero( truthParameters.rows( ) );
    parameterPerturbation.segment( 0, 3 ) = Eigen::Vector3d::Constant( 1.0 );
    parameterPerturbation.segment( 3, 3 ) = Eigen::Vector3d::Constant( 1.E-3 );

    for( int i = 6; i < parameterPerturbation.rows( ); i++ )
    {
        parameterPerturbation( i ) += 0.001;
    }
    initialParameterEstimate += parameterPerturbation;


    // Define estimation input
    std::shared_ptr< PodInput< double, double  > > podInput =
            std::make_shared< PodInput< double, double > >(
                observationsAndTimes, initialParameterEstimate.rows( ),
                Eigen::MatrixXd::Zero( truthParameters.rows( ), truthParameters.rows( ) ),
                initialParameterEstimate - truthParameters );

    // Perform estimation
    std::shared_ptr< PodOutput< double > > podOutput = orbitDeterminationManager.estimateParameters(
                podInput, std::make_shared< EstimationConvergenceChecker >( 4 ) );

    // Check estimation results
    Eigen::VectorXd estimationError = podOutput->parameterEstimate_ - truthParameters;
    for( int i = 0; i < 3; i++ )
    {
        BOOST_CHECK_SMALL( std::fabs( estimationError( i ) ), 1.0E-4 );
        BOOST_CHECK_SMALL( std::fabs( estimationError( i + 3 ) ), 1.0E-8 );
    }

    for( int i = 0; i < estimationError.rows( ) - 6; i++ )
    {
        BOOST_CHECK_SMALL( std::fabs( estimationError( i + 6 ) ), 5.0E-6 );
    }

    std::cout << "Parameter error"<< ( estimationError ).transpose( ) << std::endl;

}

BOOST_AUTO_TEST_SUITE_END( )

}

}
