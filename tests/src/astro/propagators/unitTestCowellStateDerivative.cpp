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

#include <boost/make_shared.hpp>
#include <boost/test/unit_test.hpp>

#include "tudat/math/basic/linearAlgebra.h"
#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/astro/basic_astro/orbitalElementConversions.h"

#include "tudat/interface/spice/spiceInterface.h"
#include "tudat/math/integrators/rungeKuttaCoefficients.h"
#include "tudat/math/interpolators/lagrangeInterpolator.h"
#include "tudat/astro/basic_astro/accelerationModel.h"
#include "tudat/astro/basic_astro/keplerPropagator.h"
#include "tudat/io/basicInputOutput.h"

#include "tudat/astro/basic_astro/orbitalElementConversions.h"
#include "tudat/simulation/environment_setup/body.h"
#include "tudat/astro/propagators/nBodyCowellStateDerivative.h"
#include "tudat/simulation/propagation_setup/dynamicsSimulator.h"
#include "tudat/math/integrators/createNumericalIntegrator.h"
#include "tudat/simulation/environment_setup/createBodies.h"
#include "tudat/simulation/estimation_setup/createNumericalSimulator.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"

namespace tudat
{

namespace unit_tests
{


//Using declarations.
using namespace tudat::ephemerides;
using namespace tudat::interpolators;
using namespace tudat::numerical_integrators;
using namespace tudat::spice_interface;
using namespace tudat::simulation_setup;
using namespace tudat::basic_astrodynamics;
using namespace tudat::orbital_element_conversions;
using namespace tudat::propagators;
using namespace tudat;

BOOST_AUTO_TEST_SUITE( test_cowell_propagator )

//! Test to check whether the frame origin transformations are handled properly.
BOOST_AUTO_TEST_CASE( testCowellPopagatorCentralBodies )
{
    //Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Define bodies in simulation.
    unsigned int totalNumberOfBodies = 4;
    std::vector< std::string > bodyNames;
    bodyNames.resize( totalNumberOfBodies );
    bodyNames[ 0 ] = "Earth";
    bodyNames[ 1 ] = "Sun";
    bodyNames[ 2 ] = "Moon";
    bodyNames[ 3 ] = "Mars";


    // Specify initial time
    double initialEphemerisTime = 1.0E7;
    double finalEphemerisTime = 2.0E7;
    double maximumTimeStep = 3600.0;
    double buffer = 5.0 * maximumTimeStep;

    // Create bodies needed in simulation
    BodyListSettings bodySettings = getDefaultBodySettings(
                bodyNames, initialEphemerisTime - buffer, finalEphemerisTime + buffer, "SSB", "ECLIPJ2000" );
    bodySettings.at( "Mars" )->ephemerisSettings->resetFrameOrigin( "Earth" );
    bodySettings.at( "Earth" )->ephemerisSettings->resetFrameOrigin( "Sun" );
    bodySettings.at( "Moon" )->ephemerisSettings->resetFrameOrigin( "Earth" );

    SystemOfBodies bodies = createSystemOfBodies( bodySettings );


    // Set accelerations between bodies that are to be taken into account (mutual point mass gravity between all bodies).
    SelectedAccelerationMap accelerationMap;
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfEarth;
    accelerationsOfEarth[ "Sun" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
    accelerationsOfEarth[ "Moon" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
    accelerationsOfEarth[ "Mars" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
    accelerationMap[ "Earth" ] = accelerationsOfEarth;

    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfSun;
    accelerationsOfSun[ "Moon" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
    accelerationsOfSun[ "Earth" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
    accelerationsOfSun[ "Mars" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
    accelerationMap[ "Sun" ] = accelerationsOfSun;

    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfMoon;
    accelerationsOfMoon[ "Sun" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
    accelerationsOfMoon[ "Earth" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
    accelerationsOfMoon[ "Mars" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
    accelerationMap[ "Moon" ] = accelerationsOfMoon;

    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfMars;
    accelerationsOfMars[ "Sun" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
    accelerationsOfMars[ "Earth" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
    accelerationsOfMars[ "Moon" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
    accelerationMap[ "Mars" ] = accelerationsOfMars;

    // Define list of bodies to propagate
    std::vector< std::string > bodiesToIntegrate;
    bodiesToIntegrate.push_back( "Earth" );
    bodiesToIntegrate.push_back( "Sun" );
    bodiesToIntegrate.push_back( "Moon" );
    bodiesToIntegrate.push_back( "Mars" );
    unsigned int numberOfNumericalBodies = bodiesToIntegrate.size( );

    // Define numerical integrator settings.
    std::shared_ptr< IntegratorSettings< > > integratorSettings =
            std::make_shared< IntegratorSettings< > >
            ( rungeKutta4, initialEphemerisTime, 200.0 );

    // Define central bodies to use in propagation (all w.r.t SSB).
    std::vector< std::string > centralBodies;
    std::map< std::string, std::string > centralBodyMap;
    centralBodies.resize( numberOfNumericalBodies );
    for( unsigned int i = 0; i < numberOfNumericalBodies; i++ )
    {
        centralBodies[ i ] = "SSB";
        centralBodyMap[ bodiesToIntegrate[ i ] ] = centralBodies[ i ];
    }

    // Get initial state vector as input to integration.
    Eigen::VectorXd systemInitialState = getInitialStatesOfBodies(
                bodiesToIntegrate, centralBodies, bodies, initialEphemerisTime );

    // Create acceleration models and propagation settings.
    AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodies, accelerationMap, centralBodyMap );
    std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings< double > >
            ( centralBodies, accelerationModelMap, bodiesToIntegrate, systemInitialState, finalEphemerisTime );

    // Create simulation object and propagate dynamics.
    SingleArcDynamicsSimulator< > dynamicsSimulator(
                bodies, integratorSettings, propagatorSettings, true, false, false );
    std::map< double, Eigen::VectorXd > solutionSet1 = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );

    // Define new central bodies (hierarchical system)
    centralBodies[ 0 ] = "Sun";
    centralBodies[ 1 ] = "SSB";
    centralBodies[ 2 ] = "Earth";
    centralBodies[ 3 ] = "Sun";
    for( unsigned int i = 0; i < numberOfNumericalBodies; i++ )
    {
        centralBodyMap[ bodiesToIntegrate[ i ] ] = centralBodies[ i ];
    }

    systemInitialState = getInitialStatesOfBodies(
                bodiesToIntegrate, centralBodies, bodies, initialEphemerisTime );

    // Create new acceleration models and propagation settings.
    AccelerationMap accelerationModelMap2 = createAccelerationModelsMap(
                bodies, accelerationMap, centralBodyMap );
    std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings2 =
            std::make_shared< TranslationalStatePropagatorSettings< double > >
            ( centralBodies, accelerationModelMap2, bodiesToIntegrate, systemInitialState, finalEphemerisTime );



    // Create new simulation object and propagate dynamics.
    SingleArcDynamicsSimulator< > dynamicsSimulator2(
                bodies, integratorSettings, propagatorSettings2, true, false, true );
    std::map< double, Eigen::VectorXd > solutionSet2 = dynamicsSimulator2.getEquationsOfMotionNumericalSolution( );

    // Create integration and propagation settings for reverse in time propagation
    std::map< double, Eigen::VectorXd >::iterator solutionSetIterator = (--solutionSet2.end( ) );
    Eigen::VectorXd systemFinalState = solutionSetIterator->second;
    std::shared_ptr< IntegratorSettings< > > integratorSettings2 =
            std::make_shared< IntegratorSettings< > >
            ( rungeKutta4, solutionSetIterator->first, -200.0 );
    std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings3 =
            std::make_shared< TranslationalStatePropagatorSettings< double > >
            ( centralBodies, accelerationModelMap2, bodiesToIntegrate, systemFinalState, initialEphemerisTime );

    // Create new simulation object and propagate dynamics backwards in time.
    SingleArcDynamicsSimulator< > dynamicsSimulator3(
               bodies, integratorSettings2, propagatorSettings3, true, false, false );
    std::map< double, Eigen::VectorXd > solutionSet3 = dynamicsSimulator3.getEquationsOfMotionNumericalSolution( );

    // Create interpolators from three numerical solutions (first one is inertial; second and third are non-inertial)
    std::shared_ptr< LagrangeInterpolator< double, Eigen::VectorXd > > interpolator1 =
            std::make_shared< LagrangeInterpolator< double, Eigen::VectorXd > >( solutionSet1, 8 );
    std::shared_ptr< LagrangeInterpolator< double, Eigen::VectorXd > > interpolator2 =
            std::make_shared< LagrangeInterpolator< double, Eigen::VectorXd > >( solutionSet2, 8 );
    std::shared_ptr< LagrangeInterpolator< double, Eigen::VectorXd > > interpolator3 =
            std::make_shared< LagrangeInterpolator< double, Eigen::VectorXd > >( solutionSet3, 8 );

    // Define step size to be out of sync with integration step size.
    double stepSize = 2001.1 + mathematical_constants::PI;
    double currentTime = initialEphemerisTime + stepSize;

    std::map< double, Eigen::VectorXd > analyticalSolutions;

    // Define maps to retrieve propagated orbits from interpolator.
    Eigen::VectorXd currentInertialSolution = Eigen::VectorXd::Zero( 6 * numberOfNumericalBodies );
    Eigen::VectorXd currentNonInertialSolution = Eigen::VectorXd::Zero( 6 * numberOfNumericalBodies );

    // Define map to put inertial orbit reconstructed from non-inertial orbits.
    Eigen::VectorXd reconstructedInertialSolution = Eigen::VectorXd::Zero( 6 * numberOfNumericalBodies );

    // Define error maps.
    Eigen::VectorXd stateDifference = Eigen::VectorXd::Zero( 6 * numberOfNumericalBodies );

    // Test numerical output against results with SSB as origin for ech body,
    std::shared_ptr< ephemerides::Ephemeris > sunEphemeris = bodies.at( "Sun" )->getEphemeris( );
    std::shared_ptr< ephemerides::Ephemeris > earthEphemeris = bodies.at( "Earth" )->getEphemeris( );
    std::shared_ptr< ephemerides::Ephemeris > marsEphemeris = bodies.at( "Mars" )->getEphemeris( );
    std::shared_ptr< ephemerides::Ephemeris > moonEphemeris = bodies.at( "Moon" )->getEphemeris( );

    std::shared_ptr< LagrangeInterpolator< double, Eigen::VectorXd > > currentInterpolator;

    while( currentTime < finalEphemerisTime - stepSize )
    {
        // Retrieve data from interpolators; transform to inertial frames and compare.
        currentInertialSolution = interpolator1->interpolate( currentTime );


        for( unsigned int k = 0; k < 2; k++ )
        {
            if( k == 0 )
            {
                currentInterpolator = interpolator2;
            }
            else
            {
                currentInterpolator = interpolator3;
            }

            currentNonInertialSolution = currentInterpolator->interpolate( currentTime );

            reconstructedInertialSolution.segment( 0, 6 ) = currentNonInertialSolution.segment( 0, 6 ) +
                    currentNonInertialSolution.segment( 6, 6 );
            reconstructedInertialSolution.segment( 6, 6 ) = currentNonInertialSolution.segment( 6, 6 );
            reconstructedInertialSolution.segment( 12, 6 )= currentNonInertialSolution.segment( 12, 6 ) +
                    reconstructedInertialSolution.segment( 0, 6 );
            reconstructedInertialSolution.segment( 18, 6 ) = currentNonInertialSolution.segment( 18, 6 ) +
                    reconstructedInertialSolution.segment( 6, 6 );

            // Compare states.
            stateDifference = reconstructedInertialSolution - currentInertialSolution;

            for( unsigned j = 0; j < 4; j++ )
            {
                if( j != 2 )
                {
                    BOOST_CHECK_SMALL( std::fabs( stateDifference( 0 + 6 * j ) ), 1.0E-2 );
                    BOOST_CHECK_SMALL( std::fabs( stateDifference( 1 + 6 * j ) ), 1.0E-2 );
                    BOOST_CHECK_SMALL( std::fabs( stateDifference( 2 + 6 * j ) ), 1.0E-2 );

                    BOOST_CHECK_SMALL( std::fabs( stateDifference( 3 + 6 * j ) ), 5.0E-9 );
                    BOOST_CHECK_SMALL( std::fabs( stateDifference( 4 + 6 * j ) ), 5.0E-9 );
                    BOOST_CHECK_SMALL( std::fabs( stateDifference( 5 + 6 * j ) ), 5.0E-9 );
                }
                else
                {
                    BOOST_CHECK_SMALL( std::fabs( stateDifference( 0 + 6 * j ) ), 2.0E-1 );
                    BOOST_CHECK_SMALL( std::fabs( stateDifference( 1 + 6 * j ) ), 2.0E-1 );
                    BOOST_CHECK_SMALL( std::fabs( stateDifference( 2 + 6 * j ) ), 2.0E-1 );

                    BOOST_CHECK_SMALL( std::fabs( stateDifference( 3 + 6 * j ) ), 5.0E-7 );
                    BOOST_CHECK_SMALL( std::fabs( stateDifference( 4 + 6 * j ) ), 5.0E-7 );
                    BOOST_CHECK_SMALL( std::fabs( stateDifference( 5 + 6 * j ) ), 5.0E-7 );
                }
            }
        }



        // Test whether ephemeris objects have been properly reset, i.e. whether all states have been properly transformed to the
        // ephemeris frame.

        // Retrieve data from interpolators; transform to inertial frames and compare.
        currentInertialSolution = interpolator1->interpolate( currentTime );

        reconstructedInertialSolution.segment( 0, 6 ) = earthEphemeris->getCartesianState( currentTime ) +
                sunEphemeris->getCartesianState( currentTime );
        reconstructedInertialSolution.segment( 6, 6 ) = sunEphemeris->getCartesianState( currentTime );
        reconstructedInertialSolution.segment( 12, 6 ) =
                moonEphemeris->getCartesianState( currentTime ) +
                earthEphemeris->getCartesianState( currentTime ) +
                sunEphemeris->getCartesianState( currentTime );
        reconstructedInertialSolution.segment( 18, 6 ) =
                marsEphemeris->getCartesianState( currentTime ) +
                earthEphemeris->getCartesianState( currentTime ) +
                sunEphemeris->getCartesianState( currentTime );

        // Compare states.
        stateDifference = reconstructedInertialSolution - currentInertialSolution;
        for( unsigned j = 0; j < 4; j++ )
        {
            if( j != 2 )
            {
                BOOST_CHECK_SMALL( std::fabs( stateDifference( 0 + 6 * j ) ), 1.0E-2 );
                BOOST_CHECK_SMALL( std::fabs( stateDifference( 1 + 6 * j ) ), 1.0E-2 );
                BOOST_CHECK_SMALL( std::fabs( stateDifference( 2 + 6 * j ) ), 1.0E-2 );

                BOOST_CHECK_SMALL( std::fabs( stateDifference( 3 + 6 * j ) ), 5.0E-9 );
                BOOST_CHECK_SMALL( std::fabs( stateDifference( 4 + 6 * j ) ), 5.0E-9 );
                BOOST_CHECK_SMALL( std::fabs( stateDifference( 5 + 6 * j ) ), 5.0E-9 );
            }
            else
            {
                BOOST_CHECK_SMALL( std::fabs( stateDifference( 0 + 6 * j ) ), 2.0E-1 );
                BOOST_CHECK_SMALL( std::fabs( stateDifference( 1 + 6 * j ) ), 2.0E-1 );
                BOOST_CHECK_SMALL( std::fabs( stateDifference( 2 + 6 * j ) ), 2.0E-1 );

                BOOST_CHECK_SMALL( std::fabs( stateDifference( 3 + 6 * j ) ), 5.0E-7 );
                BOOST_CHECK_SMALL( std::fabs( stateDifference( 4 + 6 * j ) ), 5.0E-7 );
                BOOST_CHECK_SMALL( std::fabs( stateDifference( 5 + 6 * j ) ), 5.0E-7 );
            }
        }

        currentTime += stepSize;
    }

}

//! Test to ensure that a point-mass acceleration on a body produces a Kepler orbit (to within
//! numerical error bounds).
template< typename TimeType, typename StateScalarType >
void testCowellPropagationOfKeplerOrbit( )
{
    //Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Define bodies in simulation.
    std::vector< std::string > bodyNames;
    bodyNames.push_back( "Earth" );
    bodyNames.push_back( "Moon" );

    // Specify initial time
    double initialEphemerisTime = 1.0E7;
    double finalEphemerisTime = 2.0E7;
    double maximumTimeStep = 3600.0;
    double buffer = 5.0 * maximumTimeStep;

    // Create bodies needed in simulation
    BodyListSettings bodySettings =
            getDefaultBodySettings(
                bodyNames, initialEphemerisTime - buffer, finalEphemerisTime + buffer,"SSB", "ECLIPJ2000" );

    if( std::is_same< long double, StateScalarType >::value )
    {
        std::dynamic_pointer_cast< InterpolatedSpiceEphemerisSettings >(
                    bodySettings.at( "Moon" )->ephemerisSettings )->setUseLongDoubleStates( 1 );
    }

    // Change ephemeris settings of Moon and Earth to make test results analysis more transparent.
    std::dynamic_pointer_cast< InterpolatedSpiceEphemerisSettings >( bodySettings.at( "Moon" )->ephemerisSettings )->
            resetFrameOrigin( "Earth" );
    bodySettings.at( "Earth" )->ephemerisSettings = std::make_shared< ConstantEphemerisSettings >(
                Eigen::Vector6d::Zero( ), "SSB", "ECLIPJ2000" );

    SystemOfBodies bodies = createSystemOfBodies( bodySettings );

    // Set accelerations between bodies that are to be taken into account.
    SelectedAccelerationMap accelerationMap;
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfMoon;
    accelerationsOfMoon[ "Earth" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
    accelerationMap[ "Moon" ] = accelerationsOfMoon;

    // Propagate the moon only
    std::vector< std::string > bodiesToIntegrate;
    bodiesToIntegrate.push_back( "Moon" );
    unsigned int numberOfNumericalBodies = bodiesToIntegrate.size( );

    // Define settings for numerical integrator.
    std::shared_ptr< IntegratorSettings< TimeType > > integratorSettings =
            std::make_shared< IntegratorSettings< TimeType > >
            ( rungeKutta4, initialEphemerisTime, 120.0 );

    // Run test where Moon gravity is/is not taken into account.
    for( unsigned testCase = 0; testCase < 2; testCase++ )
    {
        // Get initial kepler elements
        StateScalarType effectiveGravitationalParameter;
        if( testCase == 0 )
        {
            effectiveGravitationalParameter =
                    bodies.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( ) +
                    bodies.at( "Moon" )->getGravityFieldModel( )->getGravitationalParameter( );
        }
        else
        {
            effectiveGravitationalParameter =
                    bodies.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( );
        }

        // Define central bodies for integration.
        std::vector< std::string > centralBodies;
        std::map< std::string, std::string > centralBodyMap;

        if( testCase == 0 )
        {
            effectiveGravitationalParameter =
                    bodies.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( ) +
                    bodies.at( "Moon" )->getGravityFieldModel( )->getGravitationalParameter( );
            centralBodies.push_back( "Earth" );

        }
        else
        {
            effectiveGravitationalParameter =
                    bodies.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( );
            centralBodies.push_back( "SSB" );
        }
        centralBodyMap[ bodiesToIntegrate[ 0 ] ] = centralBodies[ 0 ];


        // Create system initial state.
        Eigen::Matrix< StateScalarType, 6, 1  > systemInitialState =
                Eigen::Matrix< StateScalarType, 6, 1  >( bodiesToIntegrate.size( ) * 6 );
        for( unsigned int i = 0; i < numberOfNumericalBodies ; i++ )
        {
            systemInitialState.segment( i * 6 , 6 ) =
                    spice_interface::getBodyCartesianStateAtEpoch(
                      bodiesToIntegrate[ i ], "Earth", "ECLIPJ2000", "NONE", initialEphemerisTime ).
                    template cast< StateScalarType >( );
        }

        // Create acceleration models and propagation settings.
        AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                    bodies, accelerationMap, bodiesToIntegrate, centralBodies );
        std::shared_ptr< TranslationalStatePropagatorSettings< StateScalarType > > propagatorSettings =
                std::make_shared< TranslationalStatePropagatorSettings< StateScalarType > >
                ( centralBodies, accelerationModelMap, bodiesToIntegrate, systemInitialState, finalEphemerisTime );

        // Create dynamics simulation object.
        SingleArcDynamicsSimulator< StateScalarType, TimeType > dynamicsSimulator(
                    bodies, integratorSettings, propagatorSettings, true, false, true );

        Eigen::Matrix< StateScalarType, 6, 1  > initialKeplerElements =
            orbital_element_conversions::convertCartesianToKeplerianElements< StateScalarType >(
                Eigen::Matrix< StateScalarType, 6, 1  >( systemInitialState ), effectiveGravitationalParameter );


        // Compare numerical state and kepler orbit at each time step.
        std::shared_ptr< Ephemeris > moonEphemeris = bodies.at( "Moon" )->getEphemeris( );
        double currentTime = initialEphemerisTime + buffer;
        while( currentTime < finalEphemerisTime - buffer )
        {

            Eigen::VectorXd stateDifference
                = ( orbital_element_conversions::convertKeplerianToCartesianElements(
                    propagateKeplerOrbit< StateScalarType >( initialKeplerElements, currentTime - initialEphemerisTime,
                                          effectiveGravitationalParameter ),
                    effectiveGravitationalParameter )
                - moonEphemeris->template getTemplatedStateFromEphemeris< StateScalarType >( currentTime ) ).
                    template cast< double >( );

            for( int i = 0; i < 3; i++ )
            {
                BOOST_CHECK_SMALL( stateDifference( i ), 1E-3 );
                BOOST_CHECK_SMALL( stateDifference( i  + 3 ), 2.0E-9 );

            }
            currentTime += 10000.0;
        }
    }
}

BOOST_AUTO_TEST_CASE( testCowellPropagatorKeplerCompare )
{
    testCowellPropagationOfKeplerOrbit< double, double >( );

    #if( TUDAT_BUILD_WITH_EXTENDED_PRECISION_PROPAGATION_TOOLS )
    testCowellPropagationOfKeplerOrbit< double, long double >( );
    testCowellPropagationOfKeplerOrbit< Time, double >( );
    testCowellPropagationOfKeplerOrbit< Time, long double >( );
    #endif

}

BOOST_AUTO_TEST_SUITE_END( )


}

}
