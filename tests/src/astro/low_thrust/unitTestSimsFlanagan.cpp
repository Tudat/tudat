﻿/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      Wakker, K. F. (2007), Lecture Notes astro II (Chapter 18), TU Delft course AE4-874,
 *          Delft University of technology, Delft, The Netherlands.
 *
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>
#include <Eigen/Dense>
#include <cmath>
#include <iostream>

#include "tudat/astro/ephemerides/approximatePlanetPositions.h"
#include "tudat/astro/LowThrustTrajectories/ShapeBasedMethods/createBaseFunctionHodographicShaping.h"
#include "tudat/astro/LowThrustTrajectories/ShapeBasedMethods/hodographicShaping.h"
#include "tudat/astro/LowThrustTrajectories/simsFlanagan.h"
#include "tudat/astro/LowThrustTrajectories/simsFlanaganModel.h"
#include "tudat/astro/LowThrustTrajectories/simsFlanaganOptimisationSetup.h"
#include "pagmo/algorithms/de1220.hpp"
#include "tudat/astro/basic_astro/celestialBodyConstants.h"

namespace tudat
{
namespace unit_tests
{

//! Test Sims Flanagan implementation.
BOOST_AUTO_TEST_SUITE( test_Sims_Flanagan )

using namespace tudat;
using namespace low_thrust_trajectories;
using namespace shape_based_methods;
using namespace simulation_setup;
using namespace numerical_integrators;
using namespace propagators;

//! Limit case: if the maximum thrust is set to 0, the Sims Flanagan trajectory should be a keplerian one.
BOOST_AUTO_TEST_CASE( test_Sims_Flanagan_limit_case )
{


    spice_interface::loadStandardSpiceKernels( );

    double maximumThrust = 0.0;
    double specificImpulse = 3000.0;
    double julianDate = 1000.0 * physical_constants::JULIAN_DAY;
    double timeOfFlight = 700.0 * physical_constants::JULIAN_DAY;
    double vehicleInitialMass = 1800.0;
    int numberSegments = 10;

    // Define (constant) specific impulse function.
    std::function< double( const double ) > specificImpulseFunction = [ = ]( const double ){ return specificImpulse; };



    // Define state at departure and arrival.
    ephemerides::EphemerisPointer pointerToDepartureBodyEphemeris = std::make_shared< ephemerides::ApproximatePlanetPositions>(
                ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::earthMoonBarycenter );
    ephemerides::EphemerisPointer pointerToArrivalBodyEphemeris = std::make_shared< ephemerides::ApproximatePlanetPositions >(
                ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::mars );
    Eigen::Vector6d stateAtDeparture = pointerToDepartureBodyEphemeris->getCartesianState( julianDate );
    Eigen::Vector6d stateAtArrival = pointerToArrivalBodyEphemeris->getCartesianState( julianDate + timeOfFlight );

    // Define thrust throttles.
    std::vector< Eigen::Vector3d > throttles;
    for ( int i = 0 ; i < numberSegments ; i++ )
    {
        throttles.push_back( ( Eigen::Vector3d( ) << 0.3, 0.3, 0.3 ).finished( ) );
    }

    // Create Sims Flanagan object
    double centralBodyGravitationalParameter = spice_interface::getBodyGravitationalParameter( "Sun" );
    SimsFlanaganModel simsFlanaganModel = SimsFlanaganModel(
                stateAtDeparture, stateAtArrival, spice_interface::getBodyGravitationalParameter( "Sun" ),
                vehicleInitialMass, maximumThrust, specificImpulseFunction, timeOfFlight, throttles );

    simsFlanaganModel.propagateForwardFromDepartureToMatchPoint( );

    Eigen::Vector6d simsFlanaganForwardPropagation = simsFlanaganModel.getStateAtMatchPointForwardPropagation( );
    Eigen::Vector6d keplerianForwardPropagation = orbital_element_conversions::convertKeplerianToCartesianElements(
                orbital_element_conversions::propagateKeplerOrbit(
                    orbital_element_conversions::convertCartesianToKeplerianElements( stateAtDeparture, centralBodyGravitationalParameter),
                    timeOfFlight / 2.0, centralBodyGravitationalParameter ), centralBodyGravitationalParameter );

    simsFlanaganModel.propagateBackwardFromArrivalToMatchPoint( );

    Eigen::Vector6d simsFlanaganBackwardPropagation = simsFlanaganModel.getStateAtMatchPointBackwardPropagation( );
    Eigen::Vector6d keplerianBackwardPropagation =  orbital_element_conversions::convertKeplerianToCartesianElements(
                orbital_element_conversions::propagateKeplerOrbit(
                    orbital_element_conversions::convertCartesianToKeplerianElements( stateAtArrival, centralBodyGravitationalParameter),
                    - timeOfFlight / 2.0, centralBodyGravitationalParameter ), centralBodyGravitationalParameter );

    // Check consistency between Sims-Flanagan and keplerian propagation when the thrust is set to 0, as a limit case.
    for ( int i = 0 ; i < 3 ; i++ )
    {
        BOOST_CHECK_SMALL( std::fabs( simsFlanaganForwardPropagation[ i ] - keplerianForwardPropagation[ i ] ) / keplerianForwardPropagation.segment( 0,3 ).norm( ), 1.0e-14 );
        BOOST_CHECK_SMALL( std::fabs( simsFlanaganForwardPropagation[ i + 3 ] - keplerianForwardPropagation[ i + 3 ] ) / keplerianForwardPropagation.segment( 3,3 ).norm( ), 1.0e-14 );
        BOOST_CHECK_SMALL( std::fabs( simsFlanaganBackwardPropagation[ i ] - keplerianBackwardPropagation[ i ] ) / keplerianBackwardPropagation.segment( 0,3 ).norm( ), 1.0e-14 );
        BOOST_CHECK_SMALL( std::fabs( simsFlanaganBackwardPropagation[ i + 3 ] - keplerianBackwardPropagation[ i + 3 ] ) / keplerianBackwardPropagation.segment( 3,3 ).norm( ), 1.0e-14 );
    }
}


//! Test if the propagation per segment is equivalent to the total Sims-Flanagan propagation.
BOOST_AUTO_TEST_CASE( test_Sims_Flanagan_implementation )
{


    spice_interface::loadStandardSpiceKernels( );

    double maximumThrust = 0.8;
    double specificImpulse = 3000.0;
    double julianDate = 1000.0 * physical_constants::JULIAN_DAY;
    double timeOfFlight = 700.0 * physical_constants::JULIAN_DAY;
    double vehicleInitialMass = 1800.0;
    int numberSegments = 10;

    // Define (constant) specific impulse function.
    std::function< double( const double ) > specificImpulseFunction = [ = ]( const double ){ return specificImpulse; };

    // Define state at departure and arrival.
    ephemerides::EphemerisPointer pointerToDepartureBodyEphemeris = std::make_shared< ephemerides::ApproximatePlanetPositions>(
                ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::earthMoonBarycenter );
    ephemerides::EphemerisPointer pointerToArrivalBodyEphemeris = std::make_shared< ephemerides::ApproximatePlanetPositions >(
                ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::mars );
    Eigen::Vector6d stateAtDeparture = pointerToDepartureBodyEphemeris->getCartesianState( julianDate );
    Eigen::Vector6d stateAtArrival = pointerToArrivalBodyEphemeris->getCartesianState( julianDate + timeOfFlight );

    // Define thrust throttles.
    std::vector< Eigen::Vector3d > throttles;
    for ( int i = 0 ; i < numberSegments ; i++ )
    {
        throttles.push_back( ( Eigen::Vector3d( ) << 0.3, 0.3, 0.3 ).finished( ) );
    }

    // Create Sims Flanagan object
    double centralBodyGravitationalParameter = spice_interface::getBodyGravitationalParameter( "Sun" );
    SimsFlanaganModel simsFlanaganModel = SimsFlanaganModel(
                stateAtDeparture, stateAtArrival, centralBodyGravitationalParameter,
                vehicleInitialMass, maximumThrust, specificImpulseFunction, timeOfFlight, throttles );
    simsFlanaganModel.propagateForwardFromDepartureToMatchPoint( );
    simsFlanaganModel.propagateBackwardFromArrivalToMatchPoint( );

    int numberSegmentsForwardPropagation = ( numberSegments + 1 ) / 2;
    int numberSegmentsBackwardPropagation = numberSegments / 2;
    double segmentDurationForwardPropagation = timeOfFlight / ( 2.0 * numberSegmentsForwardPropagation );
    double segmentDurationBackwardPropagation = timeOfFlight / ( 2.0 * numberSegmentsBackwardPropagation );
    std::vector< double > timesAtNodes;

    Eigen::Vector6d currentState = stateAtDeparture;
    for ( int i = 0 ; i <= numberSegmentsForwardPropagation ; i++ )
    {
        timesAtNodes.push_back( i * segmentDurationForwardPropagation );
    }
    for ( int i = 0 ; i < numberSegmentsForwardPropagation ; i++ )
    {
        currentState = simsFlanaganModel.propagateInsideForwardSegment( timesAtNodes[ i ], timesAtNodes[ i + 1 ], segmentDurationForwardPropagation,
                currentState );
    }

    // Check consistency between Sims-Flanagan forward propagation to matching point directly, and forward propagation segment by segment.
    for ( int i = 0 ; i < 3 ; i++ )
    {
        BOOST_CHECK_SMALL( std::fabs( currentState[ i ] - simsFlanaganModel.getStateAtMatchPointForwardPropagation( )[ i ] ), 1.0e-15 );
    }

    currentState = stateAtArrival;
    timesAtNodes.clear( );
    for ( int i = 0 ; i <= numberSegmentsBackwardPropagation ; i++ )
    {
        timesAtNodes.push_back( timeOfFlight / 2.0 + i * segmentDurationBackwardPropagation );
    }


    for ( int i = timesAtNodes.size( ) - 1 ; i > 0 ; i-- )
    {
        currentState = simsFlanaganModel.propagateInsideBackwardSegment( timesAtNodes[ i ], timesAtNodes[ i - 1 ], segmentDurationBackwardPropagation,
                currentState );
    }

    // Check consistency between Sims-Flanagan backward propagation to matching point directly, and backward propagation segment by segment.
    for ( int i = 0 ; i < 3 ; i++ )
    {
        BOOST_CHECK_SMALL( std::fabs( currentState[ i ] - simsFlanaganModel.getStateAtMatchPointBackwardPropagation( )[ i ] ), 1.0e-15 );
    }
}



//! Validation of the Sims-Flanagan implementation w.r.t. Pykep.
BOOST_AUTO_TEST_CASE( test_Sims_Flanagan_pykep )
{


    spice_interface::loadStandardSpiceKernels( );

    double maximumThrust = 5.0;
    double specificImpulse = 3000.0;
    double julianDate = 9264.5 * physical_constants::JULIAN_DAY;
    double timeOfFlight = 1000.0 * physical_constants::JULIAN_DAY;
    double vehicleInitialMass = 2800.0;
    int numberSegments = 20;

    int numberSegmentsForwardPropagation = ( numberSegments + 1 ) / 2;
    int numberSegmentsBackwardPropagation = numberSegments / 2;
    double segmentDurationForwardPropagation = timeOfFlight / ( 2.0 * numberSegmentsForwardPropagation );
    double segmentDurationBackwardPropagation = timeOfFlight / ( 2.0 * numberSegmentsBackwardPropagation );

    // Define (constant) specific impulse function.
    std::function< double( const double ) > specificImpulseFunction = [ = ]( const double ){  return specificImpulse;  };

    // Define state at departure and arrival.
    ephemerides::EphemerisPointer pointerToDepartureBodyEphemeris = std::make_shared< ephemerides::ApproximatePlanetPositions>(
                ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::earthMoonBarycenter );
    ephemerides::EphemerisPointer pointerToArrivalBodyEphemeris = std::make_shared< ephemerides::ApproximatePlanetPositions >(
                ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::mars );
    Eigen::Vector6d stateAtDeparture = pointerToDepartureBodyEphemeris->getCartesianState( julianDate );
    Eigen::Vector6d stateAtArrival = pointerToArrivalBodyEphemeris->getCartesianState( julianDate + timeOfFlight );


    // Define thrust throttles.
    std::vector< Eigen::Vector3d > throttles;
    throttles.push_back( ( Eigen::Vector3d( ) << 0.0224237010713910, - 0.0234163377406760, 0.000160068069831085 ).finished( ) );
    throttles.push_back( ( Eigen::Vector3d( ) << 0.0259753967158824, 0.0339093171743966, 2.32915640125360e-05 ).finished( ) );
    throttles.push_back( ( Eigen::Vector3d( ) << - 0.0648979606159804, 0.0651500009913316, - 4.56559330405514e-05 ).finished( ) );
    throttles.push_back( ( Eigen::Vector3d( ) << - 0.162983140398730, - 0.0429364225986646, 0.000465590943761620 ).finished( ) );
    throttles.push_back( ( Eigen::Vector3d( ) << - 0.0967063927336146, - 0.231538912886050, 0.00166394645382429 ).finished( ) );
    throttles.push_back( ( Eigen::Vector3d( ) << 0.144905614712472, - 0.284173666807492, 0.00259915830930070 ).finished( ) );
    throttles.push_back( ( Eigen::Vector3d( ) << 0.345846554068144, - 0.0915917190294576, 0.00192074003319761 ).finished( ) );
    throttles.push_back( ( Eigen::Vector3d( ) << 0.306184858747450, 0.193033585758311, - 0.000509612597294178 ).finished( ) );
    throttles.push_back( ( Eigen::Vector3d( ) << 0.0680496740046574, 0.329263169622182, - 0.00300390820020164 ).finished( ) );
    throttles.push_back( ( Eigen::Vector3d( ) << - 0.160486777991489, 0.243176054265806, - 0.00345409430223536 ).finished( ) );
    throttles.push_back( ( Eigen::Vector3d( ) << - 0.233334021819778, 0.0519439334849342, - 0.00131880814254809 ).finished( ) );
    throttles.push_back( ( Eigen::Vector3d( ) << - 0.162440532964659, - 0.0949724035930260, 0.00187618351107570 ).finished( ) );
    throttles.push_back( ( Eigen::Vector3d( ) << - 0.0445250941509116, - 0.136635969394165, 0.00382254463001364 ).finished( ) );
    throttles.push_back( ( Eigen::Vector3d( ) << 0.0408503328793926, - 0.0993814401580376, 0.00314782099328744 ).finished( ) );
    throttles.push_back( ( Eigen::Vector3d( ) << 0.0695885623491080, - 0.0377735482378564, 0.000253692644808700 ).finished( ) );
    throttles.push_back( ( Eigen::Vector3d( ) << 0.0569202674114514, 0.00947475661200392, - 0.00304163648593162 ).finished( ) );
    throttles.push_back( ( Eigen::Vector3d( ) << 0.0288096343906706, 0.0300320623553084, - 0.00461797222191504 ).finished( ) );
    throttles.push_back( ( Eigen::Vector3d( ) << 0.00432781094877408, 0.0294845775341630, - 0.00328024603096840 ).finished( ) );
    throttles.push_back( ( Eigen::Vector3d( ) << - 0.00919225262002618, 0.0193918543593832, 0.000488305970394490 ).finished( ) );
    throttles.push_back( ( Eigen::Vector3d( ) << - 0.0132219818258071, 0.00882700237255760, 0.00460039760472474 ).finished( ) );

    // Create sims flanagan object
    double centralBodyGravitionalParameter = spice_interface::getBodyGravitationalParameter( "Sun" );
    SimsFlanaganModel simsFlanaganModel = SimsFlanaganModel(
                stateAtDeparture, stateAtArrival, centralBodyGravitionalParameter, vehicleInitialMass,
                maximumThrust, specificImpulseFunction, timeOfFlight, throttles );

    // Compute time at each of the trajectory nodes.
    std::vector< double > timesAtNodes;
    for ( int i = 1 ; i <= numberSegmentsForwardPropagation ; i++ )
    {
        timesAtNodes.push_back( i * segmentDurationForwardPropagation );
    }
    for ( int i = 1 ; i < numberSegmentsBackwardPropagation ; i++ )
    {
        timesAtNodes.push_back( timeOfFlight / 2.0 + i * segmentDurationBackwardPropagation );
    }

    // Retrieve trajectory history
    std::map< double, Eigen::Vector6d > simsFlanaganTrajectory;
    simsFlanaganModel.propagateTrajectory( timesAtNodes, simsFlanaganTrajectory );

    // Expected value for the spacecraft mass at the Sims-Flanagan nodes (from Pykep).
    std::vector< double > expectedMassesAtNodes =
    { 2776.2969549424233, 2745.1123152911364, 2678.420656242383, 2557.4910155457815, 2379.739411828715,
      2156.689164020016, 1909.379054895334, 1661.2975125972175, 1431.9008589844336, 1233.184537707585,
      1069.592792583044, 939.9855510907765, 840.1487300900216, 764.8194149185431, 708.840224065782,
      667.6603951232578, 637.615332355369, 615.9793835860186, 600.4192020709157 };

    // Expected position vector at the Sims-Flanagan nodes (from Pykep).
    std::vector< Eigen::Vector3d > expectedPositionsAtNodes = {
        ( Eigen::Vector3d( ) << 29361822142.48377, -149720536564.364, 12625513.890192498 ).finished( ),
        ( Eigen::Vector3d( ) << 131244674959.76816, - 79709935671.14459, 11081216.394638035 ).finished( ),
        ( Eigen::Vector3d( ) << 147863359688.51364, 43601531813.244675, 1997871.8096662136 ).finished( ),
        ( Eigen::Vector3d( ) << 66395452089.14818, 139648419550.60236, - 1242449.2438607286 ).finished( ),
        ( Eigen::Vector3d( ) << - 61547245597.21583, 143052944943.39178, 33724797.08655864 ).finished( ),
        ( Eigen::Vector3d( ) << - 151138587052.9555, 48840796343.77897, 127386851.75841041 ).finished( ),
        ( Eigen::Vector3d( ) << - 143593965024.7323, - 82091360357.99661, 238529851.87103277 ).finished( ),
        ( Eigen::Vector3d( ) << - 46654714141.73993, - 169642932922.45218, 256360559.2892782 ).finished( ),
        ( Eigen::Vector3d( ) << 81506619982.13998, - 170770379742.6725, 67200021.92568266 ).finished( ),
        ( Eigen::Vector3d( ) << 179238568921.6867, - 94362440826.78069, -343395985.8438245 ).finished( ),
        ( Eigen::Vector3d( ) << 104912577427.10306, 64675894638.26674, 96289313.59812593 ).finished( ),
        ( Eigen::Vector3d( ) << 8320027427.249863, 161615506120.01874, 2679195136.6040783 ).finished( ),
        ( Eigen::Vector3d( ) << - 110027918948.2429, 156486309546.98526, 4095176698.2014832 ).finished( ),
        ( Eigen::Vector3d( ) << - 192071366858.78598, 83352031041.89726, 4452339167.99959 ).finished( ),
        ( Eigen::Vector3d( ) << - 217057910085.4913, - 20476971268.45745, 3812266501.3530645 ).finished( ),
        ( Eigen::Vector3d( ) << - 183120141902.67807, - 121439198146.61942, 2090555690.610314 ).finished( ),
        ( Eigen::Vector3d( ) << - 101541631262.10852, - 191835972941.9502, - 626986013.5014207 ).finished( ),
        ( Eigen::Vector3d( ) << 6218615400.38308, - 212409393619.73843, - 3746237914.8793163 ).finished( ),
        ( Eigen::Vector3d( ) << 111737198408.02193, - 175909675868.96454, - 6114703221.327605 ).finished( ) };

    // Expected velocity vector at the Sims-Flanagan nodes (from Pykep).
    std::vector< Eigen::Vector3d > expectedVelocitiesAtNodes = {
        ( Eigen::Vector3d( ) << 28918.33762542703, 5352.439387790828, 0.5727081570311348 ).finished( ),
        ( Eigen::Vector3d( ) << 15544.49526818364, 25115.421839016366, - 1.2401015025921627 ).finished( ),
        ( Eigen::Vector3d( ) << - 8253.308645516736, 28698.50117509112, - 2.722488814689109 ).finished( ),
        ( Eigen::Vector3d( ) << - 27307.260250963365, 13264.732717716839, 1.3072527713014581 ).finished( ),
        ( Eigen::Vector3d( ) << - 28618.73108014237, - 11708.108000325763, 13.984002704168226 ).finished( ),
        ( Eigen::Vector3d( ) << - 10786.265956564694, - 29432.077637213584, 27.158205556707973 ).finished( ),
        ( Eigen::Vector3d( ) << 13699.144504546844, - 28259.656201215643, 22.265350825225998 ).finished( ),
        ( Eigen::Vector3d( ) << 28946.943030427235, - 11031.229625788088, - 13.45225763389609 ).finished( ),
        ( Eigen::Vector3d( ) << 28417.047752970044, 9928.242838024657, - 70.52296877504925 ).finished( ),
        ( Eigen::Vector3d( ) << 15927.019447896218, 24117.469181704524, - 114.65105018889268 ).finished( ),
        ( Eigen::Vector3d( ) << - 11105.479899668795, 35033.4179444072, 682.7863864659137 ).finished( ),
        ( Eigen::Vector3d( ) << - 28485.627413538656, 9737.851568512157, 461.57152775207305 ).finished( ),
        ( Eigen::Vector3d( ) << - 24659.744784630533, - 10599.769053725684, 198.8135857660587 ).finished( ),
        ( Eigen::Vector3d( ) << - 12835.759652252656, - 21998.298073867958, - 23.971860112821446 ).finished( ),
        ( Eigen::Vector3d( ) << 1268.8403238331266, - 25018.88642908034, - 261.2095191378802 ).finished( ),
        ( Eigen::Vector3d( ) << 14092.63277287905, - 20819.008816988713, - 519.5950098529499 ).finished( ),
        ( Eigen::Vector3d( ) << 22975.29919471634, - 11071.656124219113, - 713.840718279506 ).finished( ),
        ( Eigen::Vector3d( ) << 25900.684080375395, 1871.4112174067868, - 698.7422854478164 ).finished( ),
        ( Eigen::Vector3d( ) << 21828.356239730165, 14767.931720355062, - 370.40551754476553 ).finished( ) };


    // Check consistency w.r.t. pykep for forward propagation.
    for ( int i = 0 ; i < numberSegmentsForwardPropagation ; i++ )
    {
        Eigen::Vector3d calculatedPosition = simsFlanaganTrajectory[ timesAtNodes[ i ] ].segment( 0, 3 );
        Eigen::Vector3d calculatedVelocity = simsFlanaganTrajectory[ timesAtNodes[ i ] ].segment( 3, 3 );
        Eigen::Vector3d expectedPosition = expectedPositionsAtNodes[ i ];
        Eigen::Vector3d expectedVelocity = expectedVelocitiesAtNodes[ i ];


        for ( int i = 0 ; i < 3 ; i++ )
        {
            BOOST_CHECK_SMALL( std::fabs( calculatedPosition[ i ] - expectedPosition[ i ] ) / expectedPosition.norm( ), 1.0e-10 );
            BOOST_CHECK_SMALL( std::fabs( calculatedVelocity[ i ] - expectedVelocity[ i ] ) / expectedVelocity.norm( ), 1.0e-10 );
        }
        BOOST_CHECK_SMALL( std::fabs( std::fabs( simsFlanaganModel.getMassAtSegment( i + 1 ) - expectedMassesAtNodes[ i ] ) ), 1.0e-15 );
    }

    // Check consistency w.r.t. pykep for backward propagation.
    for ( int i = 1 ; i < numberSegmentsBackwardPropagation ; i++ )
    {
        Eigen::Vector3d calculatedPosition = simsFlanaganTrajectory[ timesAtNodes[ i - 1 + numberSegmentsForwardPropagation ] ].segment( 0, 3 );
        Eigen::Vector3d calculatedVelocity = simsFlanaganTrajectory[ timesAtNodes[ i - 1 + numberSegmentsForwardPropagation ] ].segment( 3, 3 );
        Eigen::Vector3d expectedPosition = expectedPositionsAtNodes[ i + numberSegmentsForwardPropagation - 1 ];
        Eigen::Vector3d expectedVelocity = expectedVelocitiesAtNodes[ i + numberSegmentsForwardPropagation - 1 ];

        for ( int i = 0 ; i < 3 ; i++ )
        {
            BOOST_CHECK_SMALL( std::fabs( calculatedPosition[ i ] - expectedPosition[ i ] ) / expectedPosition.norm( ), 1.0e-10 );
            BOOST_CHECK_SMALL( std::fabs( calculatedVelocity[ i ] - expectedVelocity[ i ] ) / expectedVelocity.norm( ), 1.0e-10 );
        }
        BOOST_CHECK_SMALL( std::fabs( simsFlanaganModel.getMassAtSegment( i + numberSegmentsForwardPropagation ) -
                                      expectedMassesAtNodes[ i - 1 + numberSegmentsForwardPropagation ] ), 1.0e-15 );
    }

}

SystemOfBodies getTestBodyMap( )
{

    // Create central, departure and arrival bodies.
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Sun" );

    std::string frameOrigin = "SSB";
    std::string frameOrientation = "ECLIPJ2000";

    BodyListSettings bodySettings =
            getDefaultBodySettings( bodiesToCreate, frameOrigin, frameOrientation );

    // Define central body ephemeris settings.
    bodySettings.at( "Sun" )->ephemerisSettings = std::make_shared< ConstantEphemerisSettings >(
                ( Eigen::Vector6d( ) << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ).finished( ), frameOrigin, frameOrientation );


    // Create system of bodies.
    SystemOfBodies bodies = createBodies( bodySettings );

    bodies.addNewBody( "Vehicle" );
    bodies.at( "Vehicle" )->setEphemeris( std::make_shared< ephemerides::TabulatedCartesianEphemeris< > >(
                                               std::shared_ptr< interpolators::OneDimensionalInterpolator
                                               < double, Eigen::Vector6d > >( ), frameOrigin, frameOrientation ) );
    bodies.at( "Vehicle" )->setSuppressDependentOrientationCalculatorWarning( true );


    return bodies;
}

//! Test Sims-Flanagan implementation by comparing it with trajectory obtained with successive impulsive shots applied at times corresponding to
//! half of each of the Sims-Flanagan segments.
BOOST_AUTO_TEST_CASE( test_Sims_Flanagan_impulsive_shots )
{
    spice_interface::loadStandardSpiceKernels( );

    double maximumThrust = 0.8;
    double specificImpulse = 3000.0;
    double julianDate = 1000.0 * physical_constants::JULIAN_DAY;
    double timeOfFlight = 700.0 * physical_constants::JULIAN_DAY;
    double vehicleInitialMass = 1800.0;
    int numberSegments = 10;
    int numberSegmentsForwardPropagation = ( numberSegments + 1 ) / 2;
    int numberSegmentsBackwardPropagation = numberSegments / 2;
    std::function< double( const double ) > specificImpulseFunction = [ = ]( const double ){ return specificImpulse; };

    // Compute segment duration for the forward and backward propagations.
    int segmentDurationForwardPropagation = timeOfFlight / ( 2.0 * numberSegmentsForwardPropagation );
    int segmentDurationBackwardPropagation = timeOfFlight / ( 2.0 * numberSegmentsBackwardPropagation );

    // Define state at departure and arrival.
    ephemerides::EphemerisPointer pointerToDepartureBodyEphemeris = std::make_shared< ephemerides::ApproximatePlanetPositions>(
                ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::earthMoonBarycenter );
    ephemerides::EphemerisPointer pointerToArrivalBodyEphemeris = std::make_shared< ephemerides::ApproximatePlanetPositions >(
                ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::mars );
    Eigen::Vector6d stateAtDeparture = pointerToDepartureBodyEphemeris->getCartesianState( julianDate );
    Eigen::Vector6d stateAtArrival = pointerToArrivalBodyEphemeris->getCartesianState( julianDate + timeOfFlight );

    //! Create Sims-Flanagan trajectory
    std::map< double, Eigen::Vector6d > simsFlanaganTrajectory;
    std::vector< Eigen::Vector3d > bestThrottles;
    double massAtTimeOfFlight;
    Eigen::Vector6d stateAtHalfTimeOfFlightBackwardPropagation;
    {
        // Define optimisation algorithm.
        algorithm optimisationAlgorithm{ pagmo::de1220() };
        std::shared_ptr< OptimisationSettings > optimisationSettings =
                std::make_shared< OptimisationSettings >( optimisationAlgorithm, 1, 10 );

        // Create sims-flanagan object
        double centralBodyGravitionalParameter = spice_interface::getBodyGravitationalParameter( "Sun" );
        SimsFlanagan simsFlanagan = SimsFlanagan(
                    stateAtDeparture, stateAtArrival, centralBodyGravitionalParameter, vehicleInitialMass, maximumThrust,
                    specificImpulseFunction, numberSegments, timeOfFlight, optimisationSettings );

        // Retrieve trajectory matched at mid-point
        std::vector< double > fitnessVector = simsFlanagan.getBestIndividualFitness( );
        std::vector< double > bestIndividual = simsFlanagan.getBestIndividual( );
        for ( int i = 0 ; i < numberSegments ; i++ )
        {
            bestThrottles.push_back( ( Eigen::Vector3d( ) << bestIndividual[ i * 3 ], bestIndividual[ i * 3 + 1 ],
                    bestIndividual[ i * 3 + 2 ] ).finished( ) );
        }

        // Compute state and mass at half of the time of flight.
        simsFlanagan.getOptimalSimsFlanaganModel( )->propagateBackwardFromArrivalToMatchPoint( );
        stateAtHalfTimeOfFlightBackwardPropagation =
                simsFlanagan.getOptimalSimsFlanaganModel( )->getStateAtMatchPointBackwardPropagation( );
        massAtTimeOfFlight = simsFlanagan.getOptimalSimsFlanaganModel( )->getMassAtSegment( numberSegments );


        // Compute time at each node of the Sims-Flanagan trajectory.
        std::vector< double > timesAtNodes;
        for ( int i = 1 ; i <= numberSegmentsForwardPropagation ; i++ )
        {
            timesAtNodes.push_back( i * segmentDurationForwardPropagation );
        }
        for ( int i = 1 ; i <= numberSegmentsBackwardPropagation ; i++ )
        {
            timesAtNodes.push_back( timeOfFlight / 2.0 + i * segmentDurationBackwardPropagation );
        }

        // Retrieve the spacecraft trajectory at the Sims-Flanagan nodes.
        simsFlanagan.getTrajectory( timesAtNodes, simsFlanaganTrajectory );
    }


    //! DEFINE SETTINGS FOR FULL PROPAGATION

    // Compute times at half of each segment.
    std::vector< double > thrustMidTimes;
    for ( int i = 0 ; i < numberSegmentsForwardPropagation ; i++ )
    {
        thrustMidTimes.push_back( segmentDurationForwardPropagation / 2.0 + i * segmentDurationForwardPropagation );
    }
    for ( int i = 0 ; i < numberSegmentsBackwardPropagation ; i++ )
    {
        thrustMidTimes.push_back( segmentDurationBackwardPropagation / 2.0 + timeOfFlight / 2.0 + i * segmentDurationBackwardPropagation );
    }

    // Set manuever properties
    double totalManeuverTime = 9000.0;
    double maneuverRiseTime = 1500.0;
    double currentMass = vehicleInitialMass;
    std::vector< Eigen::Vector3d > deltaVs;

    // Compute maneuver deltaVs for the forward propagation subleg.
    for ( int i = 0 ; i < numberSegmentsForwardPropagation ; i++ )
    {
        Eigen::Vector3d currentDeltaVvector = maximumThrust * bestThrottles[ i ] * segmentDurationForwardPropagation / currentMass;
        deltaVs.push_back( currentDeltaVvector );
        currentMass *= std::exp( - currentDeltaVvector.norm() /
                                 ( specificImpulseFunction( 0.0 ) * physical_constants::SEA_LEVEL_GRAVITATIONAL_ACCELERATION ) );
    }
    // Compute maneuver deltaVs for the backward propagation subleg.
    for ( int i = 0 ; i < numberSegmentsBackwardPropagation ; i++ )
    {
        Eigen::Vector3d currentDeltaVvector = maximumThrust * bestThrottles[ i + numberSegmentsForwardPropagation ] *
                segmentDurationBackwardPropagation / currentMass;
        deltaVs.push_back( currentDeltaVvector );
        currentMass *= std::exp( - currentDeltaVvector.norm() /
                                 ( specificImpulseFunction( 0.0 ) * physical_constants::SEA_LEVEL_GRAVITATIONAL_ACCELERATION ) );
    }

    // Create list of bodies
    SystemOfBodies bodies = getTestBodyMap( );
    bodies.at( "Vehicle" )->setConstantBodyMass( vehicleInitialMass );

    // Defined acceleration settings
    std::string bodyToPropagate = "Vehicle";
    std::string centralBody = "Sun";
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > bodyToPropagateAccelerations;
    bodyToPropagateAccelerations[ centralBody ].push_back(
                std::make_shared< AccelerationSettings >(
                    basic_astrodynamics::central_gravity ) );
    bodyToPropagateAccelerations[ bodyToPropagate ].push_back(
                std::make_shared< MomentumWheelDesaturationAccelerationSettings >(
                    thrustMidTimes, deltaVs, totalManeuverTime, maneuverRiseTime ) );
    SelectedAccelerationMap accelerationMap;
    accelerationMap[ bodyToPropagate ] = bodyToPropagateAccelerations;


    // Create the acceleration map.
    basic_astrodynamics::AccelerationMap accelerationModelMap =
            createAccelerationModelsMap( bodies, accelerationMap,
                                         std::vector< std::string >{ bodyToPropagate },
                                         std::vector< std::string >{ centralBody } );

    // Propagate orbit forward and backward
    {
        // Define integrator settings.
        double stepSize = ( timeOfFlight ) / static_cast< double >( 250000 );
        std::shared_ptr< IntegratorSettings< double > > integratorSettings =
                std::make_shared< IntegratorSettings< double > > ( rungeKutta4, 0.0, stepSize );

        // Create termination conditions settings.
        std::shared_ptr< PropagationTerminationSettings > terminationSettings
                = std::make_shared< PropagationTimeTerminationSettings >( timeOfFlight / 2.0, true );

        //  Create propagator settings.
        std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettingsImpulsiveDeltaV =
                std::make_shared< TranslationalStatePropagatorSettings< double > >(
                    std::vector< std::string >{ centralBody }, accelerationModelMap, std::vector< std::string >{ bodyToPropagate },
                    stateAtDeparture, terminationSettings, cowell );


        Eigen::Vector6d currentState = stateAtDeparture;
        for ( int i = 0 ; i < numberSegmentsForwardPropagation ; i++ )
        {
            // Define settings for current leg
            integratorSettings->initialTime_ = i * segmentDurationForwardPropagation;
            terminationSettings = std::make_shared< PropagationTimeTerminationSettings >(
                        ( i + 1 ) * segmentDurationForwardPropagation, true );
            propagatorSettingsImpulsiveDeltaV->resetTerminationSettings( terminationSettings );
            propagatorSettingsImpulsiveDeltaV->resetInitialStates( currentState );

            // Propagate current leg
            SingleArcDynamicsSimulator< > dynamicsSimulator( bodies, integratorSettings, propagatorSettingsImpulsiveDeltaV );
            currentState = dynamicsSimulator.getEquationsOfMotionNumericalSolution().rbegin( )->second;

            // Retrieve Sims-Flanagan and numerical states
            Eigen::Vector6d currentStateSimsFlanagan = simsFlanaganTrajectory[ ( i + 1 ) * segmentDurationForwardPropagation ];
            Eigen::Vector6d currentStateImpulsiveShots = currentState;

            // Check consistency between Sims-Flanagan trajectory and pseudo-impulsive shots trajectory.
            for ( int i = 0 ; i < 3 ; i++ )
            {
                BOOST_CHECK_SMALL( std::fabs( currentStateSimsFlanagan[ i ] - currentStateImpulsiveShots[ i ] ) /
                                   currentStateImpulsiveShots.segment( 0,3 ).norm( ), 1.0e-6 );
                BOOST_CHECK_SMALL( std::fabs( currentStateSimsFlanagan[ i + 3 ] - currentStateImpulsiveShots[ i + 3 ] ) /
                        currentStateImpulsiveShots.segment( 3,3 ).norm( ), 1.0e-6 );
            }

        }

        // Reset integrator settings for backwards propagation
        integratorSettings->initialTimeStep_ = - std::fabs( integratorSettings->initialTimeStep_ );
        integratorSettings->initialTime_ = timeOfFlight;

        // Reset termination conditions settings for backwards propagation.
        terminationSettings = std::make_shared< PropagationTimeTerminationSettings >( timeOfFlight / 2.0, true );

        // Reset propagator settings for backwards propagation
        propagatorSettingsImpulsiveDeltaV = std::make_shared< TranslationalStatePropagatorSettings< double > >(
                    std::vector< std::string >{ centralBody }, accelerationModelMap, std::vector< std::string >{ bodyToPropagate },
                    stateAtArrival, terminationSettings, cowell );

        // Reset body mass for backwards propagatiom
        bodies.at( bodyToPropagate )->setConstantBodyMass( massAtTimeOfFlight );

        currentState = stateAtArrival;
        for ( int i = 0 ; i < numberSegmentsBackwardPropagation ; i++ )
        {
            // Define settings for current leg
            integratorSettings->initialTime_ = timeOfFlight - i * segmentDurationBackwardPropagation;
            terminationSettings = std::make_shared< PropagationTimeTerminationSettings >( timeOfFlight - ( i + 1 ) * segmentDurationBackwardPropagation, true );
            propagatorSettingsImpulsiveDeltaV->resetTerminationSettings( terminationSettings );
            propagatorSettingsImpulsiveDeltaV->resetInitialStates( currentState );

            // Propagate current leg
            SingleArcDynamicsSimulator< > dynamicsSimulator( bodies, integratorSettings, propagatorSettingsImpulsiveDeltaV );
            currentState = dynamicsSimulator.getEquationsOfMotionNumericalSolution().begin()->second;

            // Retrieve Sims-Flanagan and numerical states
            Eigen::Vector6d currentStateImpulsiveShots = dynamicsSimulator.getEquationsOfMotionNumericalSolution().begin()->second;
            Eigen::Vector6d currentStateSimsFlanagan;
            if ( i < numberSegmentsBackwardPropagation - 1 )
            {
                currentStateSimsFlanagan = simsFlanaganTrajectory[ timeOfFlight - ( i + 1 ) * segmentDurationBackwardPropagation ];
            }
            else
            {
                currentStateSimsFlanagan = stateAtHalfTimeOfFlightBackwardPropagation;
            }

            // Check consistency between Sims-Flanagan trajectory and pseudo-impulsive shots trajectory.
            for ( int i = 0 ; i < 3 ; i++ )
            {
                BOOST_CHECK_SMALL( std::fabs( currentStateSimsFlanagan[ i ] - currentStateImpulsiveShots[ i ] )
                                   / currentStateImpulsiveShots.segment( 0,3 ).norm( ), 1.0e-6 );
                BOOST_CHECK_SMALL( std::fabs( currentStateSimsFlanagan[ i + 3 ] - currentStateImpulsiveShots[ i + 3 ] )
                        / currentStateImpulsiveShots.segment( 3,3 ).norm( ), 1.0e-6 );
            }
        }
    }

}

//! Test full propagation for Sims Flanagan (assuming constant thrust over each segment).
BOOST_AUTO_TEST_CASE( test_Sims_Flanagan_full_propagation )
{


    spice_interface::loadStandardSpiceKernels( );

    double maximumThrust = 0.8;
    double specificImpulse = 3000.0;
    double vehicleInitialMass = 1800.0;
    int numberSegments = 500;

    // Define (constant) specific impulse function.
    std::function< double( const double ) > specificImpulseFunction = [ = ]( const double currentTime )
    {
        return specificImpulse;
    };

    double julianDate = 1000.0 * physical_constants::JULIAN_DAY;
    double timeOfFlight = 700.0 * physical_constants::JULIAN_DAY;


    // Define state at departure and arrival.
    ephemerides::EphemerisPointer pointerToDepartureBodyEphemeris = std::make_shared< ephemerides::ApproximatePlanetPositions>(
                ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::earthMoonBarycenter );
    ephemerides::EphemerisPointer pointerToArrivalBodyEphemeris = std::make_shared< ephemerides::ApproximatePlanetPositions >(
                ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::mars );
    Eigen::Vector6d stateAtDeparture = pointerToDepartureBodyEphemeris->getCartesianState( julianDate );
    Eigen::Vector6d stateAtArrival = pointerToArrivalBodyEphemeris->getCartesianState( julianDate + timeOfFlight );


    // Define optimisation algorithm.
    algorithm optimisationAlgorithm{ pagmo::de1220() };
    std::shared_ptr< OptimisationSettings > optimisationSettings =
            std::make_shared< OptimisationSettings >( optimisationAlgorithm, 1, 10 );

    // Create Sims-Flanagan object
    double centralBodyGravitationalParameter = spice_interface::getBodyGravitationalParameter( "Sun" );
    SimsFlanagan simsFlanagan = SimsFlanagan(
                stateAtDeparture, stateAtArrival, centralBodyGravitationalParameter, vehicleInitialMass, maximumThrust,
                specificImpulseFunction, numberSegments, timeOfFlight, optimisationSettings );

    std::vector< double > bestIndividual = simsFlanagan.getBestIndividual( );

    std::vector< Eigen::Vector3d > bestThrottles;
    for ( int i = 0 ; i < numberSegments ; i++ )
    {
        bestThrottles.push_back( ( Eigen::Vector3d( ) << bestIndividual[ i * 3 ], bestIndividual[ i * 3 + 1 ],
                bestIndividual[ i * 3 + 2 ] ).finished( ) );
    }


    // Define integrator settings.
    double stepSize = ( timeOfFlight ) / 30000.0;
    std::shared_ptr< IntegratorSettings< double > > integratorSettings =
            std::make_shared< IntegratorSettings< double > > ( rungeKutta4, 0.0, stepSize );

    std::string bodyToPropagate = "Vehicle";
    std::string centralBody = "Sun";

    // Create system of bodies.
    SystemOfBodies bodies = getTestBodyMap( );
    bodies.at( bodyToPropagate )->setSuppressDependentOrientationCalculatorWarning( true );
    bodies.at( bodyToPropagate )->setConstantBodyMass( vehicleInitialMass );

    ///! Set-up Sims-Flanagan full propagation.
    // Create termination conditions settings.
    std::pair< std::shared_ptr< PropagationTerminationSettings >,
            std::shared_ptr< PropagationTerminationSettings > > terminationConditions;
    terminationConditions.first = std::make_shared< PropagationTimeTerminationSettings >( 0.0, true );
    terminationConditions.second = std::make_shared< PropagationTimeTerminationSettings >( timeOfFlight, true );


    // Create object with list of dependent variables
    std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariablesList;
    dependentVariablesList.push_back( std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                                          basic_astrodynamics::thrust_acceleration, "Vehicle", "Vehicle", 0 ) );
    std::shared_ptr< DependentVariableSaveSettings > dependentVariablesToSave =
            std::make_shared< DependentVariableSaveSettings >( dependentVariablesList, false );

    // Create complete propagation settings (backward and forward propagations).
    basic_astrodynamics::AccelerationMap perturbingAccelerationsMap;
    std::pair< std::shared_ptr< PropagatorSettings< double > >,
            std::shared_ptr< PropagatorSettings< double > > > propagatorSettings =
            simsFlanagan.createLowThrustPropagatorSettings(
                bodies, bodyToPropagate, centralBody, specificImpulseFunction, perturbingAccelerationsMap,
                integratorSettings, dependentVariablesToSave );

    // Compute full propagation.
    std::map< double, Eigen::VectorXd > fullPropagationResults;
    std::map< double, Eigen::Vector6d > simsFlanaganResults;
    std::map< double, Eigen::VectorXd > dependentVariablesHistory;
    simsFlanagan.computeSemiAnalyticalAndFullPropagation(
                bodies, integratorSettings, propagatorSettings, fullPropagationResults,
                simsFlanaganResults, dependentVariablesHistory );

    // Check consistency between Sims-Flanagan semi-analytical results and full propagation, at departure and arrival.
    Eigen::Vector6d stateAtDepartureSimsFlanagan = simsFlanaganResults.begin( )->second;
    Eigen::Vector6d stateAtDepartureFullPropagation = fullPropagationResults.begin( )->second;
    Eigen::Vector6d stateAtArrivalSimsFlanagan = simsFlanaganResults.rbegin( )->second;
    Eigen::Vector6d stateAtArrivalFullPropagation = fullPropagationResults.rbegin( )->second;

    for ( int i = 0 ; i < 3 ; i++ )
    {
        BOOST_CHECK_SMALL( ( std::fabs( stateAtDepartureSimsFlanagan[ i ] - stateAtDepartureFullPropagation[ i ] ) / stateAtDeparture.segment( 0, 3 ).norm( ) ), 1.0e-3 );
        BOOST_CHECK_SMALL( ( std::fabs( stateAtDepartureSimsFlanagan[ i + 3 ] - stateAtDepartureFullPropagation[ i + 3 ] ) / stateAtDeparture.segment( 3, 3 ).norm( ) ), 1.0e-2 );
        BOOST_CHECK_SMALL( ( std::fabs( stateAtArrivalSimsFlanagan[ i ] - stateAtArrivalFullPropagation[ i ] ) / stateAtArrival.segment( 0, 3 ).norm( ) ), 1.0e-3 );
        BOOST_CHECK_SMALL( ( std::fabs( stateAtArrivalSimsFlanagan[ i + 3 ] - stateAtArrivalFullPropagation[ i + 3 ] ) / stateAtArrival.segment( 3, 3 ).norm( ) ), 1.0e-2 );
    }

    /// Test trajectory, mass, thrust, and thrust acceleration profiles for Sims-Flanagan.

    std::vector< double > epochsVector;
    epochsVector.push_back( 0.0 );
    epochsVector.push_back( timeOfFlight / 4.0 );
    epochsVector.push_back( timeOfFlight / 2.0 - 10.0 );
    epochsVector.push_back( 3.0 * timeOfFlight / 4.0 );
    epochsVector.push_back( timeOfFlight );

    std::map< double, Eigen::Vector6d > trajectory;
    std::map< double, Eigen::VectorXd > massProfile;
    std::map< double, Eigen::VectorXd > thrustProfile;
    std::map< double, Eigen::VectorXd > thrustAccelerationProfile;

    simsFlanagan.getTrajectory( epochsVector, trajectory );
    simsFlanagan.getMassProfile( epochsVector, massProfile, specificImpulseFunction, integratorSettings );
    simsFlanagan.getThrustForceProfile( epochsVector, thrustProfile, specificImpulseFunction, integratorSettings );
    simsFlanagan.getThrustAccelerationProfile( epochsVector, thrustAccelerationProfile, specificImpulseFunction, integratorSettings );


    // Retrieve low-thrust trajectory nominal accelerations (thrust + central gravity accelerations).
    basic_astrodynamics::AccelerationMap lowThrustTrajectoryAccelerations =
            simsFlanagan.retrieveLowThrustAccelerationMap(
                bodies, bodyToPropagate, centralBody, specificImpulseFunction, integratorSettings );

    // Create mass rate models
    std::map< std::string, std::shared_ptr< basic_astrodynamics::MassRateModel > > massRateModels;
    massRateModels[ bodyToPropagate ] = createMassRateModel( bodyToPropagate, std::make_shared< FromThrustMassModelSettings >( 1 ),
                                                             bodies, lowThrustTrajectoryAccelerations );

    // Define list of dependent variables to save.
    dependentVariablesList.push_back( std::make_shared< SingleDependentVariableSaveSettings >(
                                          total_mass_rate_dependent_variables, bodyToPropagate ) );

    // Create object with list of dependent variables
    dependentVariablesToSave = std::make_shared< DependentVariableSaveSettings >( dependentVariablesList, false );


    for ( std::map< double, Eigen::Vector6d >::iterator itr = trajectory.begin( ) ; itr != trajectory.end( ) ; itr++ )
    {
        double currentEpoch = itr->first;

        // Create termination conditions settings.
        std::shared_ptr< PropagationTerminationSettings > terminationSettings =
                std::make_shared< PropagationTimeTerminationSettings >( currentEpoch, true );

        // Define translational state propagation settings
        std::shared_ptr< TranslationalStatePropagatorSettings< double > > translationalStatePropagatorSettings =
                std::make_shared< TranslationalStatePropagatorSettings< double > >
                ( std::vector< std::string >{ centralBody }, lowThrustTrajectoryAccelerations,
                  std::vector< std::string >{ bodyToPropagate }, stateAtDeparture,
                  terminationSettings, gauss_modified_equinoctial, dependentVariablesToSave );

        // Create settings for propagating the mass of the vehicle.
        std::shared_ptr< MassPropagatorSettings< double > > massPropagatorSettings =
                std::make_shared< MassPropagatorSettings< double > >(
                    std::vector< std::string >{ bodyToPropagate }, massRateModels, ( Eigen::Matrix< double, 1, 1 >( ) << vehicleInitialMass ).finished( ),
                    terminationSettings );

        integratorSettings->initialTimeStep_ = std::fabs( integratorSettings->initialTimeStep_ );
        integratorSettings->initialTime_ = 0.0;

        // Create list of propagation settings.
        std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > > propagatorSettingsVector;

        // Backward propagator settings vector.
        propagatorSettingsVector.push_back( translationalStatePropagatorSettings );
        propagatorSettingsVector.push_back( massPropagatorSettings );

        // Define propagator settings.
        std::shared_ptr< PropagatorSettings< double > > propagatorSettings = std::make_shared< MultiTypePropagatorSettings< double > >(
                    propagatorSettingsVector, terminationSettings, dependentVariablesToSave );

        bodies.at( bodyToPropagate )->setConstantBodyMass( vehicleInitialMass );

        // Perform the backward propagation.
        SingleArcDynamicsSimulator< > dynamicsSimulator( bodies, integratorSettings, propagatorSettings );
        std::map< double, Eigen::VectorXd > stateHistory = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
        std::map< double, Eigen::VectorXd > dependentVariableHistory = dynamicsSimulator.getDependentVariableHistory( );

        double currentExpectedMass = stateHistory.rbegin( )->second[ 6 ];
        Eigen::Vector3d currentExpectedThrustAcceleration =  dependentVariableHistory.rbegin( )->second.segment( 0, 3 );
        Eigen::Vector3d currentExpectedThrust = currentExpectedThrustAcceleration * currentExpectedMass;

        Eigen::Vector6d calculatedState = simsFlanagan.computeCurrentStateVector(
                    currentEpoch );
        double calculatedMass = simsFlanagan.computeCurrentMass( currentEpoch, specificImpulseFunction, integratorSettings );
        Eigen::Vector3d calculatedThrust = simsFlanagan.computeCurrentThrustForce(
                    currentEpoch, specificImpulseFunction, integratorSettings ).transpose( );
        Eigen::Vector3d calculatedThrustAcceleration = simsFlanagan.computeCurrentThrustAcceleration(
                    currentEpoch, specificImpulseFunction, integratorSettings );

        for ( int i = 0 ; i < 3 ; i++ )
        {
            // Check consistency between output of getTrajectory function and of computeCurrentStateVector.
            BOOST_CHECK_SMALL( ( std::fabs( trajectory[ itr->first ][ i ] - calculatedState[ i ] ) /
                               calculatedState.segment( 0, 3 ).norm( ) ), 1.0e-6 );
            BOOST_CHECK_SMALL( ( std::fabs( trajectory[ itr->first ][ i + 3 ] - calculatedState[ i + 3 ] ) /
                    calculatedState.segment( 3, 3 ).norm( ) ), 1.0e-6 );

            // Check consistency between Sims-Flanagan thrust profile and full propagation current thrust value.
            BOOST_CHECK_SMALL( ( std::fabs( thrustProfile[ itr->first ][ i ] - currentExpectedThrust[ i ] ) ), 1.0e-8 );

            // Idem with thrust acceleration.
            BOOST_CHECK_SMALL( ( std::fabs( thrustAccelerationProfile[ itr->first ][ i ] - currentExpectedThrustAcceleration[ i ] ) ), 1.0e-7 );

            // Check consistency between Sims-Flanagan computeThrustFunction and full propagation current thrust value.
            BOOST_CHECK_SMALL( ( std::fabs( calculatedThrust[ i ] - currentExpectedThrust[ i ] ) ), 1.0e-8 );

            // Idem with thrust acceleration.
            BOOST_CHECK_SMALL( ( std::fabs( calculatedThrustAcceleration[ i ] - currentExpectedThrustAcceleration[ i ] ) ), 1.0e-10 );
        }

        // Check consistency between Sims-Flanagan mass profile and the propagated mass.
        BOOST_CHECK_SMALL( std::fabs( massProfile[ itr->first ][ 0 ] - currentExpectedMass ) / currentExpectedMass, 1.0e-3 );

        // Check consistency between computeCurrentMass function and the propagated mass.
        BOOST_CHECK_SMALL( std::fabs( calculatedMass - currentExpectedMass ), 1.0e-15 );

    }


}


//! Test Sims-Flanagan implementation by comparing it with trajectory obtained with successive impulsive shots applied at times corresponding to
//! half of each of the Sims-Flanagan segments.
BOOST_AUTO_TEST_CASE( test_Sims_Flanagan_Shape_Based )
{
    spice_interface::loadStandardSpiceKernels( );

    double julianDate = 9264.5 * physical_constants::JULIAN_DAY;
    double timeOfFlight = 1000.0 * physical_constants::JULIAN_DAY;
    int numberOfRevolutions = 2;
    double maximumThrust = 5.0;
    double specificImpulse = 3000.0;
    double vehicleInitialMass = 2800.0;
    int numberSegments = 200;

    // Define specific impulse function.
    std::function< double( const double ) > specificImpulseFunction = [ = ]( const double ){ return specificImpulse; };

    // Retrieve cartesian state at departure and arrival.
    ephemerides::EphemerisPointer pointerToDepartureBodyEphemeris = std::make_shared< ephemerides::ApproximatePlanetPositions>(
                ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::earthMoonBarycenter );
    ephemerides::EphemerisPointer pointerToArrivalBodyEphemeris = std::make_shared< ephemerides::ApproximatePlanetPositions >(
                ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::mars );
    Eigen::Vector6d cartesianStateDepartureBody = pointerToDepartureBodyEphemeris->getCartesianState( julianDate );
    Eigen::Vector6d cartesianStateArrivalBody = pointerToArrivalBodyEphemeris->getCartesianState( julianDate + timeOfFlight );



    /// Shape-based trajectory Earth-Mars transfer.

    // Initialize free coefficients vectors
    Eigen::VectorXd freeCoefficientsRadialVelocityFunction = Eigen::VectorXd::Zero( 0 );
    Eigen::VectorXd freeCoefficientsNormalVelocityFunction = Eigen::VectorXd::Zero( 0 );
    Eigen::VectorXd freeCoefficientsAxialVelocityFunction = Eigen::VectorXd::Zero( 0 );


    double frequency = 2.0 * mathematical_constants::PI / timeOfFlight;
    double scaleFactor = 1.0 / timeOfFlight;

    // Create base function settings for the components of the radial velocity composite function.
    std::shared_ptr< BaseFunctionHodographicShapingSettings > firstRadialVelocityBaseFunctionSettings =
            std::make_shared< BaseFunctionHodographicShapingSettings >( );
    std::shared_ptr< BaseFunctionHodographicShapingSettings > secondRadialVelocityBaseFunctionSettings =
            std::make_shared< PowerFunctionHodographicShapingSettings >( 1.0, scaleFactor );
    std::shared_ptr< BaseFunctionHodographicShapingSettings > thirdRadialVelocityBaseFunctionSettings =
            std::make_shared< TrigonometricFunctionHodographicShapingSettings >( frequency );

    // Create components of the radial velocity composite function.
    std::vector< std::shared_ptr< BaseFunctionHodographicShaping > > radialVelocityFunctionComponents;
    radialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( constant, firstRadialVelocityBaseFunctionSettings ) );
    radialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPower, secondRadialVelocityBaseFunctionSettings ) );
    radialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( cosine, thirdRadialVelocityBaseFunctionSettings ) );

    // Create base function settings for the components of the normal velocity composite function.
    std::shared_ptr< BaseFunctionHodographicShapingSettings > firstNormalVelocityBaseFunctionSettings =
            std::make_shared< BaseFunctionHodographicShapingSettings >( );
    std::shared_ptr< BaseFunctionHodographicShapingSettings > secondNormalVelocityBaseFunctionSettings =
            std::make_shared< PowerFunctionHodographicShapingSettings >( 1.0, scaleFactor );
    std::shared_ptr< BaseFunctionHodographicShapingSettings > thirdNormalVelocityBaseFunctionSettings =
            std::make_shared< TrigonometricFunctionHodographicShapingSettings >( frequency );

    // Create components of the normal velocity composite function.
    std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > > normalVelocityFunctionComponents;
    normalVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( constant, firstNormalVelocityBaseFunctionSettings ) );
    normalVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPower, secondNormalVelocityBaseFunctionSettings ) );
    normalVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( cosine, thirdNormalVelocityBaseFunctionSettings ) );

    // Create base function settings for the components of the axial velocity composite function.
    std::shared_ptr< BaseFunctionHodographicShapingSettings > firstAxialVelocityBaseFunctionSettings =
            std::make_shared< TrigonometricFunctionHodographicShapingSettings >( ( numberOfRevolutions + 0.5 ) * frequency );
    std::shared_ptr< BaseFunctionHodographicShapingSettings > secondAxialVelocityBaseFunctionSettings =
            std::make_shared< PowerTimesTrigonometricFunctionHodographicShapingSettings >
            ( 3.0, ( numberOfRevolutions + 0.5 ) * frequency, scaleFactor );
    std::shared_ptr< BaseFunctionHodographicShapingSettings > thirdAxialVelocityBaseFunctionSettings =
            std::make_shared< PowerTimesTrigonometricFunctionHodographicShapingSettings >(
                3.0, ( numberOfRevolutions + 0.5 ) * frequency, scaleFactor );

    // Set components for the axial velocity function.
    std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > > axialVelocityFunctionComponents;
    axialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( cosine, firstAxialVelocityBaseFunctionSettings ) );
    axialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPowerCosine, secondAxialVelocityBaseFunctionSettings ) );
    axialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPowerSine, thirdAxialVelocityBaseFunctionSettings ) );


    // Create hodographic-shaping object with defined velocity functions and boundary conditions.
    double centralBodyGravitationalParameter = spice_interface::getBodyGravitationalParameter( "Sun" );
    std::shared_ptr< HodographicShaping >  hodographicShaping = std::make_shared< HodographicShaping >(
                cartesianStateDepartureBody, cartesianStateArrivalBody, timeOfFlight, centralBodyGravitationalParameter, numberOfRevolutions,
                radialVelocityFunctionComponents, normalVelocityFunctionComponents, axialVelocityFunctionComponents,
                freeCoefficientsRadialVelocityFunction, freeCoefficientsNormalVelocityFunction, freeCoefficientsAxialVelocityFunction,
                vehicleInitialMass );

    int numberOfSteps = 10000;
    double stepSize = timeOfFlight / static_cast< double >( numberOfSteps );

    // Define integrator settings.
    std::shared_ptr< IntegratorSettings< double > > integratorSettings =
            std::make_shared< IntegratorSettings< double > > ( rungeKutta4, 0.0, stepSize );


    // Calculate number of segments for both the forward propagation (from departure to match point)
    // and the backward propagation (from arrival to match point).
    int numberSegmentsForwardPropagation = ( numberSegments + 1 ) / 2;
    int numberSegmentsBackwardPropagation = numberSegments / 2;
    int segmentDurationForwardPropagation = timeOfFlight / ( 2.0 * numberSegmentsForwardPropagation );
    int segmentDurationBackwardPropagation = timeOfFlight / ( 2.0 * numberSegmentsBackwardPropagation );



    /// Set-up Sims-Flanagan trajectory using the thrust profile from hodographic shaping as initial guess.
    std::function< Eigen::Vector3d( const double ) > initialGuessThrustFromShaping =
            getInitialGuessFunctionFromShaping( hodographicShaping, numberSegments, timeOfFlight, specificImpulseFunction, integratorSettings );

    std::vector< double > initialGuessVector = convertToSimsFlanaganThrustModel(
                initialGuessThrustFromShaping, maximumThrust, timeOfFlight, numberSegmentsForwardPropagation, numberSegmentsBackwardPropagation );


    // Define optimisation algorithm (here we ensure the Sims-Flanagan solution stays close to the hodographic shaping initial guess).
    algorithm optimisationAlgorithm{ pagmo::de1220() };

    std::shared_ptr< OptimisationSettings > optimisationSettings =
            std::make_shared< OptimisationSettings >(
                optimisationAlgorithm, 1, 10, 1.0e-6, std::make_pair( initialGuessVector, 0.002 ) );

    // Create Sims-Flanagan trajectory.
    SimsFlanagan simsFlanagan = SimsFlanagan(
                cartesianStateDepartureBody, cartesianStateArrivalBody, centralBodyGravitationalParameter,
                vehicleInitialMass, maximumThrust, specificImpulseFunction, numberSegments, timeOfFlight, optimisationSettings );

    // Compare Sims-Flanagan and hodographic shaping for a given set of epochs.
    std::vector< double > epochsVector;
    int numberSteps = 10;
    for ( int i = 1 ; i <= numberSteps ; i++ )
    {
        epochsVector.push_back( timeOfFlight / numberSteps * i );
    }

    std::map< double, Eigen::Vector6d > simsFlanaganTrajectory;
    simsFlanagan.getTrajectory( epochsVector, simsFlanaganTrajectory );

    std::map< double, Eigen::Vector6d > hodographicShapingTrajectory;
    hodographicShaping->getTrajectory( epochsVector, hodographicShapingTrajectory );

    int counter = 0;
    // Check consistency between shaping method and Sims-Flanagan.
    for ( std::map< double, Eigen::Vector6d >::iterator itr = simsFlanaganTrajectory.begin( ) ; itr != simsFlanaganTrajectory.end( ) ; itr++ )
    {
        counter++;
        for ( int i = 0 ; i < 3 ; i++ )
        {
            BOOST_CHECK_SMALL( ( std::fabs( simsFlanaganTrajectory[ itr->first ][ i ] - hodographicShapingTrajectory[ itr->first ][ i ] )
                    / hodographicShapingTrajectory[ itr->first ].segment( 0, 3 ).norm( ) ), 5.0e-2 );
            BOOST_CHECK_SMALL( ( std::fabs( simsFlanaganTrajectory[ itr->first ][ i + 3 ] - hodographicShapingTrajectory[ itr->first ][ i + 3 ] )
                    / hodographicShapingTrajectory[ itr->first ].segment( 3, 3 ).norm( ) ), 1.0e-1 );
        }
    }

    std::vector< double > fitnessVector = simsFlanagan.getBestIndividualFitness( );
    std::vector< double > bestIndividual = simsFlanagan.getBestIndividual( );

    std::vector< Eigen::Vector3d > bestThrottles;
    for ( int i = 0 ; i < numberSegments ; i++ )
    {
        bestThrottles.push_back( ( Eigen::Vector3d( ) << bestIndividual[ i * 3 ], bestIndividual[ i * 3 + 1 ],
                bestIndividual[ i * 3 + 2 ] ).finished( ) );
    }

    // Check that the Sims-Flanagan throttles are within 0.5% of the initial guess.
    for ( int i = 0 ; i < numberSegmentsForwardPropagation ; i++ )
    {
        Eigen::Vector3d simsFlanaganThrust = bestThrottles[ i ] * maximumThrust;
        Eigen::Vector3d initialGuessThrust = initialGuessThrustFromShaping( ( i + 0.5 ) * segmentDurationForwardPropagation );
        for ( int j = 0 ; j < 3 ; j++ )
        {
            BOOST_CHECK_SMALL( std::fabs( ( simsFlanaganThrust[ j ] - initialGuessThrust[ j ] ) / initialGuessThrust[ j ] ), 5.0e-3 );
        }
    }
    for ( int i = 0 ; i < numberSegmentsBackwardPropagation ; i++ )
    {
        Eigen::Vector3d simsFlanaganThrust = bestThrottles[ i + numberSegmentsForwardPropagation ] * maximumThrust;
        Eigen::Vector3d initialGuessThrust = initialGuessThrustFromShaping( timeOfFlight / 2.0 + ( i + 0.5 ) * segmentDurationBackwardPropagation );
        for ( int j = 0 ; j < 3 ; j++ )
        {
            BOOST_CHECK_SMALL( std::fabs( ( simsFlanaganThrust[ j ] - initialGuessThrust[ j ] ) / initialGuessThrust[ j ] ), 5.0e-3 );
        }
    }
}


BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
