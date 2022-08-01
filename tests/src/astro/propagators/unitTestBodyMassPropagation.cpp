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

#include "tudat/astro/basic_astro/massRateModel.h"
#include "tudat/simulation/propagation_setup//propagationSettings.h"
#include "tudat/simulation/environment_setup/body.h"
#include "tudat/simulation/propagation_setup/dynamicsSimulator.h"


namespace tudat
{

namespace unit_tests
{

using namespace tudat::propagators;
using namespace tudat::simulation_setup;
using namespace tudat::numerical_integrators;

BOOST_AUTO_TEST_SUITE( test_body_mass_propagation )

double getDummyMassRate1(
        const SystemOfBodies bodies )
{
    return ( bodies.at( "Vehicle1" )->getBodyMass( ) + 2.0 * bodies.at( "Vehicle2" )->getBodyMass( ) ) / 1.0E4;
}

double getDummyMassRate2(
        const SystemOfBodies bodies )
{
    return ( 3.0 * bodies.at( "Vehicle1" )->getBodyMass( ) + 2.0 * bodies.at( "Vehicle2" )->getBodyMass( ) ) / 1.0E4;
}

// Test mass rate of single body, linearly decreasing with time
BOOST_AUTO_TEST_CASE( testBodyMassPropagation )
{
    // Crate bodies
    SystemOfBodies bodies;
    bodies.createEmptyBody( "Vehicle" );

    // Create mass rate model.
    std::map< std::string, std::shared_ptr< basic_astrodynamics::MassRateModel > > massRateModels;
    massRateModels[ "Vehicle" ] = std::make_shared< basic_astrodynamics::CustomMassRateModel >(
                [ ]( const double ){ return -0.01; } );

    // Create settings for propagation
    Eigen::VectorXd initialMass = Eigen::VectorXd( 1 );
    initialMass( 0 ) = 500.0;
    std::shared_ptr< PropagatorSettings< double > > propagatorSettings =
            std::make_shared< MassPropagatorSettings< double > >(
                std::vector< std::string >{ "Vehicle" }, massRateModels, initialMass,
                std::make_shared< PropagationTimeTerminationSettings >( 1000.0 ) );

    // Define numerical integrator settings.
    std::shared_ptr< IntegratorSettings< > > integratorSettings =
            std::make_shared< IntegratorSettings< > >( rungeKutta4, 0.0, 1.0 );

    // Create dynamics simulation object.
    SingleArcDynamicsSimulator< double, double > dynamicsSimulator(
                bodies, integratorSettings, propagatorSettings, true, false, false );

    // Test propagated solution.
    std::map< double, Eigen::VectorXd > integratedState = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
    for( std::map< double, Eigen::VectorXd >::const_iterator stateIterator = integratedState.begin( );
         stateIterator != integratedState.end( ); stateIterator++ )
    {
        BOOST_CHECK_CLOSE_FRACTION( stateIterator->second( 0 ), 500.0 - 0.01 * stateIterator->first, 1.0E-13 );
    }
}

// Test coupled mass rate of two bodies. Model ius unphysical, but has an analytical solution, and allows the internal
// workings of the mass propagation to be more rigorously tested.
BOOST_AUTO_TEST_CASE( testTwoBodyMassPropagation )
{
    // Crate bodies
    SystemOfBodies bodies;
    bodies.createEmptyBody( "Earth" );
    bodies.createEmptyBody( "Vehicle1" );
    bodies.createEmptyBody( "Vehicle2" );

    // Create mass rate models.
    std::map< std::string, std::shared_ptr< basic_astrodynamics::MassRateModel > > massRateModels;
    massRateModels[ "Vehicle1" ] = std::make_shared< basic_astrodynamics::CustomMassRateModel >(
                std::bind( &getDummyMassRate1, bodies ) );
    massRateModels[ "Vehicle2" ] = std::make_shared< basic_astrodynamics::CustomMassRateModel >(
                std::bind( &getDummyMassRate2, bodies ) );
    bodies.at( "Earth" )->setEphemeris( std::make_shared< ephemerides::ConstantEphemeris >(
                                          [ ]( ){ return Eigen::Vector6d::Zero( ); } ) );
    bodies.at( "Vehicle1" )->setEphemeris( std::make_shared< ephemerides::ConstantEphemeris >(
                                          [ ]( ){ return Eigen::Vector6d::Zero( ); }, "Earth" ) );
    bodies.at( "Vehicle2" )->setEphemeris( std::make_shared< ephemerides::ConstantEphemeris >(
                                         [ ]( ){ return Eigen::Vector6d::Zero( ); }, "Earth" ) );

    // Create settings for propagation
    Eigen::VectorXd initialMass = Eigen::VectorXd( 2 );
    initialMass( 0 ) = 500.0;
    initialMass( 1 ) = 1000.0;
    std::shared_ptr< PropagatorSettings< double > > propagatorSettings =
            std::make_shared< MassPropagatorSettings< double > >(
                std::vector< std::string >{ "Vehicle1", "Vehicle2" }, massRateModels, initialMass,
                std::make_shared< PropagationTimeTerminationSettings >( 1000.0 ) );

    // Define numerical integrator settings.
    std::shared_ptr< IntegratorSettings< > > integratorSettings =
            std::make_shared< IntegratorSettings< > >( rungeKutta4, 0.0, 1.0 );

    // Create dynamics simulation object.
    SingleArcDynamicsSimulator< double, double > dynamicsSimulator(
                bodies, integratorSettings, propagatorSettings, true, false, true );

    // Test propagated solution.
    std::map< double, Eigen::VectorXd > integratedState = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
    for( std::map< double, Eigen::VectorXd >::const_iterator stateIterator = integratedState.begin( );
         stateIterator != integratedState.end( ); stateIterator++ )
    {
        // Test directly
        BOOST_CHECK_CLOSE_FRACTION(
                    stateIterator->second( 0 ),
                    100.0 * ( -std::exp( -stateIterator->first / 1E4 ) +
                              6.0 * std::exp( 4.0 * stateIterator->first / 1E4 ) ), 1.0E-13 );
        BOOST_CHECK_CLOSE_FRACTION(
                    stateIterator->second( 1 ),
                    100.0 * ( std::exp( -stateIterator->first / 1E4 ) +
                              9.0 * std::exp( 4.0 * stateIterator->first / 1E4 ) ), 1.0E-13 );

        // Test reset mass solution of vehicles.
        bodies.at( "Vehicle1" )->updateMass( stateIterator->first );
        bodies.at( "Vehicle2" )->updateMass( stateIterator->first );

        BOOST_CHECK_CLOSE_FRACTION( stateIterator->second( 0 ), bodies.at( "Vehicle1" )->getBodyMass( ),
                std::numeric_limits< double >::epsilon( ) );
        BOOST_CHECK_CLOSE_FRACTION( stateIterator->second( 1 ), bodies.at( "Vehicle2" )->getBodyMass( ),
                std::numeric_limits< double >::epsilon( ) );
    }
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat
