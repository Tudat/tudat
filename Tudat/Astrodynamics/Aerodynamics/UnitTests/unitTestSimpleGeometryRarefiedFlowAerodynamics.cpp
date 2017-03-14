/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    semiMajorAixscopy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_MAIN

#include <boost/test/unit_test.hpp>

#include <Tudat/SimulationSetup/tudatSimulationHeader.h>

#include "applicationOutput.h"
#include "Tudat/InputOutput/basicInputOutput.h"


namespace tudat
{

namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_simple_geometry_rarefied_flow_aerodynamic_coefficients )

BOOST_AUTO_TEST_CASE( testSimpleGeometryRarefiedFlowCoefficients )
{
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            USING STATEMENTS              //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    using namespace tudat;
    using namespace simulation_setup;
    using namespace propagators;
    using namespace numerical_integrators;
    using namespace orbital_element_conversions;
    using namespace basic_mathematics;
    using namespace basic_astrodynamics;
    using namespace gravitation;
    using namespace numerical_integrators;
    using namespace interpolators;
    using namespace input_output;
    using namespace unit_conversions;


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////     CREATE ENVIRONMENT AND VEHICLE       //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Load Spice kernels.
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "pck00009.tpc" );
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "de-403-masses.tpc" );
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "de421.bsp" );

    // Set simulation time settings.
    const double simulationStartEpoch = 0.0;
    const double simulationEndEpoch = 200.0;

    for( unsigned int useSphereShape = 0; useSphereShape < 2; useSphereShape++ )
    {
        // Define body settings for simulation.
        std::vector< std::string > bodiesToCreate;
        bodiesToCreate.push_back( "Earth" );

        // Create body objects.
        std::map< std::string, boost::shared_ptr< BodySettings > > bodySettings =
                getDefaultBodySettings( bodiesToCreate );
        bodySettings[ "Earth" ]->atmosphereSettings = boost::make_shared< AtmosphereSettings >( nrlmsise00 );
        NamedBodyMap bodyMap = createBodies( bodySettings );

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////             CREATE VEHICLE            /////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // Create spacecraft object.
        bodyMap[ "Vehicle" ] = boost::make_shared< simulation_setup::Body >( );
        bodyMap[ "Vehicle" ]->setConstantBodyMass( 400.0 );
        double referenceArea = 10.0;

        // Aerodynamics interface
        bodyMap[ "Vehicle" ]->setAerodynamicCoefficientInterface(
                    createAerodynamicCoefficientInterface(
                        boost::make_shared< RarefiedFlowSimpleGeometryAerodynamicCoefficientSettings >(
                            referenceArea, useSphereShape ), "Vehicle" ) );

        // Finalize body creation.
        setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );


        /////////////////////////////  Define input varibales for coefficients here;
        std::vector< double > coefficientIndependentVariables;// = ...
        bodyMap[ "Vehicle" ]->getAerodynamicCoefficientInterface( )->updateFullCurrentCoefficients(
                    coefficientIndependentVariables );
        Eigen::Vector6d currentAerodynamicCoefficients =
                bodyMap[ "Vehicle" ]->getAerodynamicCoefficientInterface( )->getCurrentAerodynamicCoefficients( );

        ///////////////// Check currentAerodynamicCoefficients: entry 0 is drag coefficient, rest should be zero.

        {
            SelectedAccelerationMap accelerationMap;
            std::vector< std::string > bodiesToPropagate;
            std::vector< std::string > centralBodies;


            // Define acceleration model settings.
            std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > accelerationsOfVehicle;
            accelerationsOfVehicle[ "Earth" ].push_back( boost::make_shared< AccelerationSettings >( central_gravity ) );
            accelerationsOfVehicle[ "Earth" ].push_back( boost::make_shared< AccelerationSettings >( aerodynamic ) );

            accelerationMap[ "Vehicle" ] = accelerationsOfVehicle;

            bodiesToPropagate.push_back( "Vehicle" );
            centralBodies.push_back( "Earth" );

            // Create acceleration models
            basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                        bodyMap, accelerationMap, bodiesToPropagate, centralBodies );

            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            ///////////////////////             CREATE PROPAGATION SETTINGS            ////////////////////////////////////////////
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            // Set spherical elements for Vehicle.
            Eigen::Vector6d VehicleInitialState;
            VehicleInitialState( SphericalOrbitalStateElementIndices::radiusIndex ) =
                    spice_interface::getAverageRadius( "Earth" ) + 600.0E3;
            VehicleInitialState( SphericalOrbitalStateElementIndices::latitudeIndex ) = convertDegreesToRadians( 10.0 );
            VehicleInitialState( SphericalOrbitalStateElementIndices::longitudeIndex ) = convertDegreesToRadians( 20.0 );
            VehicleInitialState( SphericalOrbitalStateElementIndices::speedIndex ) = 7500.0;
            VehicleInitialState( SphericalOrbitalStateElementIndices::flightPathIndex ) = convertDegreesToRadians( 20.0 );
            VehicleInitialState( SphericalOrbitalStateElementIndices::headingAngleIndex ) = convertDegreesToRadians( 90.0 );

            // Convert Vehicle state from spherical elements to Cartesian elements.
            Eigen::Vector6d systemInitialState = convertSphericalOrbitalToCartesianState(
                        VehicleInitialState );
            boost::shared_ptr< ephemerides::RotationalEphemeris > earthRotationalEphemeris =
                    bodyMap.at( "Earth" )->getRotationalEphemeris( );
            systemInitialState = transformStateToGlobalFrame( systemInitialState, simulationStartEpoch, earthRotationalEphemeris );

            // Define list of dependent variables to save.
            std::vector< boost::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariablesToSave;

            /////////////////////////// Extend list of variables to save (latitude, longitude, .... needed for test.
            dependentVariablesToSave.push_back( boost::make_shared< SingleDependentVariableSaveSettings >(
                                                    airspeed_dependent_variable, "Vehicle" ) );
            dependentVariablesToSave.push_back( boost::make_shared< SingleDependentVariableSaveSettings >(
                                                    altitude_dependent_variable, "Vehicle" ) );
            dependentVariablesToSave.push_back( boost::make_shared< SingleDependentVariableSaveSettings >(
                                                    aerodynamic_moment_coefficients_dependent_variable, "Vehicle" ) );

            boost::shared_ptr< DependentVariableSaveSettings > dependentVariableSaveSettings =
                    boost::make_shared< DependentVariableSaveSettings >( dependentVariablesToSave );



            // Create propagation settings.
            boost::shared_ptr< TranslationalStatePropagatorSettings < double > > propagatorSettings =
                    boost::make_shared< TranslationalStatePropagatorSettings< double > >
                    ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState, simulationEndEpoch,
                      cowell, dependentVariableSaveSettings );

            boost::shared_ptr< IntegratorSettings< > > integratorSettings =
                    boost::make_shared< IntegratorSettings< > >
                    ( rungeKutta4, simulationStartEpoch, 30.0 );


            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            ///////////////////////             PROPAGATE ORBIT            ////////////////////////////////////////////////////////
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


            // Create simulation object and propagate dynamics.
            SingleArcDynamicsSimulator< > dynamicsSimulator(
                        bodyMap, integratorSettings, propagatorSettings, true, false, false );
        }
    }

}

BOOST_AUTO_TEST_SUITE_END( )

}

}
