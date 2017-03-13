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

BOOST_AUTO_TEST_SUITE( test_aerodynamic_acceleration_force_moment_models )

BOOST_AUTO_TEST_CASE( testTabulatedDragCoefficient )
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


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////     CREATE ENVIRONMENT AND VEHICLE       //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Load Spice kernels.
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "pck00009.tpc" );
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "de-403-masses.tpc" );
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "de421.bsp" );

    // Set simulation time settings.
    const double simulationStartEpoch = 0.0;
    const double simulationEndEpoch = tudat::physical_constants::JULIAN_DAY;

    for( unsigned int useSphereShape = 0; useSphereShape < 2; useSphereShape++ )
    {
        // Define body settings for simulation.
        std::vector< std::string > bodiesToCreate;
        bodiesToCreate.push_back( "Earth" );

        // Create body objects.
        std::map< std::string, boost::shared_ptr< BodySettings > > bodySettings =
                getDefaultBodySettings( bodiesToCreate );
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
                            referenceArea, true ), "Vehicle" ) );

        // Finalize body creation.
        setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );
    }

}

BOOST_AUTO_TEST_SUITE_END( )

}

}
