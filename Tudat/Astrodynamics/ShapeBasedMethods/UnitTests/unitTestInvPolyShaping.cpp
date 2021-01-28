/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      Wakker, K. F. (2007), Lecture Notes Astrodynamics II (Chapter 18), TU Delft course AE4-874,
 *          Delft University of technology, Delft, The Netherlands.
 *
 */

#define BOOST_TEST_MAIN

#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>
#include <Eigen/Dense>
#include <math.h>
#include <iostream>

#include "Tudat/Basics/testMacros.h"

#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"
#include "Tudat/Mathematics/BasicMathematics/basicMathematicsFunctions.h"
#include "Tudat/Mathematics/BasicMathematics/coordinateConversions.h"
#include "Tudat/Astrodynamics/Ephemerides/approximatePlanetPositions.h"
#include "Tudat/SimulationSetup/tudatSimulationHeader.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/celestialBodyConstants.h"
#include "Tudat/Astrodynamics/ShapeBasedMethods/baseFunctionsInvPolyShaping.h"
#include "Tudat/Astrodynamics/ShapeBasedMethods/compositeFunctionInvPolyShaping.h"
#include "Tudat/Astrodynamics/ShapeBasedMethods/invPolyShaping.h"

namespace tudat
{
namespace unit_tests
{

//! Test invpoly shaping implementation.
BOOST_AUTO_TEST_SUITE( test_invpoly_shaping )



//! Test.
BOOST_AUTO_TEST_CASE( test_invpoly_shaping_boundary_conditions)
{

    spice_interface::loadStandardSpiceKernels( );

    int numberOfRevolutions = 1;
    double julianDate = (2461771.5 - 2400000.5 - 51544.5) * physical_constants::JULIAN_DAY; //
    double  timeOfFlight = 878.0;

    double addJD = 124.0 * physical_constants::JULIAN_DAY;
    julianDate = julianDate + addJD;


    // Create central, departure and arrival bodies.
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Sun" );

    std::map< std::string, std::shared_ptr< simulation_setup::BodySettings > > bodySettings =
            simulation_setup::getDefaultBodySettings( bodiesToCreate );

    std::string frameOrigin = "SSB";
    std::string frameOrientation = "ECLIPJ2000";

    // Define central body ephemeris settings.
    bodySettings[ "Sun" ]->ephemerisSettings = std::make_shared< simulation_setup::ConstantEphemerisSettings >(
                ( Eigen::Vector6d( ) << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ).finished( ), frameOrigin, frameOrientation );

    bodySettings[ "Sun" ]->ephemerisSettings->resetFrameOrientation( frameOrientation );
    bodySettings[ "Sun" ]->rotationModelSettings->resetOriginalFrame( frameOrientation );

    // Create body map.
    simulation_setup::NamedBodyMap bodyMap = createBodies( bodySettings );

    bodyMap[ "Vehicle" ] = std::make_shared< simulation_setup::Body >( );
    bodyMap.at( "Vehicle" )->setEphemeris( std::make_shared< ephemerides::TabulatedCartesianEphemeris< > >(
                                                         std::shared_ptr< interpolators::OneDimensionalInterpolator
                                                         < double, Eigen::Vector6d > >( ), frameOrigin, frameOrientation ) );

    setGlobalFrameBodyEphemerides( bodyMap, frameOrigin, frameOrientation );

    // Set vehicle mass.
    bodyMap[ "Vehicle" ]->setConstantBodyMass( 400.0 );

    // Ephemeris departure body.
    ephemerides::EphemerisPointer pointerToDepartureBodyEphemeris = std::make_shared< ephemerides::ApproximatePlanetPositions>(
                ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::earthMoonBarycenter );

    // Ephemeris arrival body.
    ephemerides::EphemerisPointer pointerToArrivalBodyEphemeris = std::make_shared< ephemerides::ApproximatePlanetPositions >(
                ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::mars );


    // Define root finder settings (used to update the updated value of the free coefficient, so that it matches the required time of flight).
        std::shared_ptr< root_finders::RootFinderSettings > rootFinderSettings =
                std::make_shared< root_finders::RootFinderSettings >( root_finders::bisection_root_finder, 1.0e-6, 50 );


    // Initialise input arrays
    Eigen::VectorXd julianDateInput(9);
    Eigen::VectorXd timeOfFlightInput(9);
    Eigen::VectorXd nRevolutionsInput(9);

    // Fill input arrays
    julianDateInput   << 7355,  7805,   8255,   8055,   7605,   7305,   7805,   7305,   7905;
    timeOfFlightInput << 460,   1580,   580,    900,    3740,   540,    1180,   4260,   1460;
    nRevolutionsInput << 0,     0,      0,      1,      1,      1,      2,      2,      2;

    // Initialise computed variables
    Eigen::Vector6d computedInitialState;
    Eigen::Vector6d computedFinalState;
    double computedTimeOfFlight;
    double computedTravelledAngle;

    // Main loop - computing shape-base trajectories per input set and checn boundary conditions
    for (int idx = 0; idx < 9; ++idx)
    {
        julianDate          = julianDateInput(idx) * physical_constants::JULIAN_DAY;
        timeOfFlight        = timeOfFlightInput(idx);
        numberOfRevolutions = nRevolutionsInput(idx);


        // Define integrator settings.
        std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings =
                std::make_shared< numerical_integrators::IntegratorSettings< double > > (
                    numerical_integrators::rungeKutta4, 0.0, timeOfFlight * tudat::physical_constants::JULIAN_DAY / 500.0 );


        Eigen::Vector6d initialState = pointerToDepartureBodyEphemeris->getCartesianState( julianDate );
        Eigen::Vector6d finalState   = pointerToArrivalBodyEphemeris->getCartesianState( julianDate + timeOfFlight * tudat::physical_constants::JULIAN_DAY );


        // Compute shaped trajectory.
        shape_based_methods::InvPolyShaping invPolyShaping = shape_based_methods::InvPolyShaping(
                    initialState, finalState, timeOfFlight * tudat::physical_constants::JULIAN_DAY,
                    numberOfRevolutions, bodyMap, "Vehicle", "Sun", 0.0,
                    rootFinderSettings, -1.0e-3, 1.0, integratorSettings );

        // Compute boundary values
        computedTravelledAngle = invPolyShaping.getTravelledAzimuthAngle();
        computedInitialState   = invPolyShaping.computeCurrentStateVector(0.0);
        computedFinalState     = invPolyShaping.computeCurrentStateVector(computedTravelledAngle);
        computedTimeOfFlight   = invPolyShaping.computeTimeOfFlight()/tudat::physical_constants::JULIAN_DAY;

        // Check Difference
        // Position and velocity should be sub-meter(/per second) accurate
        // Time of flight should be sub-day accurate
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( computedInitialState, initialState, 0.1 ); // m,m/s
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( computedFinalState,   finalState,   0.1 ); // m,m/s

        BOOST_CHECK_SMALL( std::fabs(  computedTimeOfFlight - timeOfFlight ), 0.1 );  // day

        std::cout << "Difference Initial State [m,m/s]: " << std::endl << computedInitialState-initialState << std::endl;
        std::cout << "Difference Final State [m,m/s]:   " << std::endl << computedFinalState-finalState << std::endl;
        std::cout << "Difference Time Of Flight [days]: " << computedTimeOfFlight-timeOfFlight << std::endl;

    }


}


//! Test.
BOOST_AUTO_TEST_CASE( test_invpoly_shaping_earth_mars_transfer )
{

    spice_interface::loadStandardSpiceKernels( );

    int numberOfRevolutions = 1;
    double julianDate = (2461771.5 - 2400000.5 - 51544.5) * physical_constants::JULIAN_DAY; //
    double  timeOfFlight = 878.0;

    double addJD = 124.0 * physical_constants::JULIAN_DAY;
    julianDate = julianDate + addJD;


    // Create central, departure and arrival bodies.
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Sun" );

    std::map< std::string, std::shared_ptr< simulation_setup::BodySettings > > bodySettings =
            simulation_setup::getDefaultBodySettings( bodiesToCreate );

    std::string frameOrigin = "SSB";
    std::string frameOrientation = "ECLIPJ2000";

    // Define central body ephemeris settings.
    bodySettings[ "Sun" ]->ephemerisSettings = std::make_shared< simulation_setup::ConstantEphemerisSettings >(
                ( Eigen::Vector6d( ) << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ).finished( ), frameOrigin, frameOrientation );

    bodySettings[ "Sun" ]->ephemerisSettings->resetFrameOrientation( frameOrientation );
    bodySettings[ "Sun" ]->rotationModelSettings->resetOriginalFrame( frameOrientation );

    // Create body map.
    simulation_setup::NamedBodyMap bodyMap = createBodies( bodySettings );

    bodyMap[ "Vehicle" ] = std::make_shared< simulation_setup::Body >( );
    bodyMap.at( "Vehicle" )->setEphemeris( std::make_shared< ephemerides::TabulatedCartesianEphemeris< > >(
                                                         std::shared_ptr< interpolators::OneDimensionalInterpolator
                                                         < double, Eigen::Vector6d > >( ), frameOrigin, frameOrientation ) );

    setGlobalFrameBodyEphemerides( bodyMap, frameOrigin, frameOrientation );

    // Set vehicle mass.
    bodyMap[ "Vehicle" ]->setConstantBodyMass( 400.0 );

    // Ephemeris departure body.
    ephemerides::EphemerisPointer pointerToDepartureBodyEphemeris = std::make_shared< ephemerides::ApproximatePlanetPositions>(
                ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::earthMoonBarycenter );

    // Ephemeris arrival body.
    ephemerides::EphemerisPointer pointerToArrivalBodyEphemeris = std::make_shared< ephemerides::ApproximatePlanetPositions >(
                ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::mars );


    // Define root finder settings (used to update the updated value of the free coefficient, so that it matches the required time of flight).
        std::shared_ptr< root_finders::RootFinderSettings > rootFinderSettings =
                std::make_shared< root_finders::RootFinderSettings >( root_finders::bisection_root_finder, 1.0e-6, 50 );

    // Initialise reference variables
    double referenceDeltaV;
    double referenceTravelledAngle;
    double referenceThrustAcceleration;
    Eigen::Vector6d referenceState;
    Eigen::Vector3d referenceThrustAccelerationVector;

    // Initlialise computed variables
    double computedTravelledAngle;
    double computedDeltaV;
    double computedThrustAcceleration;
    Eigen::Vector6d computedState;
    Eigen::Vector3d computedThrustAccelerationVector;

    // Initialise input arrays
    Eigen::VectorXd julianDateInput(9);
    Eigen::VectorXd timeOfFlightInput(9);
    Eigen::VectorXd nRevolutionsInput(9);

    // Fill input arrays
    julianDateInput   << 7355,  7805,   8255,   8055,   7605,   7305,   7805,   7305,   7905;
    timeOfFlightInput << 460,   1580,   580,    900,    3740,   540,    1180,   4260,   1460;
    nRevolutionsInput << 0,     0,      0,      1,      1,      1,      2,      2,      2;

    // Initialise reference and computed value arrays
    Eigen::VectorXd deltaVOutput(9);
    Eigen::VectorXd travelledAzimuthAngleOutput(9);
    Eigen::MatrixXd stateVectorOutput(9,6);
    Eigen::MatrixXd thrustAccVectorOutput(9,3);
    Eigen::VectorXd thrustAccelerationOutput(9);

    Eigen::MatrixXd computedStateMatrix(9,6);
    Eigen::MatrixXd computedThrustAccMatrix(9,3);

    // Fill reference arrays
    stateVectorOutput << -87300823973.59428,   87612058918.1752,   3981727047.13507, -8125.693581317898, -7518.741331042373,  42.82533629447528,
                         -85232192822.12103, -80484352414.53401,  414525754.2481433,  8771.014457906886, -7771.375694441038, -378.4177738508685,
                          51196907479.03497, -91777253831.15215, -3179710538.490308,  11036.58498591008,  6941.697321447871, -126.4972111895747,
                          101917580218.7518,  26225206844.57833, -1960261136.034638,  -2555.96163676985,  12767.86843684987,  329.9468198410945,
                         -123560888959.3799,   7312465806.80676,  3192879569.775781, -265.3933486549732, -11058.19027769334, -224.7688302498224,
                         -105622947240.2387,  65765955417.79714,  3975956914.937587, -5949.207807479625, -9250.316750253323, -46.97701048287757,
                          90650985892.17877,  57589533410.87752, -1026950952.557125, -6034.762691201258,  11260.34864189498,  384.0495674403836,
                          33332529432.83597, -100922615340.2281, -2930994692.668556,  11959.17429768256,  4840.596298621478, -192.9645743675748,
                         -99025253005.21458, -65906309422.70737,  1058922693.824453,   7164.59679739856, -9049.294527302021, -365.6066851536302;

    thrustAccVectorOutput <<  -0.0001048972330249935, -9.706188814088704e-005,  6.453285739569895e-005,
                              -0.0005347580629137853,   0.0004738113056908931,  9.939610294554176e-005,
                              -0.0001086647707631111, -6.834704295803814e-005, -7.127941861093298e-005,
                             -1.215998904499108e-005,  6.074314187132768e-005, -8.740595653104765e-005,
                               7.77420575716755e-006,   0.0003239291676162617,  7.459374261824496e-005,
                             -4.012060752405178e-005, -6.238281462339836e-005,  6.678042056921115e-005,
                             -3.509780305404573e-005,  6.548948470324217e-005, -5.661651527002955e-005,
                             -4.555077900374157e-005, -1.843713677520056e-005, -3.234730376373473e-005,
                             -3.235844403459029e-005,  4.087056100917625e-005,  3.301550681014269e-005;

    deltaVOutput << 6717.077201659434, 20940.47591060687, 6758.832524431853, 7966.577467472755, 21506.89776773543, 7068.485124002492, 9075.207303175774, 24791.67202660631, 10934.75469425099;
    travelledAzimuthAngleOutput << 6.001538705590269, 6.087142884791357, 5.995702329960109, 10.72199109153088, 8.767955222622113, 7.118162154902879, 15.32112885385222, 15.84804858021657, 16.81378506845654;
    thrustAccelerationOutput << 0.0001568085753695488, 0.0007213479912861569, 0.000146833600390442, 0.0001071326087383385, 0.0003324977749482183, 9.980432495556996e-005, 9.341390788976907e-005, 5.883153526042376e-005, 6.170490537611324e-005;

    // Main loop - compute shape-based trajectory per input set
    for (int idx = 0; idx < 9; ++idx)
    {
        // Extract input values
        julianDate          = julianDateInput(idx) * physical_constants::JULIAN_DAY;
        timeOfFlight        = timeOfFlightInput(idx);
        numberOfRevolutions = nRevolutionsInput(idx);

        // Extract reference values
        referenceState      = stateVectorOutput.row(idx);
        referenceDeltaV     = deltaVOutput(idx);
        referenceTravelledAngle           = travelledAzimuthAngleOutput(idx);
        referenceThrustAcceleration       = thrustAccelerationOutput(idx);
        referenceThrustAccelerationVector = thrustAccVectorOutput.row(idx);

        // Define integrator settings.
        std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings =
                std::make_shared< numerical_integrators::IntegratorSettings< double > > (
                    numerical_integrators::rungeKutta4, 0.0, timeOfFlight * tudat::physical_constants::JULIAN_DAY / 500.0 );


        Eigen::Vector6d initialState = pointerToDepartureBodyEphemeris->getCartesianState( julianDate );
        Eigen::Vector6d finalState   = pointerToArrivalBodyEphemeris->getCartesianState( julianDate + timeOfFlight * tudat::physical_constants::JULIAN_DAY );


        // Compute shaped trajectory.
        shape_based_methods::InvPolyShaping invPolyShaping = shape_based_methods::InvPolyShaping(
                    initialState, finalState, timeOfFlight * tudat::physical_constants::JULIAN_DAY,
                    numberOfRevolutions, bodyMap, "Vehicle", "Sun", 0.0,
                    rootFinderSettings, -1.0e-3, 1.0, integratorSettings );

        // Compute performance parameters
        computedTravelledAngle = invPolyShaping.getTravelledAzimuthAngle();
        computedState          = invPolyShaping.computeCurrentStateVector(computedTravelledAngle)/2.0;
        computedDeltaV         = invPolyShaping.computeDeltaV();

        computedThrustAcceleration       = invPolyShaping.computeCurrentThrustAccelerationMagnitude(computedTravelledAngle)/2.0;
        computedThrustAccelerationVector = invPolyShaping.computeCurrentThrustAccelerationVector(computedTravelledAngle)/2.0;

        if (invPolyShaping.getInfeasibleTOF()==false)
        {
            for (int idx2 = 0; idx2 < 6; ++idx2)
            {
                computedStateMatrix(idx,idx2) = computedState(idx2);
            }
            for (int idx2 = 0; idx2 < 3; ++idx2)
            {
                computedThrustAccMatrix(idx,idx2) = computedThrustAccelerationVector(idx2);
            }
            // Check difference between reference and computed values
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION( computedState, referenceState, 0.1 ); // m,m/s
            BOOST_CHECK_SMALL( std::fabs(  computedDeltaV - referenceDeltaV ), 0.1 );  // m/s
            BOOST_CHECK_SMALL( std::fabs(  computedTravelledAngle - referenceTravelledAngle ), 0.1 );  // rad
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION( computedThrustAccelerationVector, referenceThrustAccelerationVector, 0.1 ); //
            BOOST_CHECK_SMALL( std::fabs(  computedThrustAcceleration - referenceThrustAcceleration ), 0.1 );  //

//            std::cout << "Index: " << idx << std::endl;

//            std::cout << std::setprecision( 16 ) << "State Position comp [m]: " << std::endl << computedStateMatrix << std::endl;
//            std::cout << std::setprecision( 16 ) << std::endl << computedDeltaV << std::endl;
//            std::cout << std::setprecision( 16 ) << std::endl << computedTravelledAngle << std::endl;
//            std::cout << std::setprecision( 16 ) << std::endl << computedThrustAcceleration << std::endl;
//            std::cout << std::setprecision( 16 ) << std::endl << computedThrustAccelerationVector.norm() << std::endl;
//            std::cout << std::setprecision( 16 ) << std::endl << computedThrustAccMatrix << std::endl;


            std::cout << "Difference Position [m]: " << std::endl << computedState.segment(0, 3)-referenceState.segment(0, 3) << std::endl;
            std::cout << "Difference Velocity [m/s]:   " << std::endl << computedState.segment(3, 3)-referenceState.segment(3, 3) << std::endl;
            std::cout << "Difference Travelled Azimuth Angle [rad]: " << computedTravelledAngle-referenceTravelledAngle << std::endl;
            std::cout << "Difference Delta V [m/s]: " << computedDeltaV-referenceDeltaV << std::endl; //
            std::cout << "Difference Thrust Acc. Magn. [m/s^2]: " << computedThrustAcceleration-referenceThrustAcceleration << std::endl; //
            std::cout << "Difference Thrust Acc. Vector [m/s^2] " << std::endl << computedThrustAccelerationVector-referenceThrustAccelerationVector << std::endl; //

        }
    }
}



BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
