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

#include <iomanip>


#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"
#include "Tudat/Mathematics/BasicMathematics/basicMathematicsFunctions.h"
#include "Tudat/Mathematics/BasicMathematics/coordinateConversions.h"
#include "Tudat/Astrodynamics/Ephemerides/approximatePlanetPositions.h"
#include "Tudat/SimulationSetup/tudatSimulationHeader.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/celestialBodyConstants.h"
//#include "Tudat/Astrodynamics/ShapeBasedMethods/baseFunctionsExposinsShaping.h"
#include "Tudat/Astrodynamics/ShapeBasedMethods/compositeFunctionExposinsShaping.h"
#include "Tudat/Astrodynamics/ShapeBasedMethods/exposinsShaping.h"

namespace tudat
{
namespace unit_tests
{

//! Test exposins shaping implementation.
BOOST_AUTO_TEST_SUITE( test_exposins_shaping )

//! Test.
BOOST_AUTO_TEST_CASE( test_exposins_shaping_intermediate_results )
{
    std::cout << "Running test_exposins_shaping_intermediate_results ..." << std::endl;

    spice_interface::loadStandardSpiceKernels( );

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

    // Initialise boundary values
    double travelledAngle;
    double computedTimeOfFlight;
    Eigen::Vector6d computedInitialState;
    Eigen::Vector6d computedFinalState;

    // Shape-based Inputs
    int numberOfRevolutions;
    double julianDate;
    double  timeOfFlight;
    double windingParameter;

    // Shape-based inputs arrays
    Eigen::VectorXd julianDateInput(9);
    Eigen::VectorXd timeOfFlightInput(9);
    Eigen::VectorXd nRevolutionsInput(9);
    Eigen::VectorXd windingParameterInput(9);

    // Filling input arrays
    julianDateInput       << 7355,  7805,   8255,   8055,   7605,   7305,   7805,   7805,   7805;
    timeOfFlightInput     << 460,   1580,   580,    900,    3740,   540,    1180,   1180,   1180;
    nRevolutionsInput     << 0,     0,      0,      1,      1,      1,      2,      2,      2;
    windingParameterInput << 0.1,   0.2,    0.3,    0.4,    0.5,    0.6,    0.7,    0.8,    0.9;

    // Loop throuh 9 input sets
    for (int idx = 0; idx < 9; ++idx)
    {
        // Extract inputs from arrays
        julianDate          = julianDateInput(idx) * physical_constants::JULIAN_DAY;
        timeOfFlight        = timeOfFlightInput(idx);
        numberOfRevolutions = nRevolutionsInput(idx);
        windingParameter    = windingParameterInput(idx);

        // Define integrator settings.
        std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings =
                std::make_shared< numerical_integrators::IntegratorSettings< double > > (
                    numerical_integrators::rungeKutta4, 0.0, timeOfFlight * tudat::physical_constants::JULIAN_DAY / 500.0 );

        // Define initial and final state
        Eigen::Vector6d initialState = pointerToDepartureBodyEphemeris->getCartesianState( julianDate );
        Eigen::Vector6d finalState   = pointerToArrivalBodyEphemeris->getCartesianState( julianDate + timeOfFlight * tudat::physical_constants::JULIAN_DAY );

        // Compute shaped trajectory.
        shape_based_methods::ExposinsShaping exposinsShaping = shape_based_methods::ExposinsShaping(
                    initialState, finalState, timeOfFlight * tudat::physical_constants::JULIAN_DAY,
                    numberOfRevolutions, bodyMap, "Vehicle", "Sun", windingParameter,
                    rootFinderSettings, integratorSettings );


        // Compute boundary conditions
        travelledAngle       = exposinsShaping.getTravelledAzimuthAngle();
        computedInitialState = exposinsShaping.computeCurrentStateVector(0.0);
        computedFinalState   = exposinsShaping.computeCurrentStateVector(travelledAngle);
        computedTimeOfFlight = exposinsShaping.computeTimeOfFlight()/tudat::physical_constants::JULIAN_DAY;


        if (exposinsShaping.getInfeasibleTOF()==false)
        {
            // Check boudnary conditions
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION( computedInitialState.segment(0, 3), initialState.segment(0, 3), 0.1 ); // m,m/s
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION( computedFinalState.segment(0, 3),   finalState.segment(0, 3),   0.1 ); // m,m/s

            BOOST_CHECK_SMALL( std::fabs(  computedTimeOfFlight - timeOfFlight ), 0.1 );  // day

            // Show boundary value differences
            std::cout << "Index: " << idx << std::endl;
            std::cout << "Difference Initial Position [m]: " << std::endl << computedInitialState.segment(0, 3)-initialState.segment(0, 3) << std::endl;
            std::cout << "Difference Final Position [m]:   " << std::endl << computedFinalState.segment(0, 3)-finalState.segment(0, 3) << std::endl;
            std::cout << "Difference Time Of Flight [days]: " << computedTimeOfFlight-timeOfFlight << std::endl;
        }
    }
}

//! Test.
BOOST_AUTO_TEST_CASE( test_exposins_shaping_boundary_conditions )
{
    std::cout << "Running test_exposins_shaping_boundary_conditions ..." << std::endl;

    spice_interface::loadStandardSpiceKernels( );

    double julianDate = (4115 - 51544.5 + 343.5 -57) * physical_constants::JULIAN_DAY;
    double  timeOfFlight = 110.448; // to get -80 degrees for requried gamma

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

    // Define integrator settings.
    std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings =
            std::make_shared< numerical_integrators::IntegratorSettings< double > > (
                numerical_integrators::rungeKutta4, 0.0, timeOfFlight * tudat::physical_constants::JULIAN_DAY / 500.0 );

    // Define root finder settings (used to update the updated value of the free coefficient, so that it matches the required time of flight).
        std::shared_ptr< root_finders::RootFinderSettings > rootFinderSettings =
                std::make_shared< root_finders::RootFinderSettings >( root_finders::bisection_root_finder, 1.0e-6, 50 );


    Eigen::Vector6d initialState;
    Eigen::Vector6d finalState;

    // r_1 = 1.0, Psi = 90 deg, gamma_1 = -80 deg (requiredGamma)
    initialState(0) = 1.0 * physical_constants::ASTRONOMICAL_UNIT;
    initialState(1) = 0.0 * physical_constants::ASTRONOMICAL_UNIT;
    initialState(2) = 0.0 * physical_constants::ASTRONOMICAL_UNIT;

    initialState(3) = 0.0;
    initialState(4) = 0.0;
    initialState(5) = 0.0;

    // r_2 = 1.5, Psi = 90 deg
    finalState(0) = 0.0 * physical_constants::ASTRONOMICAL_UNIT;
    finalState(1) = 1.5 * physical_constants::ASTRONOMICAL_UNIT;
    finalState(2) = 0.0 * physical_constants::ASTRONOMICAL_UNIT;

    finalState(3) = 0.0;
    finalState(4) = 0.0;
    finalState(5) = 0.0;

    // k_2 = 1/12, Nrev = 2
    int numberOfRevolutions = 2;
    double windingParameter = 1.0/12.0;

    // Compute shaped trajectory.
    shape_based_methods::ExposinsShaping exposinsShaping = shape_based_methods::ExposinsShaping(
                initialState, finalState, timeOfFlight * tudat::physical_constants::JULIAN_DAY,
                numberOfRevolutions, bodyMap, "Vehicle", "Sun", windingParameter,
                rootFinderSettings, integratorSettings );



    // Initialise reference variables
    double referenceAzimuthAngle;
    double referenceRadius;
    double referenceAcceleration;
    double referenceInvAngAcc;

    // Initialise reference value arrays
    Eigen::VectorXd referenceAzimuthAngleData(20);
    Eigen::VectorXd referenceRadiusData(20);
    Eigen::VectorXd referenceAccelerationData(20);
    Eigen::VectorXd referenceInvAngAccData(20);

    // Fill reference value arrays
    referenceAzimuthAngleData <<  0.0,	       0.7068583,	1.4137167,	2.120575,	2.8274334,	3.5342917,	4.2411501,	4.9480084,	5.6548668,	6.3617251,	7.0685835,	7.7754418,	8.4823002,	9.1891585,	9.8960169,	10.6028752,	11.3097336,	12.0165919,	12.7234502,	13.4303086;
    referenceRadiusData       << 1.4960E+08,   3.2520E+06,	1.0223E+05,	4.7032E+03,	3.2008E+02,	3.2526E+01,	4.9742E+00,	1.1524E+00,	4.0647E-01,	2.1908E-01,	1.8083E-01,	2.2870E-01,	4.4288E-01,	1.3101E+00,	5.8982E+00,	4.0201E+01,	4.1206E+02,	6.3010E+03,	1.4238E+05,	4.7035E+06;
    referenceAccelerationData << -2.8526E-06, -5.9799E-03,	-5.9707E+00,	-2.7646E+03,	-5.7791E+05,	-5.2856E+07,	-2.0069E+09,	-2.6650E+10,	5.3653E+10,	5.0006E+12,	-4.5171E+12,	-3.6749E+12,	-6.0218E+09,	2.1864E+10,	1.4520E+09,	3.4864E+07,	3.5013E+05,	1.5440E+03,	3.0827E+00,	2.8617E-03;
    referenceInvAngAccData    << 8.1862E+14,	6.9634E+09,	1.7431E+05,	1.3226E+01,	3.1110E-03,	2.2993E-06,	5.3392E-09,	3.8057E-11,	7.7509E-13,	3.7677E-14,	6.5250E-15,	4.9493E-14,	1.0982E-12,	5.9596E-11,	9.3408E-09,	4.5116E-06,	6.8521E-03,	3.2671E+01,	4.8202E+05,	2.1497E+10;

    // Initilise computed variables
    Eigen::Vector6d computedStateVector;
    double computedRadius;
    double computedAcceleration;
    double computedInvAngAcc;

    double differenceRadius;
    double differenceAcceleration;
    double differenceInvAngAcc;

    for (int idx = 0; idx < 20; idx++)
    {
        // Inputs
        referenceAzimuthAngle   = referenceAzimuthAngleData(idx); // [rad]
        // Computed variables
        computedStateVector     = exposinsShaping.computeCurrentStateVector(referenceAzimuthAngle);   // [m]
        computedRadius          = (1/1e3)*std::sqrt(std::pow(computedStateVector(0),2.0) + std::pow(computedStateVector(1),2.0)); // [m]
        computedAcceleration    = (1/1e3)*exposinsShaping.computeCurrentThrustAccelerationMagnitude(referenceAzimuthAngle);
        computedInvAngAcc       = 1.0/std::pow(exposinsShaping.computeFirstDerivativeAzimuthAngleWrtTime(referenceAzimuthAngle),2.0);
        // Reference values
        referenceRadius         = referenceRadiusData(idx);       // [km]
        referenceAcceleration   = referenceAccelerationData(idx); // [km/s^2]
        referenceInvAngAcc      = referenceInvAngAccData(idx);    // [(s/rad)^2]
        // Differences
        differenceRadius        = 100.0*std::fabs((referenceRadius - computedRadius)/referenceRadius);
        differenceAcceleration  = 100.0*std::fabs((referenceAcceleration - computedAcceleration)/referenceAcceleration);
        differenceInvAngAcc     = 100.0*std::fabs((referenceInvAngAcc - computedInvAngAcc)/referenceInvAngAcc);

        // Difference is smaller than 1.0%
        BOOST_CHECK_SMALL( differenceRadius, 1.0 );
        BOOST_CHECK_SMALL( differenceAcceleration, 1.0 );
        BOOST_CHECK_SMALL( differenceInvAngAcc, 1.0 );

        std::cout << "differenceRadius [%]:       " << differenceRadius << std::endl;
        std::cout << "differenceAcceleration [%]: " << differenceAcceleration << std::endl;
        std::cout << "differenceInvAngAcc [%]:    " << differenceInvAngAcc << std::endl;
    }

}

//! Test.
BOOST_AUTO_TEST_CASE( test_exposins_shaping_earth_mars_transfer )
{
    std::cout << "Running test_exposins_shaping_earth_mars_transfer ..." << std::endl;

    spice_interface::loadStandardSpiceKernels( );

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


    // Initialise input variables
    int numberOfRevolutions;
    double julianDate;
    double timeOfFlight;

    // Initialise reference variables
    double referenceDeltaV;
    double referenceTravelledAngle;
    double referenceThrustAcceleration;
    Eigen::Vector6d referenceState;
    Eigen::Vector3d referenceThrustAccelerationVector;

    // Initialise computed variables
    double computedTravelledAngle;
    double windingParameter;
    double computedDeltaV;
    double computedThrustAcceleration;
    Eigen::Vector6d computedState;
    Eigen::Vector3d computedThrustAccelerationVector;

    // Initialise input variable arrays
    Eigen::VectorXd julianDateInput(9);
    Eigen::VectorXd timeOfFlightInput(9);
    Eigen::VectorXd nRevolutionsInput(9);
    Eigen::VectorXd windingParameterInput(9);

    // Fill input variables
    julianDateInput       << 7355,  7805,   8255,   8055,   7605,   7305,   7805,   7805,   7805;
    timeOfFlightInput     << 460,   1580,   580,    900,    3740,   540,    1180,   1180,   1180;
    nRevolutionsInput     << 0,     0,      0,      1,      1,      1,      2,      2,      2;
    windingParameterInput << 0.1,   0.2,    0.3,    0.4,    0.5,    0.6,    0.7,    0.8,    0.9;

    // Initialise computed and reference arrays
    Eigen::VectorXd deltaVOutput(9);
    Eigen::VectorXd travelledAzimuthAngleOutput(9);
    Eigen::MatrixXd stateVectorOutput(9,6);
    Eigen::MatrixXd thrustAccVectorOutput(9,3);
    Eigen::VectorXd thrustAccelerationOutput(9);

    Eigen::MatrixXd computedStateMatrix(9,6);
    Eigen::MatrixXd computedThrustAccMatrix(9,3);

    // Fill reference arrays
    stateVectorOutput <<   -87300823973.5873,  87612058918.16795,  3981727047.135072, -9455.313138782409, -6786.388994149452, -1219.052523675202,
                           -85232192822.10521, -80484352414.51903,  414525754.2481378,  10562.26296679809, -3538.278292867627, -190.9065221799179,
                            51196907479.03938, -91777253831.16023, -3179710538.490304,  10118.20671591858,  6854.056030833378,  1279.576316681466,
                            101917580218.7555,  26225206844.57919, -1960261136.034621, -3887.577251797896,   11786.2833112428, -49.94921647350586,
                           -123560888959.3484,  7312465806.804946,  3192879569.775778,  4233.452739990917, -11111.27386869442, -488.0123143774545,
                           -105622947240.2362,  65765955417.79555,   3975956914.93759, -5468.189235414906, -9427.096820905601,  305.0730074055093,
                            90650985892.17734,  57589533410.87675, -1026950952.557124, -7071.668570220928,  9528.421316874243,  284.4138885149502,
                            90650985892.17729,  57589533410.87673, -1026950952.557123, -9914.244376356695,  3439.616503759996,  254.2088896828113,
                            90650985892.17735,  57589533410.87677, -1026950952.557124, -4299.252411750291,  10525.97068547897,  241.8777358380139;

    thrustAccVectorOutput << -6.608221701128012e-005, -4.742937898004895e-005, -8.519833477247171e-006,
                              -0.0002882466515458332,  9.656045047944662e-005,  5.209884089195635e-006,
                             -5.270547076116196e-005,  -3.57025963069246e-005, -6.665279138785385e-006,
                              1.358301213459705e-005, -4.118071973087375e-005,  1.745202138837179e-007,
                             -5.574402173476506e-005,    0.000146307784704066,  6.425905928395506e-006,
                              6.507537743312651e-006,  1.121892197413155e-005, -3.630587795498295e-007,
                              2.075750672603594e-005,  -2.79688262550096e-005, -8.348416141401088e-007,
                               0.0003895018199948488,  -0.0001351325262360707, -9.987127755944108e-006,
                             -2.692386193527443e-005,  6.591838634457043e-005,  1.514747714537705e-006;

    deltaVOutput << 6473.21590989085, 23574.32746479469, 6243.649295691202, 5961.544126065006, 21770.32507182609, 6811.804877402086, 9617.590709084887, 817409.8413658796, 15715.00973600247;
    travelledAzimuthAngleOutput << 5.999751984480445, 6.087112330499233, 5.994155141968389, 10.72204072408985, 8.767525609220572, 7.118626322681801, 15.32101794380945, 15.32101794380945, 15.32101794380945;
    thrustAccelerationOutput << 8.178626387138032e-005, -0.0003040348592125605, -6.400756187521011e-005, -4.336335265528301e-005, -0.0001566992536395048, -1.297474737402231e-005, -3.484000987492744e-005, -0.0004123981209305584, 7.122094116494437e-005;

    // Main loop - compute shape-based trajectory per input set
    for (int idx = 0; idx < 9; ++idx)
    {
        // Extract input variables
        julianDate          = julianDateInput(idx) * physical_constants::JULIAN_DAY;
        timeOfFlight        = timeOfFlightInput(idx);
        numberOfRevolutions = nRevolutionsInput(idx);
        windingParameter    = windingParameterInput(idx);

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

        // Define initial and final state
        Eigen::Vector6d initialState = pointerToDepartureBodyEphemeris->getCartesianState( julianDate );
        Eigen::Vector6d finalState   = pointerToArrivalBodyEphemeris->getCartesianState( julianDate + timeOfFlight * tudat::physical_constants::JULIAN_DAY );

        // Compute shaped trajectory.
        shape_based_methods::ExposinsShaping exposinsShaping = shape_based_methods::ExposinsShaping(
                    initialState, finalState, timeOfFlight * tudat::physical_constants::JULIAN_DAY,
                    numberOfRevolutions, bodyMap, "Vehicle", "Sun", windingParameter,
                    rootFinderSettings, integratorSettings );

        // Compute performance parameters
        computedTravelledAngle = exposinsShaping.getTravelledAzimuthAngle();
        computedState          = exposinsShaping.computeCurrentStateVector(computedTravelledAngle)/2.0;
        computedDeltaV         = exposinsShaping.computeDeltaV();
        computedThrustAcceleration       = exposinsShaping.computeCurrentThrustAccelerationMagnitude(computedTravelledAngle)/2.0;
        computedThrustAccelerationVector = exposinsShaping.computeCurrentThrustAccelerationVector(computedTravelledAngle)/2.0;

        if (exposinsShaping.getInfeasibleTOF()==false)
        {
            for (int idx2 = 0; idx2 < 6; ++idx2)
            {
                computedStateMatrix(idx,idx2) = computedState(idx2);                
            }
            for (int idx2 = 0; idx2 < 3; ++idx2)
            {
                computedThrustAccMatrix(idx,idx2) = computedThrustAccelerationVector(idx2);
            }

            // Check difference reference and computed values
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION( computedState, referenceState, 0.1 ); // m,m/s
            BOOST_CHECK_SMALL( std::fabs(  computedDeltaV - referenceDeltaV ), 0.1 );  // m/s
            BOOST_CHECK_SMALL( std::fabs(  computedTravelledAngle - referenceTravelledAngle ), 0.1 );  // rad
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION( computedThrustAccelerationVector, referenceThrustAccelerationVector, 0.1 ); // m,m/s
            BOOST_CHECK_SMALL( std::fabs(  computedThrustAcceleration - referenceThrustAcceleration ), 0.1 );  // rad

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
