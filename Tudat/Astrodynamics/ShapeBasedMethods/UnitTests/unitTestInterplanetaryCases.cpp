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

#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"
#include "Tudat/Mathematics/BasicMathematics/basicMathematicsFunctions.h"
#include "Tudat/Mathematics/BasicMathematics/coordinateConversions.h"
#include "Tudat/Astrodynamics/Ephemerides/approximatePlanetPositions.h"
#include "Tudat/SimulationSetup/tudatSimulationHeader.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/celestialBodyConstants.h"

#include "Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/keplerPropagator.h"




#include "Tudat/Astrodynamics/ShapeBasedMethods/baseFunctionsExposinsShaping.h"
#include "Tudat/Astrodynamics/ShapeBasedMethods/compositeFunctionExposinsShaping.h"
#include "Tudat/Astrodynamics/ShapeBasedMethods/exposinsShaping.h"

#include "Tudat/Astrodynamics/ShapeBasedMethods/baseFunctionsInvPolyShaping.h"
#include "Tudat/Astrodynamics/ShapeBasedMethods/compositeFunctionInvPolyShaping.h"
#include "Tudat/Astrodynamics/ShapeBasedMethods/invPolyShaping.h"

#include "Tudat/Astrodynamics/ShapeBasedMethods/baseFunctionsSphericalShaping.h"
#include "Tudat/Astrodynamics/ShapeBasedMethods/compositeFunctionSphericalShaping.h"
#include "Tudat/Astrodynamics/ShapeBasedMethods/sphericalShaping.h"

#include "Tudat/Astrodynamics/ShapeBasedMethods/compositeFunctionHodographicShaping.h"
#include "Tudat/Astrodynamics/ShapeBasedMethods/hodographicShaping.h"
#include "Tudat/Astrodynamics/ShapeBasedMethods/baseFunctionsHodographicShaping.h"
#include "Tudat/Astrodynamics/ShapeBasedMethods/createBaseFunctionHodographicShaping.h"

namespace tudat
{
namespace unit_tests
{

//! Test spherical shaping implementation.
BOOST_AUTO_TEST_SUITE( test_interplanetary_cases )


////! Earth Tempel-1.
//BOOST_AUTO_TEST_CASE( test_interplanetary_cases_earth_tempel_transfer )
//{
//    spice_interface::loadStandardSpiceKernels( );

//    int numberOfRevolutions = 1;
//    double julianDate = (8174.5 - 15.0*365.25) * physical_constants::JULIAN_DAY;
//    double  timeOfFlight = 580.0;

//    // Create central, departure and arrival bodies.
//    std::vector< std::string > bodiesToCreate;
//    bodiesToCreate.push_back( "Sun" );

//    std::map< std::string, std::shared_ptr< simulation_setup::BodySettings > > bodySettings =
//            simulation_setup::getDefaultBodySettings( bodiesToCreate );

//    std::string frameOrigin = "SSB";
//    std::string frameOrientation = "ECLIPJ2000";

//    // Define central body ephemeris settings.
//    bodySettings[ "Sun" ]->ephemerisSettings = std::make_shared< simulation_setup::ConstantEphemerisSettings >(
//                ( Eigen::Vector6d( ) << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ).finished( ), frameOrigin, frameOrientation );

//    bodySettings[ "Sun" ]->ephemerisSettings->resetFrameOrientation( frameOrientation );
//    bodySettings[ "Sun" ]->rotationModelSettings->resetOriginalFrame( frameOrientation );

//    // Create body map.
//    simulation_setup::NamedBodyMap bodyMap = createBodies( bodySettings );

//    bodyMap[ "Vehicle" ] = std::make_shared< simulation_setup::Body >( );
//    bodyMap.at( "Vehicle" )->setEphemeris( std::make_shared< ephemerides::TabulatedCartesianEphemeris< > >(
//                                                         std::shared_ptr< interpolators::OneDimensionalInterpolator
//                                                         < double, Eigen::Vector6d > >( ), frameOrigin, frameOrientation ) );

//    setGlobalFrameBodyEphemerides( bodyMap, frameOrigin, frameOrientation );

//    // Set vehicle mass.
//    bodyMap[ "Vehicle" ]->setConstantBodyMass( 400.0 );


//    // Ephemeris departure body.
//    ephemerides::EphemerisPointer pointerToDepartureBodyEphemeris = std::make_shared< ephemerides::ApproximatePlanetPositions>(
//                ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::earthMoonBarycenter );

//    // Ephemeris arrival body.
//    ephemerides::EphemerisPointer pointerToArrivalBodyEphemeris = std::make_shared< ephemerides::ApproximatePlanetPositions >(
//                ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::mars );


//    // Define root finder settings (used to update the updated value of the free coefficient, so that it matches the required time of flight).
//        std::shared_ptr< root_finders::RootFinderSettings > rootFinderSettings =
//                std::make_shared< root_finders::RootFinderSettings >( root_finders::bisection_root_finder, 1.0e-6, 150 );




//        // /////////////////////////// HODOGRAPHIC - START ////////////////////////////
//        // Initialize free coefficients vector for radial velocity function.
//        Eigen::VectorXd freeCoefficientsRadialVelocityFunction = Eigen::VectorXd::Zero( 0 );

//        // Initialize free coefficients vector for normal velocity function.
//        Eigen::VectorXd freeCoefficientsNormalVelocityFunction = Eigen::VectorXd::Zero( 0 );

//        // Initialize free coefficients vector for axial velocity function.
//        Eigen::VectorXd freeCoefficientsAxialVelocityFunction = Eigen::VectorXd::Zero( 0 );

//        // /////////////////////////// HODOGRAPHIC - END ////////////////////////////

//        // /////////////////////////// TEMPEL1 - START ////////////////////////////

//        std::string kernelPath = input_output::getSpiceKernelPath( );
//        spice_interface::clearSpiceKernels();
//        spice_interface::loadSpiceKernelInTudat( kernelPath + "tempel1_2000_2020.bsp" );
//    //    spice_interface::loadStandardSpiceKernels();
//        frameOrigin      = "SUN";


//        // /////////////////////////// TEMPEL1 - END ////////////////////////////




//    // general variables
//        // pointers of bodies
//        // ranges of startdate
//        // ranges of time of flight
//        // Nrev


   
//        double centralBodyGravitationalParameter = bodyMap["Sun"]->getGravityFieldModel()->getGravitationalParameter( );


   

//    // ranging values - EarthTempel
//    double startDateLower =  (20.0*365.25) * physical_constants::JULIAN_DAY; // MJD2000 -> 2020 //20
////    double startDateLower = (57632 - 51544.5) * physical_constants::JULIAN_DAY;// test
////    double startDateUpper = startDateLower + (2.*445.0 * physical_constants::JULIAN_DAY); // -synodic period earth Tempel however ecc
//    double startDateUpper = startDateLower + (5.0*365.25 * physical_constants::JULIAN_DAY); // -synodic period earth Tempel however ecc
//    double startDateStep  = 40.0 * physical_constants::JULIAN_DAY; // 40

//    double TOFLower = 100.0; // days / orbit is 2038 days
//    double TOFUpper = 8900.0; // days //8900
//    double TOFStep  = 100.0;

////    // ranging values - EarthMars
////    double startDateLower = ( 51545.0 - 51544.5 ) * physical_constants::JULIAN_DAY; // MJD to MJD2000
////    double startDateUpper = startDateLower + (1.*780.0 * physical_constants::JULIAN_DAY); // -synodic period earth Tempel however ecc
////    double startDateStep  = 15.0 * physical_constants::JULIAN_DAY;

////    double TOFLower = 10.0; // days / orbit is 2038 days
////    double TOFUpper = 2000.0; // days
////    double TOFStep  = 20.0;


//    int NrevLower = 0;
//    int NrevUpper = 2;

//    // create variables
//    double deltaVexposins;
//    double deltaVinvPoly;
//    double deltaVspherical;
//    double deltaVhodographic;

//    double peakAccelerationExposins    = 0.0;
//    double peakAccelerationInvPoly     = 0.0;
//    double peakAccelerationSpherical   = 0.0;
//    double peakAccelerationHodographic = 0.0;





//    // exposins
//    double windingParameter = 0.125; // set winding parameter


//    // Testing Tempel Location Data
//    std::ofstream fileTempelPos("dataTempelPosition.txt");
//    if (fileTempelPos.is_open())
//    {
//        for (int tempelDate = 0; tempelDate <= 7000; tempelDate++)
//        {
//            double templeDateSec = tempelDate*physical_constants::JULIAN_DAY;
//            Eigen::Vector6d finalState = spice_interface::getBodyCartesianStateAtEpoch( "TEMPEL_1", frameOrigin, frameOrientation, "NONE", templeDateSec);

//            Eigen::Vector6d keplerianState;
//                    keplerianState(0) = 3.145692355253531*physical_constants::ASTRONOMICAL_UNIT;
//                    keplerianState(1) = .5096307949367873;
//                    keplerianState(2) = 10.47382730765828*mathematical_constants::PI/180.0;
//                    keplerianState(3) = 179.203181927074	*mathematical_constants::PI/180.0;
//                    keplerianState(4) = 68.74970463832437*mathematical_constants::PI/180.0;
//            double meanMotion = .1766565332681526*mathematical_constants::PI/180.0/physical_constants::JULIAN_DAY;
//            double initialEpoch = (2457533.5 - 2400000.5 - 51544.0)*tudat::physical_constants::JULIAN_DAY; // JD -> MJD -> MJD2000
//            double initialMeanAnomaly = 347.7086591090512*mathematical_constants::PI/180.0;

//            double finalMeanAnomaly = initialMeanAnomaly + meanMotion*(templeDateSec - initialEpoch);
//            double eccentricAnomaly = orbital_element_conversions::convertMeanAnomalyToEccentricAnomaly( keplerianState(1), finalMeanAnomaly );
//            double trueAnomaly      = orbital_element_conversions::convertEccentricAnomalyToTrueAnomaly( eccentricAnomaly, keplerianState(1) );
//            keplerianState(5)       = trueAnomaly;

//            Eigen::Vector6d cartesianState = orbital_element_conversions::convertKeplerianToCartesianElements( keplerianState, centralBodyGravitationalParameter );


//            Eigen::Vector6d keplerianState_spice = orbital_element_conversions::convertCartesianToKeplerianElements( finalState, centralBodyGravitationalParameter );
//            Eigen::Vector6d keplerianState_kepl  = orbital_element_conversions::convertCartesianToKeplerianElements( cartesianState, centralBodyGravitationalParameter );

//            fileTempelPos << tempelDate << ' ' << 0 << ' ' << finalState[0] << ' ' << finalState[1] << ' ' << finalState[2] << ' ' << finalState[3] << ' ' << finalState[4] << ' ' << finalState[5] << ' ' << keplerianState_spice[0] << ' ' << keplerianState_spice[1] << ' ' << keplerianState_spice[2] << ' ' << keplerianState_spice[3] << ' ' << keplerianState_spice[4] << ' ' << keplerianState_spice[5] << '\n';
//            fileTempelPos << tempelDate << ' ' << 1 << ' ' << cartesianState[0] << ' ' << cartesianState[1] << ' ' << cartesianState[2] << ' ' << cartesianState[3] << ' ' << cartesianState[4] << ' ' << cartesianState[5] << ' ' << keplerianState_kepl[0] << ' ' << keplerianState_kepl[1] << ' ' << keplerianState_kepl[2] << ' ' << keplerianState_kepl[3] << ' ' << keplerianState_kepl[4] << ' ' << keplerianState_kepl[5] << '\n';


//        }
//    }


//    // Earth and Tempel-1 Location data per day for thesis presentation
//    std::ofstream filePlanetsPos("dataPlanetsPosition_02.txt");
//    if (filePlanetsPos.is_open())
//    {
////        for (int tempelDate = (startDateLower/physical_constants::JULIAN_DAY); tempelDate <= (startDateUpper/physical_constants::JULIAN_DAY)+TOFUpper; tempelDate++)
//        for (int tempelDate = 7200; tempelDate <= 13000; tempelDate++)

//        {
//            double templeDateSec = tempelDate*physical_constants::JULIAN_DAY;

//            Eigen::Vector6d keplerianState;
//                    keplerianState(0) = 3.145692355253531*physical_constants::ASTRONOMICAL_UNIT;
//                    keplerianState(1) = .5096307949367873;
//                    keplerianState(2) = 10.47382730765828*mathematical_constants::PI/180.0;
//                    keplerianState(3) = 179.203181927074	*mathematical_constants::PI/180.0;
//                    keplerianState(4) = 68.74970463832437*mathematical_constants::PI/180.0;
//            double meanMotion = .1766565332681526*mathematical_constants::PI/180.0/physical_constants::JULIAN_DAY;
//            double initialEpoch = (2457533.5 - 2400000.5 - 51544.0)*tudat::physical_constants::JULIAN_DAY; // JD -> MJD -> MJD2000
//            double initialMeanAnomaly = 347.7086591090512*mathematical_constants::PI/180.0;

//            double finalMeanAnomaly = initialMeanAnomaly + meanMotion*(templeDateSec - initialEpoch);
//            double eccentricAnomaly = orbital_element_conversions::convertMeanAnomalyToEccentricAnomaly( keplerianState(1), finalMeanAnomaly );
//            double trueAnomaly      = orbital_element_conversions::convertEccentricAnomalyToTrueAnomaly( eccentricAnomaly, keplerianState(1) );
//            keplerianState(5)       = trueAnomaly;


//            Eigen::Vector6d cartesianState1 = pointerToDepartureBodyEphemeris->getCartesianState( templeDateSec );
//            Eigen::Vector6d cartesianState2 = orbital_element_conversions::convertKeplerianToCartesianElements( keplerianState, centralBodyGravitationalParameter );


//            filePlanetsPos << physical_constants::JULIAN_DAY << ' ' << tempelDate << ' ' << cartesianState1[0] << ' ' << cartesianState1[1] << ' ' << cartesianState1[2] << ' ' << cartesianState1[3] << ' ' << cartesianState1[4] << ' ' << cartesianState1[5] << ' ' << cartesianState2[0] << ' ' << cartesianState2[1] << ' ' << cartesianState2[2] << ' ' << cartesianState2[3] << ' ' << cartesianState2[4] << ' ' << cartesianState2[5] << '\n';


//        }
//    }



//    int Nruns       = 1;
//    int NrunsFailed = 1;
//    // grid searches
////    std::ofstream file("dataGridSearchShapeMethodsEarthTempel11.txt");
////    std::ofstream file("dataGridSearchShapeMethodsEarthTempel1_run5.txt");
//    std::ofstream file("dataGridSearchShapeMethodsEarthTempel1_run6.txt");
////    std::ofstream file("dataGridSearchShapeMethodsEarthMars.txt");
////    std::ofstream file("dataGridSearchShapeMethodsEarthMarsV2.txt");

//    std::ofstream fileTrajectories("dataGridSearchTrajectoriesEarthTempel_run5.txt");

//bool runUnitTestGridSearchFlag = false;

//if (runUnitTestGridSearchFlag)
//    {
//    if (file.is_open())
//    {
//        for (int Nrev = NrevLower; Nrev <= NrevUpper; Nrev++)
//        {
//            for (double TOF = TOFLower; TOF <= TOFUpper; TOF = TOF + TOFStep)
//            {

//                // Define integrator settings.
//                std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings =
//                        std::make_shared< numerical_integrators::IntegratorSettings< double > > (
//                            numerical_integrators::rungeKutta4, 0.0, TOF * tudat::physical_constants::JULIAN_DAY / 500.0 );

//                // /////////////////////////// HODOGRAPHIC ////////////////////////////
//                double frequency = 2.0 * mathematical_constants::PI / ( TOF * tudat::physical_constants::JULIAN_DAY );

//                double scaleFactor = 1.0 / ( TOF * tudat::physical_constants::JULIAN_DAY );

//                // Create base function settings for the components of the radial velocity composite function.
//                std::shared_ptr< shape_based_methods::BaseFunctionHodographicShapingSettings > firstRadialVelocityBaseFunctionSettings =
//                        std::make_shared< shape_based_methods::BaseFunctionHodographicShapingSettings >( );
//                std::shared_ptr< shape_based_methods::BaseFunctionHodographicShapingSettings > secondRadialVelocityBaseFunctionSettings =
//                        std::make_shared< shape_based_methods::PowerFunctionHodographicShapingSettings >( 1.0, scaleFactor );
//                std::shared_ptr< shape_based_methods::BaseFunctionHodographicShapingSettings > thirdRadialVelocityBaseFunctionSettings =
//                        std::make_shared< shape_based_methods::TrigonometricFunctionHodographicShapingSettings >( frequency );

//                // Create components of the radial velocity composite function.
//                std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > > radialVelocityFunctionComponents;
//                radialVelocityFunctionComponents.push_back(
//                            createBaseFunctionHodographicShaping( shape_based_methods::constant, firstRadialVelocityBaseFunctionSettings ) );
//                radialVelocityFunctionComponents.push_back(
//                            createBaseFunctionHodographicShaping( shape_based_methods::scaledPower, secondRadialVelocityBaseFunctionSettings ) );
//                radialVelocityFunctionComponents.push_back(
//                            createBaseFunctionHodographicShaping( shape_based_methods::cosine, thirdRadialVelocityBaseFunctionSettings ) );

//                // Create base function settings for the components of the normal velocity composite function.
//                std::shared_ptr< shape_based_methods::BaseFunctionHodographicShapingSettings > firstNormalVelocityBaseFunctionSettings =
//                        std::make_shared< shape_based_methods::BaseFunctionHodographicShapingSettings >( );
//                std::shared_ptr< shape_based_methods::BaseFunctionHodographicShapingSettings > secondNormalVelocityBaseFunctionSettings =
//                        std::make_shared< shape_based_methods::PowerFunctionHodographicShapingSettings >( 1.0, scaleFactor );
//                std::shared_ptr< shape_based_methods::BaseFunctionHodographicShapingSettings > thirdNormalVelocityBaseFunctionSettings =
//                        std::make_shared< shape_based_methods::TrigonometricFunctionHodographicShapingSettings >( frequency );

//                // Create components of the normal velocity composite function.
//                std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > > normalVelocityFunctionComponents;
//                normalVelocityFunctionComponents.push_back(
//                            createBaseFunctionHodographicShaping( shape_based_methods::constant, firstNormalVelocityBaseFunctionSettings ) );
//                normalVelocityFunctionComponents.push_back(
//                            createBaseFunctionHodographicShaping( shape_based_methods::scaledPower, secondNormalVelocityBaseFunctionSettings ) );
//                normalVelocityFunctionComponents.push_back(
//                            createBaseFunctionHodographicShaping( shape_based_methods::cosine, thirdNormalVelocityBaseFunctionSettings ) );

//                // Create base function settings for the components of the axial velocity composite function.
//                std::shared_ptr< shape_based_methods::BaseFunctionHodographicShapingSettings > firstAxialVelocityBaseFunctionSettings =
//                        std::make_shared< shape_based_methods::TrigonometricFunctionHodographicShapingSettings >( ( Nrev + 0.5 ) * frequency );
//                std::shared_ptr< shape_based_methods::BaseFunctionHodographicShapingSettings > secondAxialVelocityBaseFunctionSettings =
//                        std::make_shared< shape_based_methods::PowerTimesTrigonometricFunctionHodographicShapingSettings >
//                        ( 3.0, ( Nrev + 0.5 ) * frequency, scaleFactor );
//                std::shared_ptr< shape_based_methods::BaseFunctionHodographicShapingSettings > thirdAxialVelocityBaseFunctionSettings =
//                        std::make_shared< shape_based_methods::PowerTimesTrigonometricFunctionHodographicShapingSettings >(
//                            3.0, ( Nrev + 0.5 ) * frequency, scaleFactor );

//                // Set components for the axial velocity function.
//                std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > > axialVelocityFunctionComponents;
//                axialVelocityFunctionComponents.push_back(
//                            createBaseFunctionHodographicShaping( shape_based_methods::cosine, firstAxialVelocityBaseFunctionSettings ) );
//                axialVelocityFunctionComponents.push_back(
//                            createBaseFunctionHodographicShaping( shape_based_methods::scaledPowerCosine, secondAxialVelocityBaseFunctionSettings ) );
//                axialVelocityFunctionComponents.push_back(
//                            createBaseFunctionHodographicShaping( shape_based_methods::scaledPowerSine, thirdAxialVelocityBaseFunctionSettings ) );

//                // /////////////////////////// HODOGRAPHIC ////////////////////////////

//                for (double startDate = startDateLower; startDate <= startDateUpper; startDate = startDate + startDateStep)
//                {
//                    std::cout << "-----------------------------------------------" << std::endl;
//                    std::cout << "LOOP - nRev: " << Nrev << " startDate: " << startDate/tudat::physical_constants::JULIAN_DAY << " TOF: " << TOF << std::endl;
//                    std::cout << "-----------------------------------------------" << std::endl;


////                    //        *                          keplerianElements( 0 ) = semiMajorAxis,                   [m] \n
////                    //        *                          keplerianElements( 1 ) = eccentricity,                    [-] \n
////                    //        *                          keplerianElements( 2 ) = inclination,                   [rad] \n
////                    //        *                          keplerianElements( 3 ) = argument of periapsis,         [rad] \n
////                    //        *                          keplerianElements( 4 ) = longitude of ascending node,   [rad] \n
////                    //        *                          keplerianElements( 5 ) = true anomaly.                  [rad] \n

//                    Eigen::Vector6d keplerianState;
//                            keplerianState(0) = 3.145692355253531*physical_constants::ASTRONOMICAL_UNIT;
//                            keplerianState(1) = .5096307949367873;
//                            keplerianState(2) = 10.47382730765828*mathematical_constants::PI/180.0;
//                            keplerianState(3) = 179.203181927074	*mathematical_constants::PI/180.0;
//                            keplerianState(4) = 68.74970463832437*mathematical_constants::PI/180.0;
//                    double meanMotion = .1766565332681526*mathematical_constants::PI/180.0/physical_constants::JULIAN_DAY;
//                    double initialEpoch = (2457533.5 - 2400000.5 - 51544.0)*tudat::physical_constants::JULIAN_DAY; // JD -> MJD -> MJD2000


//                    double initialMeanAnomaly = 347.7086591090512*mathematical_constants::PI/180.0;

//                    double tofEpoch         = TOF*tudat::physical_constants::JULIAN_DAY;
//                    double finalMeanAnomaly = initialMeanAnomaly + meanMotion*(startDate + tofEpoch - initialEpoch);
//                    double eccentricAnomaly = orbital_element_conversions::convertMeanAnomalyToEccentricAnomaly( keplerianState(1), finalMeanAnomaly );
//                    double trueAnomaly      = orbital_element_conversions::convertEccentricAnomalyToTrueAnomaly( eccentricAnomaly, keplerianState(1) );
//                    keplerianState(5)       = trueAnomaly;

//                    Eigen::Vector6d cartesianState = orbital_element_conversions::convertKeplerianToCartesianElements( keplerianState, centralBodyGravitationalParameter );




//                    // states
//                    // Earth
//                    Eigen::Vector6d initialState = pointerToDepartureBodyEphemeris->getCartesianState( startDate );
////                    // Mars
////                    Eigen::Vector6d finalState = pointerToArrivalBodyEphemeris->getCartesianState(
////                                startDate + TOF * tudat::physical_constants::JULIAN_DAY );
////                     Tempel 1
////                    Eigen::Vector6d finalState = spice_interface::getBodyCartesianStateAtEpoch( "TEMPEL_1", frameOrigin, frameOrientation, "NONE", startDate + TOF*tudat::physical_constants::JULIAN_DAY);
//                    Eigen::Vector6d finalState = cartesianState;


//                    std::cout << "MJD to Calendar (Start) " << basic_astrodynamics::convertJulianDayToCalendarDate(basic_astrodynamics::convertModifiedJulianDayToJulianDay(startDate/tudat::physical_constants::JULIAN_DAY + 51544.5)) << std::endl;
//                    std::cout << "MJD to Calendar (End) " << basic_astrodynamics::convertJulianDayToCalendarDate(basic_astrodynamics::convertModifiedJulianDayToJulianDay(startDate/tudat::physical_constants::JULIAN_DAY + 51544.5 + TOF)) << std::endl;
//                    std::cout << "MJD to Calendar (End) " << basic_astrodynamics::convertJulianDayToCalendarDate(spice_interface::convertEphemerisTimeToJulianDate( (startDate) + TOF*tudat::physical_constants::JULIAN_DAY)) << std::endl;
////                    std::cout << "keplerianState" << keplerianState << std::endl;
////                    std::cout << "finalState" << finalState << std::endl;
////                    std::cout << "cartesianState " << cartesianState << std::endl;
////                    Eigen::Vector6d diffState = cartesianState-finalState;
////                    Eigen::Vector3d diffPos   = diffState.segment(0, 3);
////                    Eigen::Vector3d diffVel   = diffState.segment(3, 3);

////                    std::cout << "diffState " << (diffState) << std::endl;
////                    std::cout << "diffPos [m] " << diffPos.norm() << std::endl;
////                    std::cout << "diffPos [%] " << diffPos.norm()/finalState.segment(0, 3).norm() << std::endl;
////                    std::cout << "diffVel [m/s] " << diffVel.norm() << std::endl;
////                    std::cout << "diffVel [%] " << diffVel.norm()/finalState.segment(3, 3).norm() << std::endl;



//                    // Compute shaped trajectory. - exposins
//                    shape_based_methods::ExposinsShaping exposinsShaping = shape_based_methods::ExposinsShaping(
//                                initialState, finalState, TOF * tudat::physical_constants::JULIAN_DAY,
//                                Nrev, bodyMap, "Vehicle", "Sun", windingParameter,
//                                rootFinderSettings, integratorSettings );

//                    // Compute shaped trajectory. - invpolyshaping
//                    shape_based_methods::InvPolyShaping invPolyShaping = shape_based_methods::InvPolyShaping(
//                                initialState, finalState, TOF * tudat::physical_constants::JULIAN_DAY,
//                                Nrev, bodyMap, "Vehicle", "Sun", 0.928,
//                                rootFinderSettings, -1.0e-3, 1.0, integratorSettings );


//                    bool sphericalShapingInfeasibleTOF;
//                    peakAccelerationSpherical   = 0.0;
//                    try
//                    {
//                    // Compute shaped trajectory. - spherical shaping - constants influence not yet adjusted
//                    shape_based_methods::SphericalShaping sphericalShaping = shape_based_methods::SphericalShaping(
//                                initialState, finalState, TOF * tudat::physical_constants::JULIAN_DAY,
//                                Nrev, bodyMap, "Vehicle", "Sun", 0.0,
//                                rootFinderSettings, -1.0e-1, 1.0e-1, integratorSettings );

//                    sphericalShapingInfeasibleTOF = false;

//                    deltaVspherical = sphericalShaping.computeDeltaV();
//                    std::cout << "DeltaV - spherical " <<  deltaVspherical << std::endl;


//                    if (fileTrajectories.is_open())
//                    {

//                    double currentTravelledAngle = sphericalShaping.getFinalAzimuthAngle() - sphericalShaping.getInitialAzimuthAngle();
//                    double stepSize = ( currentTravelledAngle ) / 100.0;

//                    Eigen::Vector6d currentStateVector;
//                    // Check that the trajectory is feasible, ie curved toward the central body.
//                        for ( int i = 0 ; i <= 100 ; i++ )
//                        {
//                            double currentThetaAngle = sphericalShaping.getInitialAzimuthAngle() + i * stepSize;

//                            currentStateVector = sphericalShaping.computeCurrentStateVector(currentThetaAngle);

//                            if ( sphericalShaping.computeCurrentThrustAccelerationVector( currentThetaAngle ).norm() >  peakAccelerationSpherical )
//                            {
//                                peakAccelerationSpherical = sphericalShaping.computeCurrentThrustAccelerationVector( currentThetaAngle ).norm();

//                            }

//                            fileTrajectories << 2 << ' ' << Nrev << ' ' << startDate/tudat::physical_constants::JULIAN_DAY << ' ' << TOF << ' ' << currentStateVector[0] << ' ' << currentStateVector[1] << ' ' << currentStateVector[2] << ' ' << currentStateVector[3] << ' ' << currentStateVector[4] << ' '<< currentStateVector[5] << ' ' << sphericalShaping.convertIndependentVariableToTime(currentThetaAngle)*tudat::physical_constants::JULIAN_YEAR << ' ' << sphericalShaping.computeCurrentThrustAccelerationVector( currentThetaAngle ).norm() << '\n';
//                        }
//                    }
//                    }
//                    catch (const std::runtime_error& error)
//                    {
//                        sphericalShapingInfeasibleTOF = true;
//                        deltaVspherical = 0.0;
//                    }


//                    // Create hodographic-shaping object with defined velocity functions and boundary conditions.
//                     shape_based_methods::HodographicShaping hodographicShaping = shape_based_methods::HodographicShaping(
//                                initialState, finalState,
//                                TOF * tudat::physical_constants::JULIAN_DAY, Nrev,
//                                bodyMap, "Vehicle", "Sun",
//                                radialVelocityFunctionComponents, normalVelocityFunctionComponents, axialVelocityFunctionComponents,
//                                freeCoefficientsRadialVelocityFunction, freeCoefficientsNormalVelocityFunction,
//                                freeCoefficientsAxialVelocityFunction, integratorSettings );

//                    deltaVhodographic = hodographicShaping.computeDeltaV();
//                    std::cout << "DeltaV - hodographic " <<  hodographicShaping.computeDeltaV() << std::endl;


//                    peakAccelerationHodographic = 0.0;
//                    // Trajectory
//                    if (fileTrajectories.is_open())
//                    {
//                    double stepSize = ( TOF * tudat::physical_constants::JULIAN_DAY ) / 100.0;
//                    Eigen::Vector6d currentStateVector;
//                    // Check that the trajectory is feasible, ie curved toward the central body.
//                        for ( int i = 0 ; i <= 100 ; i++ )
//                        {
//                            double currentTime = i * stepSize;

//                            currentStateVector = hodographicShaping.computeCurrentStateVector(currentTime);

//                            if ( hodographicShaping.computeCurrentThrustAccelerationVector( currentTime ).norm() >  peakAccelerationHodographic )
//                            {
//                                peakAccelerationHodographic = hodographicShaping.computeCurrentThrustAccelerationVector( currentTime ).norm();
//                            }


//                            fileTrajectories << 3 << ' ' << Nrev << ' ' << startDate/tudat::physical_constants::JULIAN_DAY << ' ' << TOF << ' ' << currentStateVector[0] << ' ' << currentStateVector[1] << ' ' << currentStateVector[2] << ' ' << currentStateVector[3] << ' ' << currentStateVector[4] << ' '<< currentStateVector[5] << ' ' << currentTime << ' ' << hodographicShaping.computeCurrentThrustAccelerationVector( currentTime ).norm() << '\n';
//                        }
//                    }





//                    bool infeasibleTOF = invPolyShaping.getInfeasibleTOF();
//                    peakAccelerationInvPoly = 0.0;
//                    if (infeasibleTOF)
//                    {
//                        deltaVinvPoly = 0.0;
//                        NrunsFailed++;
//                    }
//                    else
//                    {
//                        deltaVinvPoly = invPolyShaping.computeDeltaV();
//                        std::cout << "DeltaV - invpoly " <<  invPolyShaping.computeDeltaV() << std::endl;


//                        // Trajectory
//                        if (fileTrajectories.is_open())
//                        {

//                        double currentTravelledAngle = invPolyShaping.getTravelledAzimuthAngle();
//                        double stepSize = ( currentTravelledAngle ) / 100.0;
//                        Eigen::Vector6d currentStateVector;
//                        // Check that the trajectory is feasible, ie curved toward the central body.
//                            for ( int i = 0 ; i <= 100 ; i++ )
//                            {
//                                double currentThetaAngle = i * stepSize;

//                                currentStateVector = invPolyShaping.computeCurrentStateVector(currentThetaAngle);

//                                if ( (std::abs(invPolyShaping.computeCurrentThrustAccelerationMagnitude( currentThetaAngle ))  * physical_constants::ASTRONOMICAL_UNIT / std::pow(tudat::physical_constants::JULIAN_YEAR, 2.0) ) >  peakAccelerationInvPoly )
//                                {
//                                    // thrust is non-dimensional so dimensionalisation is needed
//                                    peakAccelerationInvPoly = std::abs(invPolyShaping.computeCurrentThrustAccelerationMagnitude( currentThetaAngle )) * physical_constants::ASTRONOMICAL_UNIT / std::pow(tudat::physical_constants::JULIAN_YEAR, 2.0);
//                                }

//                                fileTrajectories << 1 << ' ' << Nrev << ' ' << startDate/tudat::physical_constants::JULIAN_DAY << ' ' << TOF << ' ' << currentStateVector[0] << ' ' << currentStateVector[1] << ' ' << currentStateVector[2] << ' ' << currentStateVector[3] << ' ' << currentStateVector[4] << ' '<< currentStateVector[5] << ' ' << invPolyShaping.computeCurrentTimeFromAzimuthAngle(currentThetaAngle)*tudat::physical_constants::JULIAN_YEAR  << ' '  << std::abs(invPolyShaping.computeCurrentThrustAccelerationMagnitude( currentThetaAngle )) * physical_constants::ASTRONOMICAL_UNIT / std::pow(tudat::physical_constants::JULIAN_YEAR, 2.0) << '\n';
//                            }
//                        }

//                    }


//                    peakAccelerationExposins = 0.0;
//                    infeasibleTOF = exposinsShaping.getInfeasibleTOF();
//                    if (infeasibleTOF)
//                    {
//                        deltaVexposins = 0.0;
//                        NrunsFailed++;
//                    }
//                    else
//                    {
//                        std::cout << "DeltaV - exposins " <<  exposinsShaping.computeDeltaV() << std::endl;
//                        std::cout << "DeltaVB - exposins " << exposinsShaping.computeDeltaVBoundaries() << std::endl;

//                        deltaVexposins = exposinsShaping.computeDeltaV() + exposinsShaping.computeDeltaVBoundaries();

//                        // Trajectory
//                        if (fileTrajectories.is_open())
//                        {

//                        double currentTravelledAngle = exposinsShaping.getTravelledAzimuthAngle();
//                        double stepSize = ( currentTravelledAngle ) / 100.0;
//                        Eigen::Vector6d currentStateVector;
//                        // Check that the trajectory is feasible, ie curved toward the central body.
//                            for ( int i = 0 ; i <= 100 ; i++ )
//                            {
//                                double currentThetaAngle = i * stepSize;

//                                currentStateVector = exposinsShaping.computeCurrentStateVector(currentThetaAngle);

//                                if ( std::abs(exposinsShaping.computeCurrentNormalizedThrustAccelerationMagnitude( currentThetaAngle )) >  peakAccelerationExposins )
//                                {
//                                    peakAccelerationExposins = std::abs(exposinsShaping.computeCurrentNormalizedThrustAccelerationMagnitude( currentThetaAngle ));
//                                }

//                                fileTrajectories << 0 << ' ' << Nrev << ' ' << startDate/tudat::physical_constants::JULIAN_DAY << ' ' << TOF << ' ' << currentStateVector[0] << ' ' << currentStateVector[1] << ' ' << currentStateVector[2] << ' ' << currentStateVector[3] << ' ' << currentStateVector[4] << ' ' << currentStateVector[5] << ' ' << exposinsShaping.computeCurrentTimeFromAzimuthAngle(currentThetaAngle)*tudat::physical_constants::JULIAN_YEAR  << ' ' << std::abs(exposinsShaping.computeCurrentNormalizedThrustAccelerationMagnitude( currentThetaAngle )) << '\n';
//                            }
//                        }

//                    }




//                    file << Nrev << ' ' << startDate/tudat::physical_constants::JULIAN_DAY << ' ' << TOF << ' ' << deltaVexposins << ' ' << deltaVinvPoly << ' ' << deltaVspherical << ' ' << deltaVhodographic << ' ' << peakAccelerationExposins << ' ' << peakAccelerationInvPoly << ' ' << peakAccelerationSpherical << ' ' << peakAccelerationHodographic << '\n';

//                    Nruns++;
//                }
//            }

//        }
//    }
//    std::cout << "Nruns: "       << Nruns*4 << std::endl;
//    std::cout << "NrunsFailed: " << NrunsFailed << std::endl;



//        // exposins shaping
//            // deltaV
//            // maxT
//            // trajectories
//        // inverse polynomial shaping
//        // spherical shaping
//        // hodographic shaping


//    // optimisation
//        // exposins shaping
//        // inverse polynomial shaping
//        // spherical shaping
//        // hodographic shaping




//    BOOST_CHECK_SMALL( std::fabs(  0.0 ), 100.0 );

//}
//}



////! Earth Mercury.
//BOOST_AUTO_TEST_CASE( test_interplanetary_cases_earth_mercury_transfer )
//{
//    spice_interface::loadStandardSpiceKernels( );

//    int numberOfRevolutions = 1;
//    double julianDate = (8174.5 - 15.0*365.25) * physical_constants::JULIAN_DAY;
//    double  timeOfFlight = 580.0;

//    // Create central, departure and arrival bodies.
//    std::vector< std::string > bodiesToCreate;
//    bodiesToCreate.push_back( "Sun" );

//    std::map< std::string, std::shared_ptr< simulation_setup::BodySettings > > bodySettings =
//            simulation_setup::getDefaultBodySettings( bodiesToCreate );

//    std::string frameOrigin = "SSB";
//    std::string frameOrientation = "ECLIPJ2000";

//    // Define central body ephemeris settings.
//    bodySettings[ "Sun" ]->ephemerisSettings = std::make_shared< simulation_setup::ConstantEphemerisSettings >(
//                ( Eigen::Vector6d( ) << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ).finished( ), frameOrigin, frameOrientation );

//    bodySettings[ "Sun" ]->ephemerisSettings->resetFrameOrientation( frameOrientation );
//    bodySettings[ "Sun" ]->rotationModelSettings->resetOriginalFrame( frameOrientation );

//    // Create body map.
//    simulation_setup::NamedBodyMap bodyMap = createBodies( bodySettings );

//    bodyMap[ "Vehicle" ] = std::make_shared< simulation_setup::Body >( );
//    bodyMap.at( "Vehicle" )->setEphemeris( std::make_shared< ephemerides::TabulatedCartesianEphemeris< > >(
//                                                         std::shared_ptr< interpolators::OneDimensionalInterpolator
//                                                         < double, Eigen::Vector6d > >( ), frameOrigin, frameOrientation ) );

//    setGlobalFrameBodyEphemerides( bodyMap, frameOrigin, frameOrientation );

//    // Set vehicle mass.
//    bodyMap[ "Vehicle" ]->setConstantBodyMass( 400.0 );


//    // Ephemeris departure body.
//    ephemerides::EphemerisPointer pointerToDepartureBodyEphemeris = std::make_shared< ephemerides::ApproximatePlanetPositions>(
//                ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::earthMoonBarycenter );

//    // Ephemeris arrival body.
//    ephemerides::EphemerisPointer pointerToArrivalBodyEphemeris = std::make_shared< ephemerides::ApproximatePlanetPositions >(
//                ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::mercury );


//    // Define root finder settings (used to update the updated value of the free coefficient, so that it matches the required time of flight).
//        std::shared_ptr< root_finders::RootFinderSettings > rootFinderSettings =
//                std::make_shared< root_finders::RootFinderSettings >( root_finders::bisection_root_finder, 1.0e-6, 150 );




//        // /////////////////////////// HODOGRAPHIC - START ////////////////////////////
//        // Initialize free coefficients vector for radial velocity function.
//        Eigen::VectorXd freeCoefficientsRadialVelocityFunction = Eigen::VectorXd::Zero( 0 );

//        // Initialize free coefficients vector for normal velocity function.
//        Eigen::VectorXd freeCoefficientsNormalVelocityFunction = Eigen::VectorXd::Zero( 0 );

//        // Initialize free coefficients vector for axial velocity function.
//        Eigen::VectorXd freeCoefficientsAxialVelocityFunction = Eigen::VectorXd::Zero( 0 );

//        // /////////////////////////// HODOGRAPHIC - END ////////////////////////////


//    // ranging values - EarthMercury
//    double startDateLower =  (20.0*365.25) * physical_constants::JULIAN_DAY; // MJD2000 -> 2020
//    double startDateUpper = startDateLower + (3.0*365.25* physical_constants::JULIAN_DAY); // -synodic period earth Tempel however ecc
//    double startDateStep  = 20.0 * physical_constants::JULIAN_DAY;

//    double TOFLower = 60.0; // days / orbit is 2038 days
//    double TOFUpper = 1160.0; // days
//    double TOFStep  = 20.0;

//    int NrevLower = 0;
//    int NrevUpper = 2;

//    // create variables
//    double deltaVexposins;
//    double deltaVinvPoly;
//    double deltaVspherical;
//    double deltaVhodographic;

//    double peakAccelerationExposins    = 0.0;
//    double peakAccelerationInvPoly     = 0.0;
//    double peakAccelerationSpherical   = 0.0;
//    double peakAccelerationHodographic = 0.0;



//    // exposins
//    double windingParameter = 0.125; // set winding parameter


//    int Nruns       = 1;
//    int NrunsFailed = 1;


//    // grid searches
//    std::ofstream file("dataGridSearchShapeMethodsEarthMercury_run1.txt");
//    std::ofstream fileTrajectories("dataGridSearchTrajectoriesEarthMercury_run1.txt");


//    {
//    if (file.is_open())
//    {
//        for (int Nrev = NrevLower; Nrev <= NrevUpper; Nrev++)
//        {
//            for (double TOF = TOFLower; TOF <= TOFUpper; TOF = TOF + TOFStep)
//            {

//                // Define integrator settings.
//                std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings =
//                        std::make_shared< numerical_integrators::IntegratorSettings< double > > (
//                            numerical_integrators::rungeKutta4, 0.0, TOF * tudat::physical_constants::JULIAN_DAY / 500.0 );

//                // /////////////////////////// HODOGRAPHIC ////////////////////////////
//                double frequency = 2.0 * mathematical_constants::PI / ( TOF * tudat::physical_constants::JULIAN_DAY );

//                double scaleFactor = 1.0 / ( TOF * tudat::physical_constants::JULIAN_DAY );

//                // Create base function settings for the components of the radial velocity composite function.
//                std::shared_ptr< shape_based_methods::BaseFunctionHodographicShapingSettings > firstRadialVelocityBaseFunctionSettings =
//                        std::make_shared< shape_based_methods::BaseFunctionHodographicShapingSettings >( );
//                std::shared_ptr< shape_based_methods::BaseFunctionHodographicShapingSettings > secondRadialVelocityBaseFunctionSettings =
//                        std::make_shared< shape_based_methods::PowerFunctionHodographicShapingSettings >( 1.0, scaleFactor );
//                std::shared_ptr< shape_based_methods::BaseFunctionHodographicShapingSettings > thirdRadialVelocityBaseFunctionSettings =
//                        std::make_shared< shape_based_methods::TrigonometricFunctionHodographicShapingSettings >( frequency );

//                // Create components of the radial velocity composite function.
//                std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > > radialVelocityFunctionComponents;
//                radialVelocityFunctionComponents.push_back(
//                            createBaseFunctionHodographicShaping( shape_based_methods::constant, firstRadialVelocityBaseFunctionSettings ) );
//                radialVelocityFunctionComponents.push_back(
//                            createBaseFunctionHodographicShaping( shape_based_methods::scaledPower, secondRadialVelocityBaseFunctionSettings ) );
//                radialVelocityFunctionComponents.push_back(
//                            createBaseFunctionHodographicShaping( shape_based_methods::cosine, thirdRadialVelocityBaseFunctionSettings ) );

//                // Create base function settings for the components of the normal velocity composite function.
//                std::shared_ptr< shape_based_methods::BaseFunctionHodographicShapingSettings > firstNormalVelocityBaseFunctionSettings =
//                        std::make_shared< shape_based_methods::BaseFunctionHodographicShapingSettings >( );
//                std::shared_ptr< shape_based_methods::BaseFunctionHodographicShapingSettings > secondNormalVelocityBaseFunctionSettings =
//                        std::make_shared< shape_based_methods::PowerFunctionHodographicShapingSettings >( 1.0, scaleFactor );
//                std::shared_ptr< shape_based_methods::BaseFunctionHodographicShapingSettings > thirdNormalVelocityBaseFunctionSettings =
//                        std::make_shared< shape_based_methods::TrigonometricFunctionHodographicShapingSettings >( frequency );

//                // Create components of the normal velocity composite function.
//                std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > > normalVelocityFunctionComponents;
//                normalVelocityFunctionComponents.push_back(
//                            createBaseFunctionHodographicShaping( shape_based_methods::constant, firstNormalVelocityBaseFunctionSettings ) );
//                normalVelocityFunctionComponents.push_back(
//                            createBaseFunctionHodographicShaping( shape_based_methods::scaledPower, secondNormalVelocityBaseFunctionSettings ) );
//                normalVelocityFunctionComponents.push_back(
//                            createBaseFunctionHodographicShaping( shape_based_methods::cosine, thirdNormalVelocityBaseFunctionSettings ) );

//                // Create base function settings for the components of the axial velocity composite function.
//                std::shared_ptr< shape_based_methods::BaseFunctionHodographicShapingSettings > firstAxialVelocityBaseFunctionSettings =
//                        std::make_shared< shape_based_methods::TrigonometricFunctionHodographicShapingSettings >( ( Nrev + 0.5 ) * frequency );
//                std::shared_ptr< shape_based_methods::BaseFunctionHodographicShapingSettings > secondAxialVelocityBaseFunctionSettings =
//                        std::make_shared< shape_based_methods::PowerTimesTrigonometricFunctionHodographicShapingSettings >
//                        ( 3.0, ( Nrev + 0.5 ) * frequency, scaleFactor );
//                std::shared_ptr< shape_based_methods::BaseFunctionHodographicShapingSettings > thirdAxialVelocityBaseFunctionSettings =
//                        std::make_shared< shape_based_methods::PowerTimesTrigonometricFunctionHodographicShapingSettings >(
//                            3.0, ( Nrev + 0.5 ) * frequency, scaleFactor );

//                // Set components for the axial velocity function.
//                std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > > axialVelocityFunctionComponents;
//                axialVelocityFunctionComponents.push_back(
//                            createBaseFunctionHodographicShaping( shape_based_methods::cosine, firstAxialVelocityBaseFunctionSettings ) );
//                axialVelocityFunctionComponents.push_back(
//                            createBaseFunctionHodographicShaping( shape_based_methods::scaledPowerCosine, secondAxialVelocityBaseFunctionSettings ) );
//                axialVelocityFunctionComponents.push_back(
//                            createBaseFunctionHodographicShaping( shape_based_methods::scaledPowerSine, thirdAxialVelocityBaseFunctionSettings ) );

//                // /////////////////////////// HODOGRAPHIC ////////////////////////////

//                for (double startDate = startDateLower; startDate <= startDateUpper; startDate = startDate + startDateStep)
//                {
//                    std::cout << "-----------------------------------------------" << std::endl;
//                    std::cout << "LOOP - nRev: " << Nrev << " startDate: " << startDate/tudat::physical_constants::JULIAN_DAY << " TOF: " << TOF << std::endl;
//                    std::cout << "-----------------------------------------------" << std::endl;

//                    // states
//                    // Earth
//                    Eigen::Vector6d initialState = pointerToDepartureBodyEphemeris->getCartesianState( startDate );
//                    // Mercury
//                    Eigen::Vector6d finalState = pointerToArrivalBodyEphemeris->getCartesianState(startDate + TOF * tudat::physical_constants::JULIAN_DAY );


//                    std::cout << "MJD to Calendar (Start) " << basic_astrodynamics::convertJulianDayToCalendarDate(basic_astrodynamics::convertModifiedJulianDayToJulianDay(startDate/tudat::physical_constants::JULIAN_DAY + 51544.5)) << std::endl;
//                    std::cout << "MJD to Calendar (End) " << basic_astrodynamics::convertJulianDayToCalendarDate(basic_astrodynamics::convertModifiedJulianDayToJulianDay(startDate/tudat::physical_constants::JULIAN_DAY + 51544.5 + TOF)) << std::endl;
//                    std::cout << "MJD to Calendar (End) " << basic_astrodynamics::convertJulianDayToCalendarDate(spice_interface::convertEphemerisTimeToJulianDate( (startDate) + TOF*tudat::physical_constants::JULIAN_DAY)) << std::endl;


//                    // Compute shaped trajectory. - exposins
//                    shape_based_methods::ExposinsShaping exposinsShaping = shape_based_methods::ExposinsShaping(
//                                initialState, finalState, TOF * tudat::physical_constants::JULIAN_DAY,
//                                Nrev, bodyMap, "Vehicle", "Sun", windingParameter,
//                                rootFinderSettings, integratorSettings );

//                    // Compute shaped trajectory. - invpolyshaping
//                    shape_based_methods::InvPolyShaping invPolyShaping = shape_based_methods::InvPolyShaping(
//                                initialState, finalState, TOF * tudat::physical_constants::JULIAN_DAY,
//                                Nrev, bodyMap, "Vehicle", "Sun", 0.928,
//                                rootFinderSettings, -1.0e-3, 1.0, integratorSettings );


//                    peakAccelerationSpherical   = 0.0;
//                    try
//                    {
//                    // Compute shaped trajectory. - spherical shaping - constants influence not yet adjusted
//                    shape_based_methods::SphericalShaping sphericalShaping = shape_based_methods::SphericalShaping(
//                                initialState, finalState, TOF * tudat::physical_constants::JULIAN_DAY,
//                                Nrev, bodyMap, "Vehicle", "Sun", 0.0,
//                                rootFinderSettings, -1.0e-1, 1.0e-1, integratorSettings );

//                    deltaVspherical = sphericalShaping.computeDeltaV();
//                    std::cout << "DeltaV - spherical " <<  deltaVspherical << std::endl;

//                    if (fileTrajectories.is_open())
//                    {

//                    double currentTravelledAngle = sphericalShaping.getFinalAzimuthAngle() - sphericalShaping.getInitialAzimuthAngle();
//                    double stepSize = ( currentTravelledAngle ) / 100.0;

//                    Eigen::Vector6d currentStateVector;
//                    // Check that the trajectory is feasible, ie curved toward the central body.
//                        for ( int i = 0 ; i <= 100 ; i++ )
//                        {
//                            double currentThetaAngle = sphericalShaping.getInitialAzimuthAngle() + i * stepSize;

//                            currentStateVector = sphericalShaping.computeCurrentStateVector(currentThetaAngle);

//                            if ( sphericalShaping.computeCurrentThrustAccelerationVector( currentThetaAngle ).norm() >  peakAccelerationSpherical )
//                            {
//                                peakAccelerationSpherical = sphericalShaping.computeCurrentThrustAccelerationVector( currentThetaAngle ).norm();

//                            }

//                            fileTrajectories << 2 << ' ' << Nrev << ' ' << startDate/tudat::physical_constants::JULIAN_DAY << ' ' << TOF << ' ' << currentStateVector[0] << ' ' << currentStateVector[1] << ' ' << currentStateVector[2] << ' ' << currentStateVector[3] << ' ' << currentStateVector[4] << ' '<< currentStateVector[5] << ' ' << sphericalShaping.convertIndependentVariableToTime(currentThetaAngle)*tudat::physical_constants::JULIAN_YEAR << ' ' << sphericalShaping.computeCurrentThrustAccelerationVector( currentThetaAngle ).norm() << '\n';
//                        }
//                    }

//                    }
//                    catch (const std::runtime_error& error)
//                    {
//                        deltaVspherical = 0.0;
//                        NrunsFailed++;
//                    }


//                    // Create hodographic-shaping object with defined velocity functions and boundary conditions.
//                     shape_based_methods::HodographicShaping hodographicShaping = shape_based_methods::HodographicShaping(
//                                initialState, finalState,
//                                TOF * tudat::physical_constants::JULIAN_DAY, Nrev,
//                                bodyMap, "Vehicle", "Sun",
//                                radialVelocityFunctionComponents, normalVelocityFunctionComponents, axialVelocityFunctionComponents,
//                                freeCoefficientsRadialVelocityFunction, freeCoefficientsNormalVelocityFunction,
//                                freeCoefficientsAxialVelocityFunction, integratorSettings );

//                    deltaVhodographic = hodographicShaping.computeDeltaV();
//                    std::cout << "DeltaV - hodographic " <<  hodographicShaping.computeDeltaV() << std::endl;


//                    peakAccelerationHodographic = 0.0;
//                    // Trajectory
//                    if (fileTrajectories.is_open())
//                    {
//                    double stepSize = ( TOF * tudat::physical_constants::JULIAN_DAY ) / 100.0;
//                    Eigen::Vector6d currentStateVector;
//                    // Check that the trajectory is feasible, ie curved toward the central body.
//                        for ( int i = 0 ; i <= 100 ; i++ )
//                        {
//                            double currentTime = i * stepSize;

//                            currentStateVector = hodographicShaping.computeCurrentStateVector(currentTime);

//                            if ( hodographicShaping.computeCurrentThrustAccelerationVector( currentTime ).norm() >  peakAccelerationHodographic )
//                            {
//                                peakAccelerationHodographic = hodographicShaping.computeCurrentThrustAccelerationVector( currentTime ).norm();
//                            }


//                            fileTrajectories << 3 << ' ' << Nrev << ' ' << startDate/tudat::physical_constants::JULIAN_DAY << ' ' << TOF << ' ' << currentStateVector[0] << ' ' << currentStateVector[1] << ' ' << currentStateVector[2] << ' ' << currentStateVector[3] << ' ' << currentStateVector[4] << ' '<< currentStateVector[5] << ' ' << currentTime << ' ' << hodographicShaping.computeCurrentThrustAccelerationVector( currentTime ).norm() << '\n';
//                        }
//                    }





//                    bool infeasibleTOF = invPolyShaping.getInfeasibleTOF();
//                    peakAccelerationInvPoly = 0.0;
//                    if (infeasibleTOF)
//                    {
//                        deltaVinvPoly = 0.0;
//                        NrunsFailed++;
//                    }
//                    else
//                    {
//                        deltaVinvPoly = invPolyShaping.computeDeltaV();
//                        std::cout << "DeltaV - invpoly " <<  invPolyShaping.computeDeltaV() << std::endl;


//                        // Trajectory
//                        if (fileTrajectories.is_open())
//                        {

//                        double currentTravelledAngle = invPolyShaping.getTravelledAzimuthAngle();
//                        double stepSize = ( currentTravelledAngle ) / 100.0;
//                        Eigen::Vector6d currentStateVector;
//                        // Check that the trajectory is feasible, ie curved toward the central body.
//                            for ( int i = 0 ; i <= 100 ; i++ )
//                            {
//                                double currentThetaAngle = i * stepSize;

//                                currentStateVector = invPolyShaping.computeCurrentStateVector(currentThetaAngle);

//                                if ( (std::abs(invPolyShaping.computeCurrentThrustAccelerationMagnitude( currentThetaAngle ))  * physical_constants::ASTRONOMICAL_UNIT / std::pow(tudat::physical_constants::JULIAN_YEAR, 2.0) ) >  peakAccelerationInvPoly )
//                                {
//                                    // thrust is non-dimensional so dimensionalisation is needed
//                                    peakAccelerationInvPoly = std::abs(invPolyShaping.computeCurrentThrustAccelerationMagnitude( currentThetaAngle )) * physical_constants::ASTRONOMICAL_UNIT / std::pow(tudat::physical_constants::JULIAN_YEAR, 2.0);
//                                }

//                                fileTrajectories << 1 << ' ' << Nrev << ' ' << startDate/tudat::physical_constants::JULIAN_DAY << ' ' << TOF << ' ' << currentStateVector[0] << ' ' << currentStateVector[1] << ' ' << currentStateVector[2] << ' ' << currentStateVector[3] << ' ' << currentStateVector[4] << ' '<< currentStateVector[5] << ' ' << invPolyShaping.computeCurrentTimeFromAzimuthAngle(currentThetaAngle)*tudat::physical_constants::JULIAN_YEAR  << ' '  << std::abs(invPolyShaping.computeCurrentThrustAccelerationMagnitude( currentThetaAngle )) * physical_constants::ASTRONOMICAL_UNIT / std::pow(tudat::physical_constants::JULIAN_YEAR, 2.0) << '\n';
//                            }
//                        }

//                    }


//                    peakAccelerationExposins = 0.0;
//                    infeasibleTOF = exposinsShaping.getInfeasibleTOF();
//                    if (infeasibleTOF)
//                    {
//                        deltaVexposins = 0.0;
//                        NrunsFailed++;
//                    }
//                    else
//                    {
//                        std::cout << "DeltaV - exposins " <<  exposinsShaping.computeDeltaV() << std::endl;
//                        std::cout << "DeltaVB - exposins " << exposinsShaping.computeDeltaVBoundaries() << std::endl;

//                        deltaVexposins = exposinsShaping.computeDeltaV() + exposinsShaping.computeDeltaVBoundaries();

//                        // Trajectory
//                        if (fileTrajectories.is_open())
//                        {

//                        double currentTravelledAngle = exposinsShaping.getTravelledAzimuthAngle();
//                        double stepSize = ( currentTravelledAngle ) / 100.0;
//                        Eigen::Vector6d currentStateVector;
//                        // Check that the trajectory is feasible, ie curved toward the central body.
//                            for ( int i = 0 ; i <= 100 ; i++ )
//                            {
//                                double currentThetaAngle = i * stepSize;

//                                currentStateVector = exposinsShaping.computeCurrentStateVector(currentThetaAngle);

//                                if ( std::abs(exposinsShaping.computeCurrentNormalizedThrustAccelerationMagnitude( currentThetaAngle )) >  peakAccelerationExposins )
//                                {
//                                    peakAccelerationExposins = std::abs(exposinsShaping.computeCurrentNormalizedThrustAccelerationMagnitude( currentThetaAngle ));
//                                }

//                                fileTrajectories << 0 << ' ' << Nrev << ' ' << startDate/tudat::physical_constants::JULIAN_DAY << ' ' << TOF << ' ' << currentStateVector[0] << ' ' << currentStateVector[1] << ' ' << currentStateVector[2] << ' ' << currentStateVector[3] << ' ' << currentStateVector[4] << ' ' << currentStateVector[5] << ' ' << exposinsShaping.computeCurrentTimeFromAzimuthAngle(currentThetaAngle)*tudat::physical_constants::JULIAN_YEAR  << ' ' << std::abs(exposinsShaping.computeCurrentNormalizedThrustAccelerationMagnitude( currentThetaAngle )) << '\n';
//                            }
//                        }

//                    }




//                    file << Nrev << ' ' << startDate/tudat::physical_constants::JULIAN_DAY << ' ' << TOF << ' ' << deltaVexposins << ' ' << deltaVinvPoly << ' ' << deltaVspherical << ' ' << deltaVhodographic << ' ' << peakAccelerationExposins << ' ' << peakAccelerationInvPoly << ' ' << peakAccelerationSpherical << ' ' << peakAccelerationHodographic << '\n';

//                    Nruns++;
//                }
//            }

//        }
//    }
//    std::cout << "Nruns: "       << Nruns*4 << std::endl;
//    std::cout << "NrunsFailed: " << NrunsFailed << std::endl;



//    BOOST_CHECK_SMALL( std::fabs(  0.0 ), 100.0 );

//}


//    }


////! Earth Mars.
//BOOST_AUTO_TEST_CASE( test_interplanetary_cases_earth_mars_transfer )
//{
//    spice_interface::loadStandardSpiceKernels( );

//    // Create central, departure and arrival bodies.
//    std::vector< std::string > bodiesToCreate;
//    bodiesToCreate.push_back( "Sun" );

//    std::map< std::string, std::shared_ptr< simulation_setup::BodySettings > > bodySettings =
//            simulation_setup::getDefaultBodySettings( bodiesToCreate );

//    std::string frameOrigin = "SSB";
//    std::string frameOrientation = "ECLIPJ2000";

//    // Define central body ephemeris settings.
//    bodySettings[ "Sun" ]->ephemerisSettings = std::make_shared< simulation_setup::ConstantEphemerisSettings >(
//                ( Eigen::Vector6d( ) << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ).finished( ), frameOrigin, frameOrientation );

//    bodySettings[ "Sun" ]->ephemerisSettings->resetFrameOrientation( frameOrientation );
//    bodySettings[ "Sun" ]->rotationModelSettings->resetOriginalFrame( frameOrientation );

//    // Create body map.
//    simulation_setup::NamedBodyMap bodyMap = createBodies( bodySettings );

//    bodyMap[ "Vehicle" ] = std::make_shared< simulation_setup::Body >( );
//    bodyMap.at( "Vehicle" )->setEphemeris( std::make_shared< ephemerides::TabulatedCartesianEphemeris< > >(
//                                                         std::shared_ptr< interpolators::OneDimensionalInterpolator
//                                                         < double, Eigen::Vector6d > >( ), frameOrigin, frameOrientation ) );

//    setGlobalFrameBodyEphemerides( bodyMap, frameOrigin, frameOrientation );

//    // Set vehicle mass.
//    bodyMap[ "Vehicle" ]->setConstantBodyMass( 400.0 );


//    // Ephemeris departure body.
//    ephemerides::EphemerisPointer pointerToDepartureBodyEphemeris = std::make_shared< ephemerides::ApproximatePlanetPositions>(
//                ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::earthMoonBarycenter );

//    // Ephemeris arrival body.
//    ephemerides::EphemerisPointer pointerToArrivalBodyEphemeris = std::make_shared< ephemerides::ApproximatePlanetPositions >(
//                ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::mars );


//    // Define root finder settings (used to update the updated value of the free coefficient, so that it matches the required time of flight).
//        std::shared_ptr< root_finders::RootFinderSettings > rootFinderSettings =
//                std::make_shared< root_finders::RootFinderSettings >( root_finders::bisection_root_finder, 1.0e-6, 150 );




//        // /////////////////////////// HODOGRAPHIC - START ////////////////////////////
//        // Initialize free coefficients vector for radial velocity function.
//        Eigen::VectorXd freeCoefficientsRadialVelocityFunction = Eigen::VectorXd::Zero( 0 );

//        // Initialize free coefficients vector for normal velocity function.
//        Eigen::VectorXd freeCoefficientsNormalVelocityFunction = Eigen::VectorXd::Zero( 0 );

//        // Initialize free coefficients vector for axial velocity function.
//        Eigen::VectorXd freeCoefficientsAxialVelocityFunction = Eigen::VectorXd::Zero( 0 );

//        // /////////////////////////// HODOGRAPHIC - END ////////////////////////////


//    // ranging values - EarthMars
//    double startDateLower =  (20.0*365.25) * physical_constants::JULIAN_DAY; // MJD2000 -> 2020
//    double startDateUpper = startDateLower + (3.0*365.25* physical_constants::JULIAN_DAY); // - //
//    double startDateStep  = 50.0 * physical_constants::JULIAN_DAY;

//    double TOFLower = 100.0; // days / orbit is 2038 days
//    double TOFUpper = 7800.0; // synodic period earth Mars however ecc  // 7800.0
//    double TOFStep  = 40.0;

//    int NrevLower = 0;
//    int NrevUpper = 2; // 2

//    // create variables
//    double deltaVexposins;
//    double deltaVinvPoly;
//    double deltaVspherical;
//    double deltaVhodographic;

//    double peakAccelerationExposins    = 0.0;
//    double peakAccelerationInvPoly     = 0.0;
//    double peakAccelerationSpherical   = 0.0;
//    double peakAccelerationHodographic = 0.0;



//    // exposins
//    double windingParameter = 0.125; // set winding parameter


//    int Nruns       = 1;
//    int NrunsFailed = 1;


//    // grid searches
//    std::ofstream file("dataGridSearchShapeMethodsEarthMars_run2.txt");
//    std::ofstream fileTrajectories("dataGridSearchTrajectoriesEarthMars_run2.txt");


//    {
//    if (file.is_open())
//    {
//        for (int Nrev = NrevLower; Nrev <= NrevUpper; Nrev++)
//        {
//            for (double TOF = TOFLower; TOF <= TOFUpper; TOF = TOF + TOFStep)
//            {

//                // Define integrator settings.
//                std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings =
//                        std::make_shared< numerical_integrators::IntegratorSettings< double > > (
//                            numerical_integrators::rungeKutta4, 0.0, TOF * tudat::physical_constants::JULIAN_DAY / 500.0 );

//                // /////////////////////////// HODOGRAPHIC ////////////////////////////
//                double frequency = 2.0 * mathematical_constants::PI / ( TOF * tudat::physical_constants::JULIAN_DAY );

//                double scaleFactor = 1.0 / ( TOF * tudat::physical_constants::JULIAN_DAY );

//                // Create base function settings for the components of the radial velocity composite function.
//                std::shared_ptr< shape_based_methods::BaseFunctionHodographicShapingSettings > firstRadialVelocityBaseFunctionSettings =
//                        std::make_shared< shape_based_methods::BaseFunctionHodographicShapingSettings >( );
//                std::shared_ptr< shape_based_methods::BaseFunctionHodographicShapingSettings > secondRadialVelocityBaseFunctionSettings =
//                        std::make_shared< shape_based_methods::PowerFunctionHodographicShapingSettings >( 1.0, scaleFactor );
//                std::shared_ptr< shape_based_methods::BaseFunctionHodographicShapingSettings > thirdRadialVelocityBaseFunctionSettings =
//                        std::make_shared< shape_based_methods::TrigonometricFunctionHodographicShapingSettings >( frequency );

//                // Create components of the radial velocity composite function.
//                std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > > radialVelocityFunctionComponents;
//                radialVelocityFunctionComponents.push_back(
//                            createBaseFunctionHodographicShaping( shape_based_methods::constant, firstRadialVelocityBaseFunctionSettings ) );
//                radialVelocityFunctionComponents.push_back(
//                            createBaseFunctionHodographicShaping( shape_based_methods::scaledPower, secondRadialVelocityBaseFunctionSettings ) );
//                radialVelocityFunctionComponents.push_back(
//                            createBaseFunctionHodographicShaping( shape_based_methods::cosine, thirdRadialVelocityBaseFunctionSettings ) );

//                // Create base function settings for the components of the normal velocity composite function.
//                std::shared_ptr< shape_based_methods::BaseFunctionHodographicShapingSettings > firstNormalVelocityBaseFunctionSettings =
//                        std::make_shared< shape_based_methods::BaseFunctionHodographicShapingSettings >( );
//                std::shared_ptr< shape_based_methods::BaseFunctionHodographicShapingSettings > secondNormalVelocityBaseFunctionSettings =
//                        std::make_shared< shape_based_methods::PowerFunctionHodographicShapingSettings >( 1.0, scaleFactor );
//                std::shared_ptr< shape_based_methods::BaseFunctionHodographicShapingSettings > thirdNormalVelocityBaseFunctionSettings =
//                        std::make_shared< shape_based_methods::TrigonometricFunctionHodographicShapingSettings >( frequency );

//                // Create components of the normal velocity composite function.
//                std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > > normalVelocityFunctionComponents;
//                normalVelocityFunctionComponents.push_back(
//                            createBaseFunctionHodographicShaping( shape_based_methods::constant, firstNormalVelocityBaseFunctionSettings ) );
//                normalVelocityFunctionComponents.push_back(
//                            createBaseFunctionHodographicShaping( shape_based_methods::scaledPower, secondNormalVelocityBaseFunctionSettings ) );
//                normalVelocityFunctionComponents.push_back(
//                            createBaseFunctionHodographicShaping( shape_based_methods::cosine, thirdNormalVelocityBaseFunctionSettings ) );

//                // Create base function settings for the components of the axial velocity composite function.
//                std::shared_ptr< shape_based_methods::BaseFunctionHodographicShapingSettings > firstAxialVelocityBaseFunctionSettings =
//                        std::make_shared< shape_based_methods::TrigonometricFunctionHodographicShapingSettings >( ( Nrev + 0.5 ) * frequency );
//                std::shared_ptr< shape_based_methods::BaseFunctionHodographicShapingSettings > secondAxialVelocityBaseFunctionSettings =
//                        std::make_shared< shape_based_methods::PowerTimesTrigonometricFunctionHodographicShapingSettings >
//                        ( 3.0, ( Nrev + 0.5 ) * frequency, scaleFactor );
//                std::shared_ptr< shape_based_methods::BaseFunctionHodographicShapingSettings > thirdAxialVelocityBaseFunctionSettings =
//                        std::make_shared< shape_based_methods::PowerTimesTrigonometricFunctionHodographicShapingSettings >(
//                            3.0, ( Nrev + 0.5 ) * frequency, scaleFactor );

//                // Set components for the axial velocity function.
//                std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > > axialVelocityFunctionComponents;
//                axialVelocityFunctionComponents.push_back(
//                            createBaseFunctionHodographicShaping( shape_based_methods::cosine, firstAxialVelocityBaseFunctionSettings ) );
//                axialVelocityFunctionComponents.push_back(
//                            createBaseFunctionHodographicShaping( shape_based_methods::scaledPowerCosine, secondAxialVelocityBaseFunctionSettings ) );
//                axialVelocityFunctionComponents.push_back(
//                            createBaseFunctionHodographicShaping( shape_based_methods::scaledPowerSine, thirdAxialVelocityBaseFunctionSettings ) );

//                // /////////////////////////// HODOGRAPHIC ////////////////////////////

//                for (double startDate = startDateLower; startDate <= startDateUpper; startDate = startDate + startDateStep)
//                {
//                    std::cout << "-----------------------------------------------" << std::endl;
//                    std::cout << "LOOP - nRev: " << Nrev << " startDate: " << startDate/tudat::physical_constants::JULIAN_DAY << " TOF: " << TOF << std::endl;
//                    std::cout << "-----------------------------------------------" << std::endl;

//                    // states
//                    // Earth
//                    Eigen::Vector6d initialState = pointerToDepartureBodyEphemeris->getCartesianState( startDate );
//                    // Mars
//                    Eigen::Vector6d finalState = pointerToArrivalBodyEphemeris->getCartesianState(startDate + TOF * tudat::physical_constants::JULIAN_DAY );


//                    std::cout << "MJD to Calendar (Start) " << basic_astrodynamics::convertJulianDayToCalendarDate(basic_astrodynamics::convertModifiedJulianDayToJulianDay(startDate/tudat::physical_constants::JULIAN_DAY + 51544.5)) << std::endl;
//                    std::cout << "MJD to Calendar (End) " << basic_astrodynamics::convertJulianDayToCalendarDate(basic_astrodynamics::convertModifiedJulianDayToJulianDay(startDate/tudat::physical_constants::JULIAN_DAY + 51544.5 + TOF)) << std::endl;
//                    std::cout << "MJD to Calendar (End) " << basic_astrodynamics::convertJulianDayToCalendarDate(spice_interface::convertEphemerisTimeToJulianDate( (startDate) + TOF*tudat::physical_constants::JULIAN_DAY)) << std::endl;


//                    // Compute shaped trajectory. - exposins
//                    shape_based_methods::ExposinsShaping exposinsShaping = shape_based_methods::ExposinsShaping(
//                                initialState, finalState, TOF * tudat::physical_constants::JULIAN_DAY,
//                                Nrev, bodyMap, "Vehicle", "Sun", windingParameter,
//                                rootFinderSettings, integratorSettings );

//                    // Compute shaped trajectory. - invpolyshaping
//                    shape_based_methods::InvPolyShaping invPolyShaping = shape_based_methods::InvPolyShaping(
//                                initialState, finalState, TOF * tudat::physical_constants::JULIAN_DAY,
//                                Nrev, bodyMap, "Vehicle", "Sun", 0.928,
//                                rootFinderSettings, -1.0e-3, 1.0, integratorSettings );


//                    bool sphericalShapingInfeasibleTOF;
//                    peakAccelerationSpherical   = 0.0;
//                    try
//                    {
//                    // Compute shaped trajectory. - spherical shaping - constants influence not yet adjusted
//                    shape_based_methods::SphericalShaping sphericalShaping = shape_based_methods::SphericalShaping(
//                                initialState, finalState, TOF * tudat::physical_constants::JULIAN_DAY,
//                                Nrev, bodyMap, "Vehicle", "Sun", 0.0,
//                                rootFinderSettings, -1.0e-1, 1.0e-1, integratorSettings );

//                    sphericalShapingInfeasibleTOF = false;

//                    deltaVspherical = sphericalShaping.computeDeltaV();
//                    std::cout << "DeltaV - spherical " <<  deltaVspherical << std::endl;


//                    if (fileTrajectories.is_open())
//                    {

//                    double currentTravelledAngle = sphericalShaping.getFinalAzimuthAngle() - sphericalShaping.getInitialAzimuthAngle();
//                    double stepSize = ( currentTravelledAngle ) / 100.0;

//                    Eigen::Vector6d currentStateVector;
//                    // Check that the trajectory is feasible, ie curved toward the central body.
//                        for ( int i = 0 ; i <= 100 ; i++ )
//                        {
//                            double currentThetaAngle = sphericalShaping.getInitialAzimuthAngle() + i * stepSize;

//                            currentStateVector = sphericalShaping.computeCurrentStateVector(currentThetaAngle);

//                            if ( sphericalShaping.computeCurrentThrustAccelerationVector( currentThetaAngle ).norm() >  peakAccelerationSpherical )
//                            {
//                                peakAccelerationSpherical = sphericalShaping.computeCurrentThrustAccelerationVector( currentThetaAngle ).norm();

//                            }

//                            fileTrajectories << 2 << ' ' << Nrev << ' ' << startDate/tudat::physical_constants::JULIAN_DAY << ' ' << TOF << ' ' << currentStateVector[0] << ' ' << currentStateVector[1] << ' ' << currentStateVector[2] << ' ' << currentStateVector[3] << ' ' << currentStateVector[4] << ' '<< currentStateVector[5] << ' ' << sphericalShaping.convertIndependentVariableToTime(currentThetaAngle)*tudat::physical_constants::JULIAN_YEAR << ' ' << sphericalShaping.computeCurrentThrustAccelerationVector( currentThetaAngle ).norm() << '\n';
//                        }
//                    }
//                    }
//                    catch (const std::runtime_error& error)
//                    {
//                        sphericalShapingInfeasibleTOF = true;
//                    }


//                    // Create hodographic-shaping object with defined velocity functions and boundary conditions.
//                     shape_based_methods::HodographicShaping hodographicShaping = shape_based_methods::HodographicShaping(
//                                initialState, finalState,
//                                TOF * tudat::physical_constants::JULIAN_DAY, Nrev,
//                                bodyMap, "Vehicle", "Sun",
//                                radialVelocityFunctionComponents, normalVelocityFunctionComponents, axialVelocityFunctionComponents,
//                                freeCoefficientsRadialVelocityFunction, freeCoefficientsNormalVelocityFunction,
//                                freeCoefficientsAxialVelocityFunction, integratorSettings );





//                    bool infeasibleTOF = invPolyShaping.getInfeasibleTOF();
//                    peakAccelerationInvPoly = 0.0;
//                    if (infeasibleTOF)
//                    {
//                        deltaVinvPoly = 0.0;
//                        NrunsFailed++;
//                    }
//                    else
//                    {
//                        deltaVinvPoly = invPolyShaping.computeDeltaV();
//                        std::cout << "DeltaV - invpoly " <<  invPolyShaping.computeDeltaV() << std::endl;


//                        // Trajectory
//                        if (fileTrajectories.is_open())
//                        {

//                        double currentTravelledAngle = invPolyShaping.getTravelledAzimuthAngle();
//                        double stepSize = ( currentTravelledAngle ) / 100.0;
//                        Eigen::Vector6d currentStateVector;
//                        // Check that the trajectory is feasible, ie curved toward the central body.
//                            for ( int i = 0 ; i <= 100 ; i++ )
//                            {
//                                double currentThetaAngle = i * stepSize;

//                                currentStateVector = invPolyShaping.computeCurrentStateVector(currentThetaAngle);

//                                if ( (std::abs(invPolyShaping.computeCurrentThrustAccelerationMagnitude( currentThetaAngle ))  * physical_constants::ASTRONOMICAL_UNIT / std::pow(tudat::physical_constants::JULIAN_YEAR, 2.0) ) >  peakAccelerationInvPoly )
//                                {
//                                    // thrust is non-dimensional so dimensionalisation is needed
//                                    peakAccelerationInvPoly = std::abs(invPolyShaping.computeCurrentThrustAccelerationMagnitude( currentThetaAngle )) * physical_constants::ASTRONOMICAL_UNIT / std::pow(tudat::physical_constants::JULIAN_YEAR, 2.0);
//                                }

//                                fileTrajectories << 1 << ' ' << Nrev << ' ' << startDate/tudat::physical_constants::JULIAN_DAY << ' ' << TOF << ' ' << currentStateVector[0] << ' ' << currentStateVector[1] << ' ' << currentStateVector[2] << ' ' << currentStateVector[3] << ' ' << currentStateVector[4] << ' '<< currentStateVector[5] << ' ' << invPolyShaping.computeCurrentTimeFromAzimuthAngle(currentThetaAngle)*tudat::physical_constants::JULIAN_YEAR  << ' '  << std::abs(invPolyShaping.computeCurrentThrustAccelerationMagnitude( currentThetaAngle )) * physical_constants::ASTRONOMICAL_UNIT / std::pow(tudat::physical_constants::JULIAN_YEAR, 2.0) << '\n';
//                            }
//                        }

//                    }

//                    peakAccelerationExposins = 0.0;
//                    infeasibleTOF = exposinsShaping.getInfeasibleTOF();
//                    if (infeasibleTOF)
//                    {
//                        deltaVexposins = 0.0;
//                        NrunsFailed++;
//                    }
//                    else
//                    {
//                        std::cout << "DeltaV - exposins " <<  exposinsShaping.computeDeltaV() << std::endl;
//                        std::cout << "DeltaVB - exposins " << exposinsShaping.computeDeltaVBoundaries() << std::endl;

//                        deltaVexposins = exposinsShaping.computeDeltaV() + exposinsShaping.computeDeltaVBoundaries();

//                        // Trajectory
//                        if (fileTrajectories.is_open())
//                        {

//                        double currentTravelledAngle = exposinsShaping.getTravelledAzimuthAngle();
//                        double stepSize = ( currentTravelledAngle ) / 100.0;
//                        Eigen::Vector6d currentStateVector;
//                        // Check that the trajectory is feasible, ie curved toward the central body.
//                            for ( int i = 0 ; i <= 100 ; i++ )
//                            {
//                                double currentThetaAngle = i * stepSize;

//                                currentStateVector = exposinsShaping.computeCurrentStateVector(currentThetaAngle);

//                                if ( std::abs(exposinsShaping.computeCurrentNormalizedThrustAccelerationMagnitude( currentThetaAngle )) >  peakAccelerationExposins )
//                                {
//                                    peakAccelerationExposins = std::abs(exposinsShaping.computeCurrentNormalizedThrustAccelerationMagnitude( currentThetaAngle ));
//                                }

//                                fileTrajectories << 0 << ' ' << Nrev << ' ' << startDate/tudat::physical_constants::JULIAN_DAY << ' ' << TOF << ' ' << currentStateVector[0] << ' ' << currentStateVector[1] << ' ' << currentStateVector[2] << ' ' << currentStateVector[3] << ' ' << currentStateVector[4] << ' ' << currentStateVector[5] << ' ' << exposinsShaping.computeCurrentTimeFromAzimuthAngle(currentThetaAngle)*tudat::physical_constants::JULIAN_YEAR  << ' ' << std::abs(exposinsShaping.computeCurrentNormalizedThrustAccelerationMagnitude( currentThetaAngle )) << '\n';
//                            }
//                        }

//                    }

//                    if (sphericalShapingInfeasibleTOF)
//                    {
//                        deltaVspherical = 0.0;
//                        NrunsFailed++;
//                    }

//                    peakAccelerationHodographic = 0.0;
//                    deltaVhodographic = hodographicShaping.computeDeltaV();
//                    std::cout << "DeltaV - hodographic " <<  hodographicShaping.computeDeltaV() << std::endl;
//                    // Trajectory
//                    if (fileTrajectories.is_open())
//                    {
//                    double stepSize = ( TOF * tudat::physical_constants::JULIAN_DAY ) / 100.0;
//                    Eigen::Vector6d currentStateVector;
//                    // Check that the trajectory is feasible, ie curved toward the central body.
//                        for ( int i = 0 ; i <= 100 ; i++ )
//                        {
//                            double currentTime = i * stepSize;

//                            currentStateVector = hodographicShaping.computeCurrentStateVector(currentTime);

//                            if ( hodographicShaping.computeCurrentThrustAccelerationVector( currentTime ).norm() >  peakAccelerationHodographic )
//                            {
//                                peakAccelerationHodographic = hodographicShaping.computeCurrentThrustAccelerationVector( currentTime ).norm();
//                            }


//                            fileTrajectories << 3 << ' ' << Nrev << ' ' << startDate/tudat::physical_constants::JULIAN_DAY << ' ' << TOF << ' ' << currentStateVector[0] << ' ' << currentStateVector[1] << ' ' << currentStateVector[2] << ' ' << currentStateVector[3] << ' ' << currentStateVector[4] << ' '<< currentStateVector[5] << ' ' << currentTime << ' ' << hodographicShaping.computeCurrentThrustAccelerationVector( currentTime ).norm() << '\n';
//                        }
//                    }



//                    file << Nrev << ' ' << startDate/tudat::physical_constants::JULIAN_DAY << ' ' << TOF << ' ' << deltaVexposins << ' ' << deltaVinvPoly << ' ' << deltaVspherical << ' ' << deltaVhodographic << ' ' << peakAccelerationExposins << ' ' << peakAccelerationInvPoly << ' ' << peakAccelerationSpherical << ' ' << peakAccelerationHodographic << '\n';

//                    Nruns++;
//                }
//            }

//        }
//    }
//    std::cout << "Nruns: "       << Nruns*4 << std::endl;
//    std::cout << "NrunsFailed: " << NrunsFailed << std::endl;



//    BOOST_CHECK_SMALL( std::fabs(  0.0 ), 100.0 );

//}


//    }



//! Spherical Shaping.
BOOST_AUTO_TEST_CASE( test_interplanetary_cases_spherical_transfer )
{
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


    //Tempel-1
    double centralBodyGravitationalParameter = bodyMap["Sun"]->getGravityFieldModel()->getGravitationalParameter( );



    // Define root finder settings (used to update the updated value of the free coefficient, so that it matches the required time of flight).
        std::shared_ptr< root_finders::RootFinderSettings > rootFinderSettings =
                std::make_shared< root_finders::RootFinderSettings >( root_finders::bisection_root_finder, 1.0e-6, 150 );



    // ranging values - EArth-Tempel-1
    double startDateLower = 0.0 * physical_constants::JULIAN_DAY; // MJD2000 -> 2020
    double startDateUpper = 6000.0 * physical_constants::JULIAN_DAY; // -
    double startDateStep  = 5.0 * physical_constants::JULIAN_DAY;

    double TOFLower = 400.0; // days / orbit is 2038 days
    double TOFUpper = 1500.0; // synodic period earth Mars however ecc
    double TOFStep  = 5.0;

    int NrevLower = 0;
    int NrevUpper = 2;

    // create variables
    double deltaVspherical;


    double peakAccelerationSpherical   = 0.0;

    int Nruns       = 1;
    int NrunsFailed = 1;


    // grid searches
    std::ofstream file("dataGridSearchShapeMethodsEarthTempelVerif_run2.txt");
    std::ofstream fileTrajectories("dataGridSearchTrajectoriesEarthTempelVerif_run1_run2.txt");

//     double startDate = 20.0*365.25 * physical_constants::JULIAN_DAY;
//     double TOF       = 500.0; // 500
//     int Nrev         = 1;

        for (int Nrev = NrevLower; Nrev <= NrevUpper; Nrev++)
        {
            for (double TOF = TOFLower; TOF <= TOFUpper; TOF = TOF + TOFStep)
            {
                std::cout << "-----------------------------------------------" << std::endl;
                std::cout << "LOOP - nRev: " << Nrev << " TOF: " << TOF << std::endl;
                std::cout << "-----------------------------------------------" << std::endl;
                for (double startDate = startDateLower; startDate <= startDateUpper; startDate = startDate + startDateStep)
                {
//                    std::cout << "-----------------------------------------------" << std::endl;
//                    std::cout << "LOOP - nRev: " << Nrev << " startDate: " << startDate/tudat::physical_constants::JULIAN_DAY << " TOF: " << TOF << std::endl;
//                    std::cout << "-----------------------------------------------" << std::endl;

             // Define integrator settings.
             std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings =
             std::make_shared< numerical_integrators::IntegratorSettings< double > > (
                 numerical_integrators::rungeKutta4, 0.0, TOF * tudat::physical_constants::JULIAN_DAY / 500.0 );

             // states
             // Earth
             Eigen::Vector6d initialState = pointerToDepartureBodyEphemeris->getCartesianState( startDate );
//              Mars
//             Eigen::Vector6d finalState = pointerToArrivalBodyEphemeris->getCartesianState(startDate + TOF * tudat::physical_constants::JULIAN_DAY );

            // Tempel-1
             Eigen::Vector6d keplerianState;
                     keplerianState(0) = 3.145692355253531*physical_constants::ASTRONOMICAL_UNIT;
                     keplerianState(1) = .5096307949367873;
                     keplerianState(2) = 10.47382730765828*mathematical_constants::PI/180.0;
                     keplerianState(3) = 179.203181927074	*mathematical_constants::PI/180.0;
                     keplerianState(4) = 68.74970463832437*mathematical_constants::PI/180.0;
             double meanMotion = .1766565332681526*mathematical_constants::PI/180.0/physical_constants::JULIAN_DAY;
             double initialEpoch = (2457533.5 - 2400000.5 - 51544.0)*tudat::physical_constants::JULIAN_DAY; // JD -> MJD -> MJD2000


             double initialMeanAnomaly = 347.7086591090512*mathematical_constants::PI/180.0;

             double tofEpoch         = TOF*tudat::physical_constants::JULIAN_DAY;
             double finalMeanAnomaly = initialMeanAnomaly + meanMotion*(startDate + tofEpoch - initialEpoch);
             double eccentricAnomaly = orbital_element_conversions::convertMeanAnomalyToEccentricAnomaly( keplerianState(1), finalMeanAnomaly );
             double trueAnomaly      = orbital_element_conversions::convertEccentricAnomalyToTrueAnomaly( eccentricAnomaly, keplerianState(1) );
             keplerianState(5)       = trueAnomaly;

             Eigen::Vector6d finalState = orbital_element_conversions::convertKeplerianToCartesianElements( keplerianState, centralBodyGravitationalParameter );


             try
             {
             // Compute shaped trajectory. - spherical shaping - constants influence not yet adjusted
             shape_based_methods::SphericalShaping sphericalShaping = shape_based_methods::SphericalShaping(
                     initialState, finalState, TOF * tudat::physical_constants::JULIAN_DAY,
                     Nrev, bodyMap, "Vehicle", "Sun", 0.0,
                     rootFinderSettings, -1.0e-1, 1.0e-1, integratorSettings );

             deltaVspherical = sphericalShaping.computeDeltaV();
//             std::cout << "DeltaV - spherical " <<  deltaVspherical << std::endl;

             }
             catch (const std::runtime_error& error)
             {
             deltaVspherical = 0.0;
             NrunsFailed++;
//             std::cout << "DeltaV - spherical FAIL"  << std::endl;

             }

             file << Nrev << ' ' << startDate/tudat::physical_constants::JULIAN_DAY << ' ' << TOF << ' ' << deltaVspherical << '\n';
            }
        }
        }

    std::cout << "Nruns: "       << Nruns*4 << std::endl;
    std::cout << "NrunsFailed: " << NrunsFailed << std::endl;



    BOOST_CHECK_SMALL( std::fabs(  0.0 ), 100.0 );

}



BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
