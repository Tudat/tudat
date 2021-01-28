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
BOOST_AUTO_TEST_SUITE( test_shapebased_coasting )


//! Earth Tempel-1.
BOOST_AUTO_TEST_CASE( test_shapebased_coasting_earth_tempel_transfer )
{

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            USING STATEMENTS              //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    using namespace tudat;
    using namespace tudat::simulation_setup;
    using namespace tudat::propagators;
    using namespace tudat::numerical_integrators;
    using namespace tudat::orbital_element_conversions;
    using namespace tudat::basic_mathematics;
    using namespace tudat::unit_conversions;




    spice_interface::loadStandardSpiceKernels( );

    int numberOfRevolutions = 1;
    double julianDate = (8174.5 - 15.0*365.25) * physical_constants::JULIAN_DAY;
    double  timeOfFlight = 580.0;

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
                std::make_shared< root_finders::RootFinderSettings >( root_finders::bisection_root_finder, 1.0e-6, 150 );


    // /////////////////////////// TEMPEL1 - START ////////////////////////////

    std::string kernelPath = input_output::getSpiceKernelPath( );
    spice_interface::clearSpiceKernels();
    spice_interface::loadSpiceKernelInTudat( kernelPath + "tempel1_2000_2020.bsp" );
    //    spice_interface::loadStandardSpiceKernels();
    frameOrigin      = "SUN";


    // /////////////////////////// TEMPEL1 - END ////////////////////////////


   
    double centralBodyGravitationalParameter = bodyMap["Sun"]->getGravityFieldModel()->getGravitationalParameter( );


   
    // //////////////////////////// //
    // ranging values - EarthTempel //
    // //////////////////////////// //

    // Start Date
    double startDateLower =  (20.0*365.25) * physical_constants::JULIAN_DAY; // MJD2000 -> 2020 //20
    double startDateUpper = startDateLower; // + (5.0*365.25 * physical_constants::JULIAN_DAY); // -synodic period earth Tempel however ecc
    double startDateStep  = 40.0 * physical_constants::JULIAN_DAY; // 40

//    // TOF total
//    double TOFtLower = 1000.0; // days / orbit is 2038 days
//    double TOFtUpper = 1000.0; // days //8900
//    double TOFtStep  = 100.0;

    // Number of Revolutions
    int NrevLower = 1;
    int NrevUpper = 1;

    // TOF Phase 1
    double TOF1Lower = 200.0;
    double TOF1Upper = 400.0;
    double TOF1Step  = 50.0;

    // TOF Phase 2 / Coasting
    double TOFCLower = 100.0;
    double TOFCUpper = 400.0;
    double TOFCStep  = 100.0;

    // TOF Phase 3
    double TOF3Lower = 100.0;
    double TOF3Upper = 400.0;
    double TOF3Step  = 100.0;




    // create variables
    double deltaVinvPoly1;
    double deltaVinvPoly3;
    double deltaVspherical1;
    double deltaVspherical3;

    double peakAccelerationInvPoly1      = 0.0;
    double peakAccelerationInvPoly3      = 0.0;
    double peakAccelerationSpherical1   = 0.0;
    double peakAccelerationSpherical3   = 0.0;


    double currentTravelledAngle;
    double stepSize;
    Eigen::Vector6d currentStateVector;



    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////     CREATE ENVIRONMENT AND VEHICLE       //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Load Spice kernels.
    spice_interface::loadStandardSpiceKernels( );


//                    // Create body objects.
//                    std::vector< std::string > bodiesToCreate;
//                    bodiesToCreate.push_back( "Earth" );
//                    std::map< std::string, std::shared_ptr< BodySettings > > bodySettings =
//                            getDefaultBodySettings( bodiesToCreate );
//                    bodySettings[ "Earth" ]->ephemerisSettings = std::make_shared< ConstantEphemerisSettings >(
//                                Eigen::Vector6d::Zero( ) );

//                    // Create Earth object
//                    NamedBodyMap bodyMap = createBodies( bodySettings );

//                    // Create spacecraft object.
//                    bodyMap[ "Asterix" ] = std::make_shared< simulation_setup::Body >( );

//                    // Finalize body creation.
//                    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            CREATE ACCELERATIONS          //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Define propagator settings variables.
    SelectedAccelerationMap accelerationMap;
    std::vector< std::string > bodiesToPropagate;
    std::vector< std::string > centralBodies;

    bodiesToPropagate.push_back( "Vehicle" );
    centralBodies.push_back( "Sun" );

    // Define propagation settings.
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfVehicle;
    accelerationsOfVehicle[ "Sun" ].push_back( std::make_shared< AccelerationSettings >(
                                                     basic_astrodynamics::central_gravity ) );
    accelerationMap[ "Vehicle" ] = accelerationsOfVehicle;

    // Create acceleration models and propagation settings.
    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodyMap, accelerationMap, bodiesToPropagate, centralBodies );

    // Create numerical integrator settings.
    double simulationStartEpoch = 0.0;
    const double fixedStepSize = physical_constants::JULIAN_DAY;
    std::shared_ptr< IntegratorSettings< > > integratorSettings =
            std::make_shared< IntegratorSettings< > >
            ( rungeKutta4, simulationStartEpoch, fixedStepSize );


    /////////////////////////////////////////////////////////////////////////////////

    int Nruns       = 0;
    int Nfails      = 0;

    bool gridSearchFlag                   = false;
    bool optimumTrajectoriesFlag          = false;
    bool multiOptimumTrajectoriesFlag     = true;
    bool multiIPSSOptimumTrajectoriesFlag = false;

    Eigen::MatrixXd optimumInputs;



        // multi - objective
        std::ofstream file("dataGridSearchShapeMethodsCoastingEarthTempel1_runMO30.txt");
        std::ofstream fileTrajectories("dataGridSearchCoastingTrajectoriesEarthTempel_runMO30.txt");

        optimumInputs =
                 tudat::input_output::readMatrixFromFile( "C:/tudatBundle/tudatExampleApplications/libraryExamples/PaGMOEx/SimulationOutput/population_199_Shaping_x_EarthTempelSBIPCRev0_runMO42"
                                                          ".dat"); // , " \t", "#"

//        std::ofstream file("dataGridSearchShapeMethodsCoastingEarthTempel1_runO38.txt");
//        std::ofstream fileTrajectories("dataGridSearchCoastingTrajectoriesEarthTempel_runO38.txt");

//        optimumInputs =
//                 tudat::input_output::readMatrixFromFile( "C:/tudatBundle/tudatExampleApplications/libraryExamples/PaGMOEx/SimulationOutput/population_999_Shaping_x_EarthTempelSBCRev0_runO38"
//                                                          ".dat"); // , " \t", "#"

    // grid searches



//  fileTrajectories << Method << ' ' << Nrev << ' ' << startDate << ' ' << TOFt << ' ' << TOF1 << ' ' << TOFC << ' ' << Phase << ' ' << initialState[0-5] << ' ' << Time << ' ' << Acc << '\n';


    bool SSTrajectories = false;
    bool IPTrajectories = true;

    bool mSSTrajectories = false;
    bool mIPTrajectories = false;


    std::cout << "optimumInputs: " << optimumInputs(1,1) << std::endl;





    if (gridSearchFlag)
    {
    // GRID SEARCH
    {
    if (file.is_open())
    {
        for (int Nrev = NrevLower; Nrev <= NrevUpper; Nrev++)
        {               
            for (double startDate = startDateLower; startDate <= startDateUpper; startDate = startDate + startDateStep)
            {

                for (double TOF1 = TOF1Lower; TOF1 <= TOF1Upper; TOF1 = TOF1 + TOF1Step)
                {
                    std::cout << "-----------------------------------------------" << std::endl;
                    std::cout << "LOOP - nRev: " << Nrev << " startDate: " << startDate/tudat::physical_constants::JULIAN_DAY << " TOF1 " << TOF1<< std::endl;
                    std::cout << "-----------------------------------------------" << std::endl;
                    for (double TOF3 = TOF3Lower; TOF3 <= TOF3Upper; TOF3 = TOF3 + TOF3Step)
                    {
                        for (double TOFC = TOFCLower; TOFC <= TOFCUpper; TOFC = TOFC + TOFCStep)
                            {
                                // TOFs
                                double TOFt = TOF1 + TOFC + TOF3;



                                // Delta V initialisation


                                // Tempel-1 State
                                Eigen::Vector6d keplerianState;
                                        keplerianState(0) = 3.145692355253531*physical_constants::ASTRONOMICAL_UNIT;
                                        keplerianState(1) = .5096307949367873;
                                        keplerianState(2) = 10.47382730765828*mathematical_constants::PI/180.0;
                                        keplerianState(3) = 179.203181927074	*mathematical_constants::PI/180.0;
                                        keplerianState(4) = 68.74970463832437*mathematical_constants::PI/180.0;
                                double meanMotion = .1766565332681526*mathematical_constants::PI/180.0/physical_constants::JULIAN_DAY;
                                double initialEpoch = (2457533.5 - 2400000.5 - 51544.0)*tudat::physical_constants::JULIAN_DAY; // JD -> MJD -> MJD2000


                                double initialMeanAnomaly = 347.7086591090512*mathematical_constants::PI/180.0;

                                double tofEpoch         = TOFt*tudat::physical_constants::JULIAN_DAY;
                                double finalMeanAnomaly = initialMeanAnomaly + meanMotion*(startDate + tofEpoch - initialEpoch);
                                double eccentricAnomaly = orbital_element_conversions::convertMeanAnomalyToEccentricAnomaly( keplerianState(1), finalMeanAnomaly );
                                double trueAnomaly      = orbital_element_conversions::convertEccentricAnomalyToTrueAnomaly( eccentricAnomaly, keplerianState(1) );
                                keplerianState(5)       = trueAnomaly;

                                Eigen::Vector6d cartesianState = orbital_element_conversions::convertKeplerianToCartesianElements( keplerianState, centralBodyGravitationalParameter );




                                // states
                                // Earth
                                Eigen::Vector6d initialState = pointerToDepartureBodyEphemeris->getCartesianState( startDate );

                                // Tempel 1
                                Eigen::Vector6d finalState = cartesianState;


//                                std::cout << "MJD to Calendar (Start) " << basic_astrodynamics::convertJulianDayToCalendarDate(basic_astrodynamics::convertModifiedJulianDayToJulianDay(startDate/tudat::physical_constants::JULIAN_DAY + 51544.5)) << std::endl;
//                                std::cout << "MJD to Calendar (End) " << basic_astrodynamics::convertJulianDayToCalendarDate(basic_astrodynamics::convertModifiedJulianDayToJulianDay(startDate/tudat::physical_constants::JULIAN_DAY + 51544.5 + TOFt)) << std::endl;
//                                std::cout << "MJD to Calendar (End) " << basic_astrodynamics::convertJulianDayToCalendarDate(spice_interface::convertEphemerisTimeToJulianDate( (startDate) + TOFt*tudat::physical_constants::JULIAN_DAY)) << std::endl;



                                /// ################################################################################################################ //
                                /// ////////////////////        0. COASTING SET-UP          ///////////////////////////////////////////////////////////
                                /// ################################################################################################################ //

                                // coasting set-up
                                // spherical state (r, theta, phi, Vr, Vtheta, Vphi).


                                Eigen::Vector6d initialStateSphericalCoordinates = coordinate_conversions::convertCartesianToSphericalState( initialState );
                                Eigen::Vector6d finalStateSphericalCoordinates = coordinate_conversions::convertCartesianToSphericalState( finalState );


                                // Compute initial and final value of the azimuth angle.
                                double initialAzimuthAngle = initialStateSphericalCoordinates( 1 );
                                double finalAzimuthAngle;

                                if ( ( finalStateSphericalCoordinates( 1 ) - initialStateSphericalCoordinates( 1 ) ) < 0.0 )
                                {
                                    finalAzimuthAngle = finalStateSphericalCoordinates( 1 ) + 2.0 * mathematical_constants::PI * ( Nrev + 1.0 );
                                }
                                else
                                {
                                    finalAzimuthAngle = finalStateSphericalCoordinates( 1 ) + 2.0 * mathematical_constants::PI * ( Nrev );
                                }

                                // Position C
                                // radius
                                double radiusLower = 0.8; // initial Orbit
                                double radiusUpper = 2.0; // final Orbit
                                double radiusStep  = 0.4;

                                // azimuth angle
                                double azimuthStep  = (finalAzimuthAngle - initialAzimuthAngle)/4.0; //
                                double azimuthLower = initialAzimuthAngle + azimuthStep; // initialAzimuthAngle
                                double azimuthUpper = finalAzimuthAngle - azimuthStep; // finalAzimuthAngle

                                //elevation
                                double elevationLower = 1.0; // ???
                                double elevationUpper = 1.0; // ???
                                double elevationStep  = 1.0; //

                                // Velocity C
                                // radial Velociy
                                double VrLower = 1.0;
                                double VrUpper = 1.0;
                                double VrStep  = 1.0;

                                // angular velocity
                                double VthetaLower = 0.9;
                                double VthetaUpper = 1.1;
                                double VthetaStep  = 0.1;

                                // ??? velocity
                                double VphiLower = 1.0;
                                double VphiUpper = 1.0;
                                double VphiStep  = 0.1;


                                for ( double radius = radiusLower ; radius <= radiusUpper ; radius = radius + radiusStep )
                                {
                                    for ( double azimuth = azimuthLower ; azimuth <= azimuthUpper ; azimuth = azimuth + azimuthStep )
                                    {
                                        for ( double elevation = elevationLower ; elevation <= elevationUpper ; elevation = elevation + elevationStep )
                                        {
                                            for ( double Vr = VrLower ; Vr <= VrUpper ; Vr = Vr + VrStep )
                                            {
                                                for ( double Vtheta = VthetaLower ; Vtheta <= VthetaUpper ; Vtheta = Vtheta + VthetaStep )
                                                {
                                                    for ( double Vphi = VphiLower ; Vphi <= VphiUpper ; Vphi = Vphi + VphiStep )
                                                    {
                                                        Nruns++;
                                                        Eigen::Vector6d sphericalState;
                                                        sphericalState(0) = radius*physical_constants::ASTRONOMICAL_UNIT; // r
                                                        sphericalState(1) = azimuth; // theta
                                                        sphericalState(2) = 0.0; // phi
                                                        sphericalState(3) = Vr*initialStateSphericalCoordinates(3); // Vr
                                                        sphericalState(4) = Vtheta*initialStateSphericalCoordinates(4); // Vtheta
                                                        sphericalState(5) = 0.0; // Vphi


                                                        Eigen::Vector6d coastingState1 = coordinate_conversions::convertSphericalToCartesianState( sphericalState );
                                    //                    Eigen::Vector6d coastingState1 = initialState;




                                                        fileTrajectories << Nruns << ' ' << 2 << ' ' << Nrev << ' ' << startDate/tudat::physical_constants::JULIAN_DAY << ' ' << TOFt << ' ' << TOF1 << ' ' << TOFC << ' ' << radius << ' ' << azimuth << ' ' << elevation << ' ' << Vr << ' ' << Vtheta << ' ' << Vphi << ' ' << 1 << ' ' << initialState[0] << ' ' << initialState[1] << ' ' << initialState[2] << ' ' << initialState[3] << ' ' << initialState[4] << ' '<< initialState[5] << ' ' << 0 << ' ' << 0 << '\n';


                                                        /// ################################################################################################################ //
                                                        /// ////////////////////        1. PHASE 1 - POWERED FLIGHT         ///////////////////////////////////////////////////
                                                        /// ################################################################################################################ //
                                                        bool sphericalShapingInfeasibleTOF1;
                                                        bool sphericalShapingInfeasibleTOF3;

                                                        bool invPolyShapingInfeasableTOF1;
                                                        bool invPolyShapingInfeasableTOF3;

                                                        peakAccelerationSpherical1   = 0.0;
                                                        try
                                                        {

                                                        // Define integrator settings.
                                                            std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettingsSB1 =
                                                                    std::make_shared< numerical_integrators::IntegratorSettings< double > > (
                                                                        numerical_integrators::rungeKutta4, 0.0, TOF1 * tudat::physical_constants::JULIAN_DAY / 500.0 );

//                                                        // Compute shaped trajectory. - spherical shaping - constants influence not yet adjusted
//                                                        shape_based_methods::SphericalShaping sphericalShaping = shape_based_methods::SphericalShaping(
//                                                                    initialState, coastingState1, TOF1 * tudat::physical_constants::JULIAN_DAY,
//                                                                    Nrev, bodyMap, "Vehicle", "Sun", 0.0,
//                                                                    rootFinderSettings, -1.0e-1, 1.0e-1, integratorSettingsSB1 );

                                                         // Compute shaped trajectory. - invpolyshaping
                                                        shape_based_methods::InvPolyShaping invPolyShaping1 = shape_based_methods::InvPolyShaping(
                                                                    initialState, coastingState1, TOF1 * tudat::physical_constants::JULIAN_DAY,
                                                                    Nrev, bodyMap, "Vehicle", "Sun", 0.928,
                                                                    rootFinderSettings, -1.0e-3, 1.0, integratorSettingsSB1 );




                                                        sphericalShapingInfeasibleTOF1 = false;
                                                        invPolyShapingInfeasableTOF1   = false;

//                                                        deltaVspherical1 = sphericalShaping.computeDeltaV();
                                                        deltaVinvPoly1   = invPolyShaping1.computeDeltaV();

                                                        // std::cout << "DeltaV - spherical " <<  deltaVspherical1 << std::endl;


                                                        if (fileTrajectories.is_open())
                                                        {
//                                                            // spherical
//                                                        currentTravelledAngle = sphericalShaping.getFinalAzimuthAngle() - sphericalShaping.getInitialAzimuthAngle();
//                                                        stepSize = ( currentTravelledAngle ) / 100.0;

//                                                        Eigen::Vector6d currentStateVector;
//                                                        // Check that the trajectory is feasible, ie curved toward the central body.
//                                                        for ( int i = 0 ; i <= 100 ; i++ )
//                                                        {
//                                                                double currentThetaAngle = sphericalShaping.getInitialAzimuthAngle() + i * stepSize;

//                                                                currentStateVector = sphericalShaping.computeCurrentStateVector(currentThetaAngle);

//                                                                if ( sphericalShaping.computeCurrentThrustAccelerationVector( currentThetaAngle ).norm() >  peakAccelerationSpherical1 )
//                                                                {
//                                                                    peakAccelerationSpherical1 = sphericalShaping.computeCurrentThrustAccelerationVector( currentThetaAngle ).norm();

//                                                                }

//                                                                fileTrajectories << Nruns << ' ' << 2 << ' ' << Nrev << ' ' << startDate/tudat::physical_constants::JULIAN_DAY << ' ' << TOFt << ' ' << TOF1 << ' ' << TOFC << ' ' << radius << ' ' << azimuth << ' ' << elevation << ' ' << Vr << ' ' << Vtheta << ' ' << Vphi << ' ' << 1 << ' ' << currentStateVector[0] << ' ' << currentStateVector[1] << ' ' << currentStateVector[2] << ' ' << currentStateVector[3] << ' ' << currentStateVector[4] << ' '<< currentStateVector[5] << ' ' << sphericalShaping.convertIndependentVariableToTime(currentThetaAngle)*tudat::physical_constants::JULIAN_YEAR << ' ' << sphericalShaping.computeCurrentThrustAccelerationVector( currentThetaAngle ).norm() << '\n';
//                                                        }

                                                        currentTravelledAngle = invPolyShaping1.getTravelledAzimuthAngle();
                                                        stepSize = ( currentTravelledAngle ) / 100.0;
                                                        // Check that the trajectory is feasible, ie curved toward the central body.
                                                            for ( int i = 0 ; i <= 100 ; i++ )
                                                            {
                                                                double currentThetaAngle = i * stepSize;

                                                                currentStateVector = invPolyShaping1.computeCurrentStateVector(currentThetaAngle);

                                                                if ( (std::abs(invPolyShaping1.computeCurrentThrustAccelerationMagnitude( currentThetaAngle ))  * physical_constants::ASTRONOMICAL_UNIT / std::pow(tudat::physical_constants::JULIAN_YEAR, 2.0) ) >  peakAccelerationInvPoly1 )
                                                                {
                                                                    // thrust is non-dimensional so dimensionalisation is needed
                                                                    peakAccelerationInvPoly1 = std::abs(invPolyShaping1.computeCurrentThrustAccelerationMagnitude( currentThetaAngle )) * physical_constants::ASTRONOMICAL_UNIT / std::pow(tudat::physical_constants::JULIAN_YEAR, 2.0);
                                                                }

                                                                fileTrajectories << Nruns << ' ' << 1 << ' ' << Nrev << ' ' << startDate/tudat::physical_constants::JULIAN_DAY << ' ' << TOFt << ' ' << TOF1 << ' ' << TOFC << ' ' << radius << ' ' << azimuth << ' ' << elevation << ' ' << Vr << ' ' << Vtheta << ' ' << Vphi << ' ' << 1 << ' ' << currentStateVector[0] << ' ' << currentStateVector[1] << ' ' << currentStateVector[2] << ' ' << currentStateVector[3] << ' ' << currentStateVector[4] << ' '<< currentStateVector[5] << ' ' << invPolyShaping1.convertIndependentVariableToTime(currentThetaAngle)*tudat::physical_constants::JULIAN_YEAR << ' ' << std::abs(invPolyShaping1.computeCurrentThrustAccelerationMagnitude( currentThetaAngle )) * physical_constants::ASTRONOMICAL_UNIT / std::pow(tudat::physical_constants::JULIAN_YEAR, 2.0) << '\n';
                                                            } // for loop stateVector

                                                        } // if file sopen
                                                        } // try
                                                        catch (const std::runtime_error& error)
                                                        {
                                                            sphericalShapingInfeasibleTOF1 = true;
                                                            deltaVspherical1 = 0.0;
                                                            // std::cout << "DeltaV - spherical " <<  deltaVspherical1 << std::endl;
                                                        } // catch

                                                        /// ################################################################################################################ //
                                                        /// ////////////////////        2. PHASE 2 - COASTING               ///////////////////////////////////////////////////
                                                        /// ################################################################################################################ //


                                                        if (sphericalShapingInfeasibleTOF1==false) // && deltaVspherical1<100e3
                                                        {

                                                            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                                                            ///////////////////////             CREATE PROPAGATION SETTINGS            ////////////////////////////////////////////
                                                            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                                                            // Set initial conditions for the Asterix satellite that will be propagated in this simulation.
                                                            // The initial conditions are given in Keplerian elements and later on converted to Cartesian
                                                            // elements.

                                                            // Set simulation end epoch.
                                                            const double simulationEndEpoch = TOFC*tudat::physical_constants::JULIAN_DAY;

                                                            // Create propagator settings.
                                                            std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
                                                                    std::make_shared< TranslationalStatePropagatorSettings< double > >
                                                                    ( centralBodies, accelerationModelMap, bodiesToPropagate, coastingState1, simulationEndEpoch );



                                                            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                                                            ///////////////////////             PROPAGATE ORBIT            ////////////////////////////////////////////////////////
                                                            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                                                            // Create simulation object and propagate dynamics.
                                                            SingleArcDynamicsSimulator< > dynamicsSimulator( bodyMap, integratorSettings, propagatorSettings );
                                                            std::map< double, Eigen::VectorXd > integrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );



                                                            Eigen::VectorXd coastingState2 = ( --integrationResult.end( ) )->second;

                            //                                std::cout << "initialState " << initialState << std::endl;
                            //                                std::cout << "coastingState1 " << coastingState1 << std::endl;
                            //                                std::cout << "coastingState2 " << coastingState2 << std::endl;
                            //                                std::cout << "finalState " << finalState << std::endl;

                                                            fileTrajectories << Nruns << ' ' << 2 << ' ' << Nrev << ' ' << startDate/tudat::physical_constants::JULIAN_DAY << ' ' << TOFt << ' ' << TOF1 << ' ' << TOFC << ' ' << radius << ' ' << azimuth << ' ' << elevation << ' ' << Vr << ' ' << Vtheta << ' ' << Vphi << ' ' << 2 << ' ' << coastingState1[0] << ' ' << coastingState1[1] << ' ' << coastingState1[2] << ' ' << coastingState1[3] << ' ' << coastingState1[4] << ' '<< coastingState1[5] << ' ' << 0 << ' ' << 0 << '\n';
                                                            fileTrajectories << Nruns << ' ' << 2 << ' ' << Nrev << ' ' << startDate/tudat::physical_constants::JULIAN_DAY << ' ' << TOFt << ' ' << TOF1 << ' ' << TOFC << ' ' << radius << ' ' << azimuth << ' ' << elevation << ' ' << Vr << ' ' << Vtheta << ' ' << Vphi << ' ' << 2 << ' ' << coastingState2[0] << ' ' << coastingState2[1] << ' ' << coastingState2[2] << ' ' << coastingState2[3] << ' ' << coastingState2[4] << ' '<< coastingState2[5] << ' ' << 0 << ' ' << 0 << '\n';


                            //                                // Write satellite propagation history to file.
                            //                                input_output::writeDataMapToTextFile( integrationResult,
                            //                                                                      "shapeBasedCoastingPropagationHistory_run2.dat",
                            //                                                                      "C:/tudatBundle/tudat/bin/unit_tests/",
                            //                                                                      "",
                            //                                                                      std::numeric_limits< double >::digits10,
                            //                                                                      std::numeric_limits< double >::digits10,
                            //                                                                      "," );




                                                            /// ################################################################################################################ //
                                                            /// ////////////////////        3. PHASE 3 - POWERED FLIGHT               /////////////////////////////////////////////
                                                            /// ################################################################################################################ //

                                                            peakAccelerationSpherical3   = 0.0;
                                                            try
                                                            {

                                                            // Define integrator settings.
                                                                    std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettingsSB2 =
                                                                            std::make_shared< numerical_integrators::IntegratorSettings< double > > (
                                                                                numerical_integrators::rungeKutta4, 0.0, TOF3 * tudat::physical_constants::JULIAN_DAY / 500.0 );

//                                                            // Compute shaped trajectory. - spherical shaping - constants influence not yet adjusted
//                                                            shape_based_methods::SphericalShaping sphericalShaping = shape_based_methods::SphericalShaping(
//                                                                        coastingState2, finalState, TOF3 * tudat::physical_constants::JULIAN_DAY,
//                                                                        Nrev, bodyMap, "Vehicle", "Sun", 0.0,
//                                                                        rootFinderSettings, -1.0e-1, 1.0e-1, integratorSettingsSB2 );

                                                            // Compute shaped trajectory. - invpolyshaping
                                                            shape_based_methods::InvPolyShaping invPolyShaping3 = shape_based_methods::InvPolyShaping(
                                                                        coastingState2, finalState, TOF3 * tudat::physical_constants::JULIAN_DAY,
                                                                        Nrev, bodyMap, "Vehicle", "Sun", 0.928,
                                                                        rootFinderSettings, -1.0e-3, 1.0, integratorSettingsSB2 );




                                                            sphericalShapingInfeasibleTOF3 = false;
                                                            invPolyShapingInfeasableTOF3   = false;


//                                                            deltaVspherical3 = sphericalShaping.computeDeltaV();
                                                            deltaVinvPoly3   = invPolyShaping3.computeDeltaV();

                                                            // std::cout << "DeltaV - spherical " <<  deltaVspherical3 << std::endl;


                                                            if (fileTrajectories.is_open())
                                                            {
//                                                            // spherical
//                                                            currentTravelledAngle = sphericalShaping.getFinalAzimuthAngle() - sphericalShaping.getInitialAzimuthAngle();
//                                                            stepSize = ( currentTravelledAngle ) / 100.0;

//                                                            // Check that the trajectory is feasible, ie curved toward the central body.
//                                                            for ( int i = 0 ; i <= 100 ; i++ )
//                                                            {
//                                                                    double currentThetaAngle = sphericalShaping.getInitialAzimuthAngle() + i * stepSize;

//                                                                    currentStateVector = sphericalShaping.computeCurrentStateVector(currentThetaAngle);

//                                                                    if ( sphericalShaping.computeCurrentThrustAccelerationVector( currentThetaAngle ).norm() >  peakAccelerationSpherical3 )
//                                                                    {
//                                                                        peakAccelerationSpherical3 = sphericalShaping.computeCurrentThrustAccelerationVector( currentThetaAngle ).norm();

//                                                                    }

//                                                                    fileTrajectories << Nruns << ' ' << 2 << ' ' << Nrev << ' ' << startDate/tudat::physical_constants::JULIAN_DAY << ' ' << TOFt << ' ' << TOF1 << ' ' << TOFC << ' ' << radius << ' ' << azimuth << ' ' << elevation << ' ' << Vr << ' ' << Vtheta << ' ' << Vphi << ' ' << 3 << ' ' << currentStateVector[0] << ' ' << currentStateVector[1] << ' ' << currentStateVector[2] << ' ' << currentStateVector[3] << ' ' << currentStateVector[4] << ' '<< currentStateVector[5] << ' ' << sphericalShaping.convertIndependentVariableToTime(currentThetaAngle)*tudat::physical_constants::JULIAN_YEAR << ' ' << sphericalShaping.computeCurrentThrustAccelerationVector( currentThetaAngle ).norm() << '\n';
//                                                            }
                                                            currentTravelledAngle = invPolyShaping3.getTravelledAzimuthAngle();
                                                            stepSize = ( currentTravelledAngle ) / 100.0;
                                                            // Check that the trajectory is feasible, ie curved toward the central body.
                                                                for ( int i = 0 ; i <= 100 ; i++ )
                                                                {
                                                                    double currentThetaAngle = i * stepSize;

                                                                    currentStateVector = invPolyShaping3.computeCurrentStateVector(currentThetaAngle);

                                                                    if ( (std::abs(invPolyShaping3.computeCurrentThrustAccelerationMagnitude( currentThetaAngle ))  * physical_constants::ASTRONOMICAL_UNIT / std::pow(tudat::physical_constants::JULIAN_YEAR, 2.0) ) >  peakAccelerationInvPoly1 )
                                                                    {
                                                                        // thrust is non-dimensional so dimensionalisation is needed
                                                                        peakAccelerationInvPoly3 = std::abs(invPolyShaping3.computeCurrentThrustAccelerationMagnitude( currentThetaAngle )) * physical_constants::ASTRONOMICAL_UNIT / std::pow(tudat::physical_constants::JULIAN_YEAR, 2.0);
                                                                    }

                                                                    fileTrajectories << Nruns << ' ' << 1 << ' ' << Nrev << ' ' << startDate/tudat::physical_constants::JULIAN_DAY << ' ' << TOFt << ' ' << TOF1 << ' ' << TOFC << ' ' << radius << ' ' << azimuth << ' ' << elevation << ' ' << Vr << ' ' << Vtheta << ' ' << Vphi << ' ' << 1 << ' ' << currentStateVector[0] << ' ' << currentStateVector[1] << ' ' << currentStateVector[2] << ' ' << currentStateVector[3] << ' ' << currentStateVector[4] << ' '<< currentStateVector[5] << ' ' << invPolyShaping3.convertIndependentVariableToTime(currentThetaAngle)*tudat::physical_constants::JULIAN_YEAR << ' ' << std::abs(invPolyShaping3.computeCurrentThrustAccelerationMagnitude( currentThetaAngle )) * physical_constants::ASTRONOMICAL_UNIT / std::pow(tudat::physical_constants::JULIAN_YEAR, 2.0) << '\n';
                                                                } // for loop stateVector

                                                            }
                                                            }
                                                            catch (const std::runtime_error& error)
                                                            {
                                                                sphericalShapingInfeasibleTOF3 = true;
                                                                deltaVspherical3 = 0.0;
                                                                // std::cout << "DeltaV - spherical " <<  deltaVspherical3 << std::endl;

                                                            }


//                                                            file << Nruns << ' ' << 2 << ' ' << Nrev << ' ' << startDate/tudat::physical_constants::JULIAN_DAY << ' ' << TOFt << ' ' << TOF1 << ' ' << TOFC << ' ' << radius << ' ' << azimuth << ' ' << elevation << ' ' << Vr << ' ' << Vtheta << ' ' << Vphi << ' ' << deltaVspherical1 << ' ' << deltaVspherical3 << ' ' << peakAccelerationSpherical1 << ' ' << peakAccelerationSpherical3 << ' ' << '\n';
                                                            file << Nruns << ' ' << 1 << ' ' << Nrev << ' ' << startDate/tudat::physical_constants::JULIAN_DAY << ' ' << TOFt << ' ' << TOF1 << ' ' << TOFC << ' ' << radius << ' ' << azimuth << ' ' << elevation << ' ' << Vr << ' ' << Vtheta << ' ' << Vphi << ' ' << deltaVinvPoly1 << ' ' << deltaVinvPoly3 << ' ' << peakAccelerationInvPoly1 << ' ' << peakAccelerationInvPoly3 << ' ' << '\n';



                                                        }// if deltaVspherical1
                                                        else
                                                        {
                                                            Nfails++;
                                                        }




                                                    }// Vphi LOOP
                                                }//Vtheta LOOP
                                            }// Vr LOOP
                                        }// elevation LOOP
                                    }// azimuth LOOP
                                }// radius LOOP



                            } // TOFC LOOP
                    }// TOF3 LOOP
                }// TOF1 LOOP
            }// Start Date LOOP
        }// Nrev LOOP
    }// file.is_open()



    std::cout << "Coasting program done: Nruns = " << Nruns << " Nfails = " << Nfails << std::endl;



    BOOST_CHECK_SMALL( std::fabs(  0.0 ), 100.0 );



    }

    }


    if (optimumTrajectoriesFlag)
    {
        for ( int nOptimum = 0 ; nOptimum < 100 ; nOptimum++ )
        {

            std::cout << "nPop: " << nOptimum << std::endl;


            Nruns++;
            // TOFs
            double Nrev = 0;

            double startDate = optimumInputs(nOptimum, 0 );

            double TOF1      = optimumInputs(nOptimum, 1 )/tudat::physical_constants::JULIAN_DAY;
            double TOFC      = optimumInputs(nOptimum, 2 )/tudat::physical_constants::JULIAN_DAY;
            double TOF3      = optimumInputs(nOptimum, 3 )/tudat::physical_constants::JULIAN_DAY;

            double radius    = optimumInputs(nOptimum, 4 );
            double azimuth   = optimumInputs(nOptimum, 5 );
            double elevation = optimumInputs(nOptimum, 6 );
            double Vr        = optimumInputs(nOptimum, 7 );
            double Vtheta    = optimumInputs(nOptimum, 8 );
            double Vphi      = optimumInputs(nOptimum, 9 );

            double TOFt = TOF1 + TOFC + TOF3;

            // Tempel-1 State
            Eigen::Vector6d keplerianState;
                    keplerianState(0) = 3.145692355253531*physical_constants::ASTRONOMICAL_UNIT;
                    keplerianState(1) = .5096307949367873;
                    keplerianState(2) = 10.47382730765828*mathematical_constants::PI/180.0;
                    keplerianState(3) = 179.203181927074	*mathematical_constants::PI/180.0;
                    keplerianState(4) = 68.74970463832437*mathematical_constants::PI/180.0;
            double meanMotion = .1766565332681526*mathematical_constants::PI/180.0/physical_constants::JULIAN_DAY;
            double initialEpoch = (2457533.5 - 2400000.5 - 51544.0)*tudat::physical_constants::JULIAN_DAY; // JD -> MJD -> MJD2000


            double initialMeanAnomaly = 347.7086591090512*mathematical_constants::PI/180.0;

            double tofEpoch         = TOFt*tudat::physical_constants::JULIAN_DAY;
            double finalMeanAnomaly = initialMeanAnomaly + meanMotion*(startDate + tofEpoch - initialEpoch);
            double eccentricAnomaly = orbital_element_conversions::convertMeanAnomalyToEccentricAnomaly( keplerianState(1), finalMeanAnomaly );
            double trueAnomaly      = orbital_element_conversions::convertEccentricAnomalyToTrueAnomaly( eccentricAnomaly, keplerianState(1) );
            keplerianState(5)       = trueAnomaly;

            Eigen::Vector6d cartesianState = orbital_element_conversions::convertKeplerianToCartesianElements( keplerianState, centralBodyGravitationalParameter );

            // States
            // Earth
            Eigen::Vector6d initialState = pointerToDepartureBodyEphemeris->getCartesianState( startDate );
            // Tempel 1
            Eigen::Vector6d finalState = cartesianState;



            /// ################################################################################################################ //
            /// ////////////////////        0. COASTING SET-UP          ///////////////////////////////////////////////////////////
            /// ################################################################################################################ //


            Eigen::Vector6d sphericalState;
            sphericalState(0) = radius*physical_constants::ASTRONOMICAL_UNIT; // r
            sphericalState(1) = azimuth; // theta
            sphericalState(2) = elevation; // phi
            sphericalState(3) = Vr; // Vr
            sphericalState(4) = Vtheta; // Vtheta
            sphericalState(5) = Vphi; // Vphi


            Eigen::Vector6d coastingState1 = coordinate_conversions::convertSphericalToCartesianState( sphericalState );


            fileTrajectories << Nruns << ' ' << 2 << ' ' << Nrev << ' ' << startDate/tudat::physical_constants::JULIAN_DAY << ' ' << TOFt << ' ' << TOF1 << ' ' << TOFC << ' ' << radius << ' ' << azimuth << ' ' << elevation << ' ' << Vr << ' ' << Vtheta << ' ' << Vphi << ' ' << 1 << ' ' << initialState[0] << ' ' << initialState[1] << ' ' << initialState[2] << ' ' << initialState[3] << ' ' << initialState[4] << ' '<< initialState[5] << ' ' << 0 << ' ' << 0 << '\n';


            /// ################################################################################################################ //
            /// ////////////////////        1. PHASE 1 - POWERED FLIGHT         ///////////////////////////////////////////////////
            /// ################################################################################################################ //
            bool sphericalShapingInfeasibleTOF1;
            bool sphericalShapingInfeasibleTOF3;

            bool invPolyShapingInfeasableTOF1;
            bool invPolyShapingInfeasableTOF3;

            try
            {

            // Define integrator settings.
                std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettingsSB1 =
                        std::make_shared< numerical_integrators::IntegratorSettings< double > > (
                            numerical_integrators::rungeKutta4, 0.0, TOF1 * tudat::physical_constants::JULIAN_DAY / 500.0 );


            if (SSTrajectories)
            {
                // Compute shaped trajectory. - spherical shaping - constants influence not yet adjusted
                shape_based_methods::SphericalShaping sphericalShaping1 = shape_based_methods::SphericalShaping(
                            initialState, coastingState1, TOF1 * tudat::physical_constants::JULIAN_DAY,
                            Nrev, bodyMap, "Vehicle", "Sun", 0.0,
                            rootFinderSettings, -1.0e-1, 1.0e-1, integratorSettingsSB1 );

                deltaVspherical1 = sphericalShaping1.computeDeltaV();

                currentTravelledAngle = sphericalShaping1.getFinalAzimuthAngle() - sphericalShaping1.getInitialAzimuthAngle();
                stepSize = ( currentTravelledAngle ) / 100.0;

                // trajectory and trust
                Eigen::Vector6d currentStateVector;
                peakAccelerationSpherical1      = 0.0;
                // Check that the trajectory is feasible, ie curved toward the central body.
                for ( int i = 0 ; i <= 100 ; i++ )
                {
                        double currentThetaAngle = sphericalShaping1.getInitialAzimuthAngle() + i * stepSize;

                        currentStateVector = sphericalShaping1.computeCurrentStateVector(currentThetaAngle);

                        if ( sphericalShaping1.computeCurrentThrustAccelerationVector( currentThetaAngle ).norm() >  peakAccelerationSpherical1 )
                        {
                            peakAccelerationSpherical1 = sphericalShaping1.computeCurrentThrustAccelerationVector( currentThetaAngle ).norm();

                        }

                        fileTrajectories << Nruns << ' ' << 2 << ' ' << Nrev << ' ' << startDate/tudat::physical_constants::JULIAN_DAY << ' ' << TOFt << ' ' << TOF1 << ' ' << TOFC << ' ' << radius << ' ' << azimuth << ' ' << elevation << ' ' << Vr << ' ' << Vtheta << ' ' << Vphi << ' ' << 1 << ' ' << currentStateVector[0] << ' ' << currentStateVector[1] << ' ' << currentStateVector[2] << ' ' << currentStateVector[3] << ' ' << currentStateVector[4] << ' '<< currentStateVector[5] << ' ' << sphericalShaping1.convertIndependentVariableToTime(currentThetaAngle)*tudat::physical_constants::JULIAN_YEAR << ' ' << sphericalShaping1.computeCurrentThrustAccelerationVector( currentThetaAngle ).norm() << '\n';
                }

                sphericalShapingInfeasibleTOF1 = false;
            } // if SS Trajectories

            if (IPTrajectories)
            {
                 // Compute shaped trajectory. - invpolyshaping
                shape_based_methods::InvPolyShaping invPolyShaping1 = shape_based_methods::InvPolyShaping(
                            initialState, coastingState1, TOF1 * tudat::physical_constants::JULIAN_DAY,
                            Nrev, bodyMap, "Vehicle", "Sun", 0.928,
                            rootFinderSettings, -1.0e-3, 1.0, integratorSettingsSB1 );

                 deltaVinvPoly1   = invPolyShaping1.computeDeltaV();

                 // trajectory and thrust
                 currentTravelledAngle = invPolyShaping1.getTravelledAzimuthAngle();
                 stepSize = ( currentTravelledAngle ) / 100.0;
                 // Check that the trajectory is feasible, ie curved toward the central body.
                 peakAccelerationInvPoly1 = 0.0;
                 for ( int i = 0 ; i <= 100 ; i++ )
                     {
                         double currentThetaAngle = i * stepSize;

                         currentStateVector = invPolyShaping1.computeCurrentStateVector(currentThetaAngle);

                         if ( (std::abs(invPolyShaping1.computeCurrentThrustAccelerationMagnitude( currentThetaAngle ))  * physical_constants::ASTRONOMICAL_UNIT / std::pow(tudat::physical_constants::JULIAN_YEAR, 2.0) ) >  peakAccelerationInvPoly1 )
                         {
                             // thrust is non-dimensional so dimensionalisation is needed
                             peakAccelerationInvPoly1 = std::abs(invPolyShaping1.computeCurrentThrustAccelerationMagnitude( currentThetaAngle )) * physical_constants::ASTRONOMICAL_UNIT / std::pow(tudat::physical_constants::JULIAN_YEAR, 2.0);
                         }

                         fileTrajectories << Nruns << ' ' << 1 << ' ' << Nrev << ' ' << startDate/tudat::physical_constants::JULIAN_DAY << ' ' << TOFt << ' ' << TOF1 << ' ' << TOFC << ' ' << radius << ' ' << azimuth << ' ' << elevation << ' ' << Vr << ' ' << Vtheta << ' ' << Vphi << ' ' << 1 << ' ' << currentStateVector[0] << ' ' << currentStateVector[1] << ' ' << currentStateVector[2] << ' ' << currentStateVector[3] << ' ' << currentStateVector[4] << ' '<< currentStateVector[5] << ' ' << invPolyShaping1.convertIndependentVariableToTime(currentThetaAngle)*tudat::physical_constants::JULIAN_YEAR << ' ' << std::abs(invPolyShaping1.computeCurrentThrustAccelerationMagnitude( currentThetaAngle )) * physical_constants::ASTRONOMICAL_UNIT / std::pow(tudat::physical_constants::JULIAN_YEAR, 2.0) << '\n';
                     } // for loop stateVector

                 invPolyShapingInfeasableTOF1   = false;

            } // if IPTrajectories


//            if (fileTrajectories.is_open())
//            {

//            } // if file sopen
            } // try
            catch (const std::runtime_error& error)
            {
                if (SSTrajectories)
                {
                    sphericalShapingInfeasibleTOF1 = true;
                    deltaVspherical1 = 0.0;
                }

                if (IPTrajectories)
                {
                    invPolyShapingInfeasableTOF1 = true;
                    deltaVinvPoly1 = 0.0;
                }

                // std::cout << "DeltaV - spherical " <<  deltaVspherical1 << std::endl;
            } // catch

            /// ################################################################################################################ //
            /// ////////////////////        2. PHASE 2 - COASTING               ///////////////////////////////////////////////////
            /// ################################################################################################################ //




            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            ///////////////////////             CREATE PROPAGATION SETTINGS            ////////////////////////////////////////////
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            // Set initial conditions for the Asterix satellite that will be propagated in this simulation.
            // The initial conditions are given in Keplerian elements and later on converted to Cartesian
            // elements.

            // Set simulation end epoch.
            const double simulationEndEpoch = TOFC*tudat::physical_constants::JULIAN_DAY;

            // Create propagator settings.
            std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
                    std::make_shared< TranslationalStatePropagatorSettings< double > >
                    ( centralBodies, accelerationModelMap, bodiesToPropagate, coastingState1, simulationEndEpoch );



            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            ///////////////////////             PROPAGATE ORBIT            ////////////////////////////////////////////////////////
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            // Create simulation object and propagate dynamics.
            SingleArcDynamicsSimulator< > dynamicsSimulator( bodyMap, integratorSettings, propagatorSettings );
            std::map< double, Eigen::VectorXd > integrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );



            Eigen::VectorXd coastingState2 = ( --integrationResult.end( ) )->second;



//            fileTrajectories << Nruns << ' ' << 2 << ' ' << Nrev << ' ' << startDate/tudat::physical_constants::JULIAN_DAY << ' ' << TOFt << ' ' << TOF1 << ' ' << TOFC << ' ' << radius << ' ' << azimuth << ' ' << elevation << ' ' << Vr << ' ' << Vtheta << ' ' << Vphi << ' ' << 2 << ' ' << coastingState1[0] << ' ' << coastingState1[1] << ' ' << coastingState1[2] << ' ' << coastingState1[3] << ' ' << coastingState1[4] << ' '<< coastingState1[5] << ' ' << 0 << ' ' << 0 << '\n';

            Eigen::Vector6d currentCoastingStateVector;
            for( std::map< double, Eigen::VectorXd >::const_iterator stateIterator = integrationResult.begin( );
                 stateIterator != integrationResult.end( ); stateIterator++ )
            {
                currentCoastingStateVector = stateIterator->second;
                fileTrajectories << Nruns << ' ' << 2 << ' ' << Nrev << ' ' << startDate/tudat::physical_constants::JULIAN_DAY << ' ' << TOFt << ' ' << TOF1 << ' ' << TOFC << ' ' << radius << ' ' << azimuth << ' ' << elevation << ' ' << Vr << ' ' << Vtheta << ' ' << Vphi << ' ' << 2 << ' ' << currentCoastingStateVector[0] << ' ' << currentCoastingStateVector[1] << ' ' << currentCoastingStateVector[2] << ' ' << currentCoastingStateVector[3] << ' ' << currentCoastingStateVector[4] << ' '<< currentCoastingStateVector[5] << ' ' << stateIterator->first/tudat::physical_constants::JULIAN_DAY << ' ' << 0 << '\n';
//                            std::cout << "currentCoastingState " << stateIterator->first/tudat::physical_constants::JULIAN_DAY << std::endl;

            }
//                        std::cout << "TOFC " << TOFC << std::endl;


//            std::cout << "initialState " << initialState << std::endl;
//            std::cout << "coastingState1 " << coastingState1 << std::endl;
//            std::cout << "coastingState2 " << coastingState2 << std::endl;
//            std::cout << "diff coastingState2 " << (coastingState1-coastingState2) << std::endl;
//            std::cout << "finalState " << finalState << std::endl;

            fileTrajectories << Nruns << ' ' << 2 << ' ' << Nrev << ' ' << startDate/tudat::physical_constants::JULIAN_DAY << ' ' << TOFt << ' ' << TOF1 << ' ' << TOFC << ' ' << radius << ' ' << azimuth << ' ' << elevation << ' ' << Vr << ' ' << Vtheta << ' ' << Vphi << ' ' << 2 << ' ' << coastingState2[0] << ' ' << coastingState2[1] << ' ' << coastingState2[2] << ' ' << coastingState2[3] << ' ' << coastingState2[4] << ' '<< coastingState2[5] << ' ' << TOFC*tudat::physical_constants::JULIAN_DAY << ' ' << 0 << '\n';


//            // Write satellite propagation history to file.
//            input_output::writeDataMapToTextFile( integrationResult,
//                                                  "shapeBasedCoastingPropagationHistory_run2.dat",
//                                                  "C:/tudatBundle/tudat/bin/unit_tests/",
//                                                  "",
//                                                  std::numeric_limits< double >::digits10,
//                                                  std::numeric_limits< double >::digits10,
//                                                  "," );




            /// ################################################################################################################ //
            /// ////////////////////        3. PHASE 3 - POWERED FLIGHT               /////////////////////////////////////////////
            /// ################################################################################################################ //

            peakAccelerationSpherical3   = 0.0;
            try
            {

            // Define integrator settings.
                    std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettingsSB3 =
                            std::make_shared< numerical_integrators::IntegratorSettings< double > > (
                                numerical_integrators::rungeKutta4, 0.0, TOF3 * tudat::physical_constants::JULIAN_DAY / 500.0 );

            if (SSTrajectories)
            {
                // Compute shaped trajectory. - spherical shaping - constants influence not yet adjusted
                shape_based_methods::SphericalShaping sphericalShaping3 = shape_based_methods::SphericalShaping(
                            coastingState2, finalState, TOF3 * tudat::physical_constants::JULIAN_DAY,
                            Nrev, bodyMap, "Vehicle", "Sun", 0.0,
                            rootFinderSettings, -1.0e-1, 1.0e-1, integratorSettingsSB3 );

                sphericalShapingInfeasibleTOF3 = false;

                deltaVspherical3 = sphericalShaping3.computeDeltaV();

                // spherical
                currentTravelledAngle = sphericalShaping3.getFinalAzimuthAngle() - sphericalShaping3.getInitialAzimuthAngle();
                stepSize = ( currentTravelledAngle ) / 100.0;

                // Check that the trajectory is feasible, ie curved toward the central body.
                peakAccelerationSpherical3 = 0.0;
                for ( int i = 0 ; i <= 100 ; i++ )
                {
                        double currentThetaAngle = sphericalShaping3.getInitialAzimuthAngle() + i * stepSize;

                        currentStateVector = sphericalShaping3.computeCurrentStateVector(currentThetaAngle);

                        if ( sphericalShaping3.computeCurrentThrustAccelerationVector( currentThetaAngle ).norm() >  peakAccelerationSpherical3 )
                        {
                            peakAccelerationSpherical3 = sphericalShaping3.computeCurrentThrustAccelerationVector( currentThetaAngle ).norm();

                        }

                        fileTrajectories << Nruns << ' ' << 2 << ' ' << Nrev << ' ' << startDate/tudat::physical_constants::JULIAN_DAY << ' ' << TOFt << ' ' << TOF1 << ' ' << TOFC << ' ' << radius << ' ' << azimuth << ' ' << elevation << ' ' << Vr << ' ' << Vtheta << ' ' << Vphi << ' ' << 3 << ' ' << currentStateVector[0] << ' ' << currentStateVector[1] << ' ' << currentStateVector[2] << ' ' << currentStateVector[3] << ' ' << currentStateVector[4] << ' '<< currentStateVector[5] << ' ' << sphericalShaping3.convertIndependentVariableToTime(currentThetaAngle)*tudat::physical_constants::JULIAN_YEAR << ' ' << sphericalShaping3.computeCurrentThrustAccelerationVector( currentThetaAngle ).norm() << '\n';
                }



            } // SSTRAJECTORIES

            if (IPTrajectories)
            {
                // Compute shaped trajectory. - invpolyshaping
                shape_based_methods::InvPolyShaping invPolyShaping3 = shape_based_methods::InvPolyShaping(
                            coastingState2, finalState, TOF3 * tudat::physical_constants::JULIAN_DAY,
                            Nrev, bodyMap, "Vehicle", "Sun", 0.928,
                            rootFinderSettings, -1.0e-3, 1.0, integratorSettingsSB3 );

                invPolyShapingInfeasableTOF3   = false;

                deltaVinvPoly3   = invPolyShaping3.computeDeltaV();

                currentTravelledAngle = invPolyShaping3.getTravelledAzimuthAngle();
                stepSize = ( currentTravelledAngle ) / 100.0;
                // Check that the trajectory is feasible, ie curved toward the central body.
                peakAccelerationInvPoly3 = 0.0;
                    for ( int i = 0 ; i <= 100 ; i++ )
                    {
                        double currentThetaAngle = i * stepSize;

                        currentStateVector = invPolyShaping3.computeCurrentStateVector(currentThetaAngle);

                        if ( (std::abs(invPolyShaping3.computeCurrentThrustAccelerationMagnitude( currentThetaAngle ))  * physical_constants::ASTRONOMICAL_UNIT / std::pow(tudat::physical_constants::JULIAN_YEAR, 2.0) ) >  peakAccelerationInvPoly3 )
                        {
                            // thrust is non-dimensional so dimensionalisation is needed
                            peakAccelerationInvPoly3 = std::abs(invPolyShaping3.computeCurrentThrustAccelerationMagnitude( currentThetaAngle )) * physical_constants::ASTRONOMICAL_UNIT / std::pow(tudat::physical_constants::JULIAN_YEAR, 2.0);
                        }

                        fileTrajectories << Nruns << ' ' << 1 << ' ' << Nrev << ' ' << startDate/tudat::physical_constants::JULIAN_DAY << ' ' << TOFt << ' ' << TOF1 << ' ' << TOFC << ' ' << radius << ' ' << azimuth << ' ' << elevation << ' ' << Vr << ' ' << Vtheta << ' ' << Vphi << ' ' << 3 << ' ' << currentStateVector[0] << ' ' << currentStateVector[1] << ' ' << currentStateVector[2] << ' ' << currentStateVector[3] << ' ' << currentStateVector[4] << ' '<< currentStateVector[5] << ' ' << invPolyShaping3.convertIndependentVariableToTime(currentThetaAngle)*tudat::physical_constants::JULIAN_YEAR << ' ' << std::abs(invPolyShaping3.computeCurrentThrustAccelerationMagnitude( currentThetaAngle )) * physical_constants::ASTRONOMICAL_UNIT / std::pow(tudat::physical_constants::JULIAN_YEAR, 2.0) << '\n';
                    } // for loop stateVector

            } // IP TRAJECTORIES
//            if (fileTrajectories.is_open())
//            {
//            }
            }

            catch (const std::runtime_error& error)
            {

                if (SSTrajectories)
                {
                    sphericalShapingInfeasibleTOF3 = true;
                    deltaVspherical3 = 0.0;
                }
                if (IPTrajectories)
                {
                    invPolyShapingInfeasableTOF3 = true;
                    deltaVinvPoly3 = 0.0;

                }

            }

            if (SSTrajectories)
            {
                file << Nruns << ' ' << 2 << ' ' << Nrev << ' ' << startDate/tudat::physical_constants::JULIAN_DAY << ' ' << TOFt << ' ' << TOF1 << ' ' << TOFC << ' ' << radius << ' ' << azimuth << ' ' << elevation << ' ' << Vr << ' ' << Vtheta << ' ' << Vphi << ' ' << deltaVspherical1 << ' ' << deltaVspherical3 << ' ' << peakAccelerationSpherical1 << ' ' << peakAccelerationSpherical3 << ' ' << '\n';

            }

            if (IPTrajectories)
            {
                file << Nruns << ' ' << 1 << ' ' << Nrev << ' ' << startDate/tudat::physical_constants::JULIAN_DAY << ' ' << TOFt << ' ' << TOF1 << ' ' << TOFC << ' ' << radius << ' ' << azimuth << ' ' << elevation << ' ' << Vr << ' ' << Vtheta << ' ' << Vphi << ' ' << deltaVinvPoly1 << ' ' << deltaVinvPoly3 << ' ' << peakAccelerationInvPoly1 << ' ' << peakAccelerationInvPoly3 << ' ' << '\n';
            }



            // if deltaVspherical1


    }


    }


    if (multiOptimumTrajectoriesFlag)
    {
        for ( int nOptimum = 0 ; nOptimum < 1000 ; nOptimum++ )
        {

            std::cout << "nPop: " << nOptimum << std::endl;


            Nruns++;
            // TOFs
            double Nrev = 0;

            double startDate = optimumInputs(nOptimum, 0 );

            double TOF1      = optimumInputs(nOptimum, 1 )/tudat::physical_constants::JULIAN_DAY;
            double TOFC      = optimumInputs(nOptimum, 2 )/tudat::physical_constants::JULIAN_DAY;
            double TOF3      = optimumInputs(nOptimum, 3 )/tudat::physical_constants::JULIAN_DAY;

            double radius    = optimumInputs(nOptimum, 4 );
            double azimuth   = optimumInputs(nOptimum, 5 );
            double elevation = optimumInputs(nOptimum, 6 );
            double Vr        = optimumInputs(nOptimum, 7 );
            double Vtheta    = optimumInputs(nOptimum, 8 );
            double Vphi      = optimumInputs(nOptimum, 9 );

            double TOFt = TOF1 + TOFC + TOF3;

            // Tempel-1 State
            Eigen::Vector6d keplerianState;
                    keplerianState(0) = 3.145692355253531*physical_constants::ASTRONOMICAL_UNIT;
                    keplerianState(1) = .5096307949367873;
                    keplerianState(2) = 10.47382730765828*mathematical_constants::PI/180.0;
                    keplerianState(3) = 179.203181927074	*mathematical_constants::PI/180.0;
                    keplerianState(4) = 68.74970463832437*mathematical_constants::PI/180.0;
            double meanMotion = .1766565332681526*mathematical_constants::PI/180.0/physical_constants::JULIAN_DAY;
            double initialEpoch = (2457533.5 - 2400000.5 - 51544.0)*tudat::physical_constants::JULIAN_DAY; // JD -> MJD -> MJD2000


            double initialMeanAnomaly = 347.7086591090512*mathematical_constants::PI/180.0;

            double tofEpoch         = TOFt*tudat::physical_constants::JULIAN_DAY;
            double finalMeanAnomaly = initialMeanAnomaly + meanMotion*(startDate + tofEpoch - initialEpoch);
            double eccentricAnomaly = orbital_element_conversions::convertMeanAnomalyToEccentricAnomaly( keplerianState(1), finalMeanAnomaly );
            double trueAnomaly      = orbital_element_conversions::convertEccentricAnomalyToTrueAnomaly( eccentricAnomaly, keplerianState(1) );
            keplerianState(5)       = trueAnomaly;

            Eigen::Vector6d cartesianState = orbital_element_conversions::convertKeplerianToCartesianElements( keplerianState, centralBodyGravitationalParameter );

            // States
            // Earth
            Eigen::Vector6d initialState = pointerToDepartureBodyEphemeris->getCartesianState( startDate );
            // Tempel 1
            Eigen::Vector6d finalState = cartesianState;



            /// ################################################################################################################ //
            /// ////////////////////        0. COASTING SET-UP          ///////////////////////////////////////////////////////////
            /// ################################################################################################################ //


            Eigen::Vector6d sphericalState;
            sphericalState(0) = radius*physical_constants::ASTRONOMICAL_UNIT; // r
            sphericalState(1) = azimuth; // theta
            sphericalState(2) = elevation; // phi
            sphericalState(3) = Vr; // Vr
            sphericalState(4) = Vtheta; // Vtheta
            sphericalState(5) = Vphi; // Vphi


            Eigen::Vector6d coastingState1 = coordinate_conversions::convertSphericalToCartesianState( sphericalState );


            fileTrajectories << Nruns << ' ' << 2 << ' ' << Nrev << ' ' << startDate/tudat::physical_constants::JULIAN_DAY << ' ' << TOFt << ' ' << TOF1 << ' ' << TOFC << ' ' << radius << ' ' << azimuth << ' ' << elevation << ' ' << Vr << ' ' << Vtheta << ' ' << Vphi << ' ' << 1 << ' ' << initialState[0] << ' ' << initialState[1] << ' ' << initialState[2] << ' ' << initialState[3] << ' ' << initialState[4] << ' '<< initialState[5] << ' ' << 0 << ' ' << 0 << '\n';


            /// ################################################################################################################ //
            /// ////////////////////        1. PHASE 1 - POWERED FLIGHT         ///////////////////////////////////////////////////
            /// ################################################################################################################ //
            bool sphericalShapingInfeasibleTOF1;
            bool sphericalShapingInfeasibleTOF3;

            bool invPolyShapingInfeasableTOF1;
            bool invPolyShapingInfeasableTOF3;

            try
            {

            // Define integrator settings.
                std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettingsSB1 =
                        std::make_shared< numerical_integrators::IntegratorSettings< double > > (
                            numerical_integrators::rungeKutta4, 0.0, TOF1 * tudat::physical_constants::JULIAN_DAY / 500.0 );


            if (SSTrajectories)
            {
                // Compute shaped trajectory. - spherical shaping - constants influence not yet adjusted
                shape_based_methods::SphericalShaping sphericalShaping1 = shape_based_methods::SphericalShaping(
                            initialState, coastingState1, TOF1 * tudat::physical_constants::JULIAN_DAY,
                            Nrev, bodyMap, "Vehicle", "Sun", 0.0,
                            rootFinderSettings, -1.0e-1, 1.0e-1, integratorSettingsSB1 );

                deltaVspherical1 = sphericalShaping1.computeDeltaV();

                currentTravelledAngle = sphericalShaping1.getFinalAzimuthAngle() - sphericalShaping1.getInitialAzimuthAngle();
                stepSize = ( currentTravelledAngle ) / 100.0;

                // trajectory and trust
                Eigen::Vector6d currentStateVector;
                peakAccelerationSpherical1      = 0.0;
                // Check that the trajectory is feasible, ie curved toward the central body.
                for ( int i = 0 ; i <= 100 ; i++ )
                {
                        double currentThetaAngle = sphericalShaping1.getInitialAzimuthAngle() + i * stepSize;

                        currentStateVector = sphericalShaping1.computeCurrentStateVector(currentThetaAngle);

                        if ( sphericalShaping1.computeCurrentThrustAccelerationVector( currentThetaAngle ).norm() >  peakAccelerationSpherical1 )
                        {
                            peakAccelerationSpherical1 = sphericalShaping1.computeCurrentThrustAccelerationVector( currentThetaAngle ).norm();

                        }

                        fileTrajectories << Nruns << ' ' << 2 << ' ' << Nrev << ' ' << startDate/tudat::physical_constants::JULIAN_DAY << ' ' << TOFt << ' ' << TOF1 << ' ' << TOFC << ' ' << radius << ' ' << azimuth << ' ' << elevation << ' ' << Vr << ' ' << Vtheta << ' ' << Vphi << ' ' << 1 << ' ' << currentStateVector[0] << ' ' << currentStateVector[1] << ' ' << currentStateVector[2] << ' ' << currentStateVector[3] << ' ' << currentStateVector[4] << ' '<< currentStateVector[5] << ' ' << sphericalShaping1.convertIndependentVariableToTime(currentThetaAngle)*tudat::physical_constants::JULIAN_YEAR << ' ' << sphericalShaping1.computeCurrentThrustAccelerationVector( currentThetaAngle ).norm() << '\n';
                }

                sphericalShapingInfeasibleTOF1 = false;
            } // if SS Trajectories

            if (IPTrajectories)
            {
                 // Compute shaped trajectory. - invpolyshaping
                shape_based_methods::InvPolyShaping invPolyShaping1 = shape_based_methods::InvPolyShaping(
                            initialState, coastingState1, TOF1 * tudat::physical_constants::JULIAN_DAY,
                            Nrev, bodyMap, "Vehicle", "Sun", 0.928,
                            rootFinderSettings, -1.0e-3, 1.0, integratorSettingsSB1 );

                 deltaVinvPoly1   = invPolyShaping1.computeDeltaV();

                 // trajectory and thrust
                 currentTravelledAngle = invPolyShaping1.getTravelledAzimuthAngle();
                 stepSize = ( currentTravelledAngle ) / 100.0;
                 // Check that the trajectory is feasible, ie curved toward the central body.
                 peakAccelerationInvPoly1 = 0.0;
                 for ( int i = 0 ; i <= 100 ; i++ )
                     {
                         double currentThetaAngle = i * stepSize;

                         currentStateVector = invPolyShaping1.computeCurrentStateVector(currentThetaAngle);

                         if ( (std::abs(invPolyShaping1.computeCurrentThrustAccelerationMagnitude( currentThetaAngle ))  * physical_constants::ASTRONOMICAL_UNIT / std::pow(tudat::physical_constants::JULIAN_YEAR, 2.0) ) >  peakAccelerationInvPoly1 )
                         {
                             // thrust is non-dimensional so dimensionalisation is needed
                             peakAccelerationInvPoly1 = std::abs(invPolyShaping1.computeCurrentThrustAccelerationMagnitude( currentThetaAngle )) * physical_constants::ASTRONOMICAL_UNIT / std::pow(tudat::physical_constants::JULIAN_YEAR, 2.0);
                         }

                         fileTrajectories << Nruns << ' ' << 1 << ' ' << Nrev << ' ' << startDate/tudat::physical_constants::JULIAN_DAY << ' ' << TOFt << ' ' << TOF1 << ' ' << TOFC << ' ' << radius << ' ' << azimuth << ' ' << elevation << ' ' << Vr << ' ' << Vtheta << ' ' << Vphi << ' ' << 1 << ' ' << currentStateVector[0] << ' ' << currentStateVector[1] << ' ' << currentStateVector[2] << ' ' << currentStateVector[3] << ' ' << currentStateVector[4] << ' '<< currentStateVector[5] << ' ' << invPolyShaping1.convertIndependentVariableToTime(currentThetaAngle)*tudat::physical_constants::JULIAN_YEAR << ' ' << std::abs(invPolyShaping1.computeCurrentThrustAccelerationMagnitude( currentThetaAngle )) * physical_constants::ASTRONOMICAL_UNIT / std::pow(tudat::physical_constants::JULIAN_YEAR, 2.0) << '\n';
                     } // for loop stateVector

                 invPolyShapingInfeasableTOF1   = false;

            } // if IPTrajectories


//            if (fileTrajectories.is_open())
//            {

//            } // if file sopen
            } // try
            catch (const std::runtime_error& error)
            {
                if (SSTrajectories)
                {
                    sphericalShapingInfeasibleTOF1 = true;
                    deltaVspherical1 = 0.0;
                }

                if (IPTrajectories)
                {
                    invPolyShapingInfeasableTOF1 = true;
                    deltaVinvPoly1 = 0.0;
                }

                // std::cout << "DeltaV - spherical " <<  deltaVspherical1 << std::endl;
            } // catch

            /// ################################################################################################################ //
            /// ////////////////////        2. PHASE 2 - COASTING               ///////////////////////////////////////////////////
            /// ################################################################################################################ //




            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            ///////////////////////             CREATE PROPAGATION SETTINGS            ////////////////////////////////////////////
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            // Set initial conditions for the Asterix satellite that will be propagated in this simulation.
            // The initial conditions are given in Keplerian elements and later on converted to Cartesian
            // elements.

            // Set simulation end epoch.
            const double simulationEndEpoch = TOFC*tudat::physical_constants::JULIAN_DAY;

            // Create propagator settings.
            std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
                    std::make_shared< TranslationalStatePropagatorSettings< double > >
                    ( centralBodies, accelerationModelMap, bodiesToPropagate, coastingState1, simulationEndEpoch );



            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            ///////////////////////             PROPAGATE ORBIT            ////////////////////////////////////////////////////////
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            // Create simulation object and propagate dynamics.
            SingleArcDynamicsSimulator< > dynamicsSimulator( bodyMap, integratorSettings, propagatorSettings );
            std::map< double, Eigen::VectorXd > integrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );



            Eigen::VectorXd coastingState2 = ( --integrationResult.end( ) )->second;



//            fileTrajectories << Nruns << ' ' << 2 << ' ' << Nrev << ' ' << startDate/tudat::physical_constants::JULIAN_DAY << ' ' << TOFt << ' ' << TOF1 << ' ' << TOFC << ' ' << radius << ' ' << azimuth << ' ' << elevation << ' ' << Vr << ' ' << Vtheta << ' ' << Vphi << ' ' << 2 << ' ' << coastingState1[0] << ' ' << coastingState1[1] << ' ' << coastingState1[2] << ' ' << coastingState1[3] << ' ' << coastingState1[4] << ' '<< coastingState1[5] << ' ' << 0 << ' ' << 0 << '\n';

            Eigen::Vector6d currentCoastingStateVector;
            for( std::map< double, Eigen::VectorXd >::const_iterator stateIterator = integrationResult.begin( );
                 stateIterator != integrationResult.end( ); stateIterator++ )
            {
                currentCoastingStateVector = stateIterator->second;
                fileTrajectories << Nruns << ' ' << 2 << ' ' << Nrev << ' ' << startDate/tudat::physical_constants::JULIAN_DAY << ' ' << TOFt << ' ' << TOF1 << ' ' << TOFC << ' ' << radius << ' ' << azimuth << ' ' << elevation << ' ' << Vr << ' ' << Vtheta << ' ' << Vphi << ' ' << 2 << ' ' << currentCoastingStateVector[0] << ' ' << currentCoastingStateVector[1] << ' ' << currentCoastingStateVector[2] << ' ' << currentCoastingStateVector[3] << ' ' << currentCoastingStateVector[4] << ' '<< currentCoastingStateVector[5] << ' ' << stateIterator->first/tudat::physical_constants::JULIAN_DAY << ' ' << 0 << '\n';
//                            std::cout << "currentCoastingState " << stateIterator->first/tudat::physical_constants::JULIAN_DAY << std::endl;

            }
//                        std::cout << "TOFC " << TOFC << std::endl;


//            std::cout << "initialState " << initialState << std::endl;
//            std::cout << "coastingState1 " << coastingState1 << std::endl;
//            std::cout << "coastingState2 " << coastingState2 << std::endl;
//            std::cout << "diff coastingState2 " << (coastingState1-coastingState2) << std::endl;
//            std::cout << "finalState " << finalState << std::endl;

            fileTrajectories << Nruns << ' ' << 2 << ' ' << Nrev << ' ' << startDate/tudat::physical_constants::JULIAN_DAY << ' ' << TOFt << ' ' << TOF1 << ' ' << TOFC << ' ' << radius << ' ' << azimuth << ' ' << elevation << ' ' << Vr << ' ' << Vtheta << ' ' << Vphi << ' ' << 2 << ' ' << coastingState2[0] << ' ' << coastingState2[1] << ' ' << coastingState2[2] << ' ' << coastingState2[3] << ' ' << coastingState2[4] << ' '<< coastingState2[5] << ' ' << TOFC*tudat::physical_constants::JULIAN_DAY << ' ' << 0 << '\n';


//            // Write satellite propagation history to file.
//            input_output::writeDataMapToTextFile( integrationResult,
//                                                  "shapeBasedCoastingPropagationHistory_run2.dat",
//                                                  "C:/tudatBundle/tudat/bin/unit_tests/",
//                                                  "",
//                                                  std::numeric_limits< double >::digits10,
//                                                  std::numeric_limits< double >::digits10,
//                                                  "," );




            /// ################################################################################################################ //
            /// ////////////////////        3. PHASE 3 - POWERED FLIGHT               /////////////////////////////////////////////
            /// ################################################################################################################ //

            peakAccelerationSpherical3   = 0.0;
            try
            {

            // Define integrator settings.
                    std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettingsSB3 =
                            std::make_shared< numerical_integrators::IntegratorSettings< double > > (
                                numerical_integrators::rungeKutta4, 0.0, TOF3 * tudat::physical_constants::JULIAN_DAY / 500.0 );

            if (SSTrajectories)
            {
                // Compute shaped trajectory. - spherical shaping - constants influence not yet adjusted
                shape_based_methods::SphericalShaping sphericalShaping3 = shape_based_methods::SphericalShaping(
                            coastingState2, finalState, TOF3 * tudat::physical_constants::JULIAN_DAY,
                            Nrev, bodyMap, "Vehicle", "Sun", 0.0,
                            rootFinderSettings, -1.0e-1, 1.0e-1, integratorSettingsSB3 );

                sphericalShapingInfeasibleTOF3 = false;

                deltaVspherical3 = sphericalShaping3.computeDeltaV();

                // spherical
                currentTravelledAngle = sphericalShaping3.getFinalAzimuthAngle() - sphericalShaping3.getInitialAzimuthAngle();
                stepSize = ( currentTravelledAngle ) / 100.0;

                // Check that the trajectory is feasible, ie curved toward the central body.
                peakAccelerationSpherical3 = 0.0;
                for ( int i = 0 ; i <= 100 ; i++ )
                {
                        double currentThetaAngle = sphericalShaping3.getInitialAzimuthAngle() + i * stepSize;

                        currentStateVector = sphericalShaping3.computeCurrentStateVector(currentThetaAngle);

                        if ( sphericalShaping3.computeCurrentThrustAccelerationVector( currentThetaAngle ).norm() >  peakAccelerationSpherical3 )
                        {
                            peakAccelerationSpherical3 = sphericalShaping3.computeCurrentThrustAccelerationVector( currentThetaAngle ).norm();

                        }

                        fileTrajectories << Nruns << ' ' << 2 << ' ' << Nrev << ' ' << startDate/tudat::physical_constants::JULIAN_DAY << ' ' << TOFt << ' ' << TOF1 << ' ' << TOFC << ' ' << radius << ' ' << azimuth << ' ' << elevation << ' ' << Vr << ' ' << Vtheta << ' ' << Vphi << ' ' << 3 << ' ' << currentStateVector[0] << ' ' << currentStateVector[1] << ' ' << currentStateVector[2] << ' ' << currentStateVector[3] << ' ' << currentStateVector[4] << ' '<< currentStateVector[5] << ' ' << sphericalShaping3.convertIndependentVariableToTime(currentThetaAngle)*tudat::physical_constants::JULIAN_YEAR << ' ' << sphericalShaping3.computeCurrentThrustAccelerationVector( currentThetaAngle ).norm() << '\n';
                }



            } // SSTRAJECTORIES

            if (IPTrajectories)
            {
                // Compute shaped trajectory. - invpolyshaping
                shape_based_methods::InvPolyShaping invPolyShaping3 = shape_based_methods::InvPolyShaping(
                            coastingState2, finalState, TOF3 * tudat::physical_constants::JULIAN_DAY,
                            Nrev, bodyMap, "Vehicle", "Sun", 0.928,
                            rootFinderSettings, -1.0e-3, 1.0, integratorSettingsSB3 );

                invPolyShapingInfeasableTOF3   = false;

                deltaVinvPoly3   = invPolyShaping3.computeDeltaV();

                currentTravelledAngle = invPolyShaping3.getTravelledAzimuthAngle();
                stepSize = ( currentTravelledAngle ) / 100.0;
                // Check that the trajectory is feasible, ie curved toward the central body.
                peakAccelerationInvPoly3 = 0.0;
                    for ( int i = 0 ; i <= 100 ; i++ )
                    {
                        double currentThetaAngle = i * stepSize;

                        currentStateVector = invPolyShaping3.computeCurrentStateVector(currentThetaAngle);

                        if ( (std::abs(invPolyShaping3.computeCurrentThrustAccelerationMagnitude( currentThetaAngle ))  * physical_constants::ASTRONOMICAL_UNIT / std::pow(tudat::physical_constants::JULIAN_YEAR, 2.0) ) >  peakAccelerationInvPoly3 )
                        {
                            // thrust is non-dimensional so dimensionalisation is needed
                            peakAccelerationInvPoly3 = std::abs(invPolyShaping3.computeCurrentThrustAccelerationMagnitude( currentThetaAngle )) * physical_constants::ASTRONOMICAL_UNIT / std::pow(tudat::physical_constants::JULIAN_YEAR, 2.0);
                        }

                        fileTrajectories << Nruns << ' ' << 1 << ' ' << Nrev << ' ' << startDate/tudat::physical_constants::JULIAN_DAY << ' ' << TOFt << ' ' << TOF1 << ' ' << TOFC << ' ' << radius << ' ' << azimuth << ' ' << elevation << ' ' << Vr << ' ' << Vtheta << ' ' << Vphi << ' ' << 3 << ' ' << currentStateVector[0] << ' ' << currentStateVector[1] << ' ' << currentStateVector[2] << ' ' << currentStateVector[3] << ' ' << currentStateVector[4] << ' '<< currentStateVector[5] << ' ' << invPolyShaping3.convertIndependentVariableToTime(currentThetaAngle)*tudat::physical_constants::JULIAN_YEAR << ' ' << std::abs(invPolyShaping3.computeCurrentThrustAccelerationMagnitude( currentThetaAngle )) * physical_constants::ASTRONOMICAL_UNIT / std::pow(tudat::physical_constants::JULIAN_YEAR, 2.0) << '\n';
                    } // for loop stateVector

            } // IP TRAJECTORIES
//            if (fileTrajectories.is_open())
//            {
//            }
            }

            catch (const std::runtime_error& error)
            {

                if (SSTrajectories)
                {
                    sphericalShapingInfeasibleTOF3 = true;
                    deltaVspherical3 = 0.0;
                }
                if (IPTrajectories)
                {
                    invPolyShapingInfeasableTOF3 = true;
                    deltaVinvPoly3 = 0.0;

                }

            }

            if (SSTrajectories)
            {
                file << Nruns << ' ' << 2 << ' ' << Nrev << ' ' << startDate/tudat::physical_constants::JULIAN_DAY << ' ' << TOFt << ' ' << TOF1 << ' ' << TOFC << ' ' << radius << ' ' << azimuth << ' ' << elevation << ' ' << Vr << ' ' << Vtheta << ' ' << Vphi << ' ' << deltaVspherical1 << ' ' << deltaVspherical3 << ' ' << peakAccelerationSpherical1 << ' ' << peakAccelerationSpherical3 << ' ' << '\n';

            }

            if (IPTrajectories)
            {
                file << Nruns << ' ' << 1 << ' ' << Nrev << ' ' << startDate/tudat::physical_constants::JULIAN_DAY << ' ' << TOFt << ' ' << TOF1 << ' ' << TOFC << ' ' << radius << ' ' << azimuth << ' ' << elevation << ' ' << Vr << ' ' << Vtheta << ' ' << Vphi << ' ' << deltaVinvPoly1 << ' ' << deltaVinvPoly3 << ' ' << peakAccelerationInvPoly1 << ' ' << peakAccelerationInvPoly3 << ' ' << '\n';
            }



            // if deltaVspherical1


    }


    }


    if (multiIPSSOptimumTrajectoriesFlag)
    {
        for ( int nOptimum = 0 ; nOptimum < 1000 ; nOptimum++ )
        {

            std::cout << "nPop: " << nOptimum << std::endl;


            Nruns++;
            // TOFs
            double Nrev = 0;

            double startDate = optimumInputs(nOptimum, 0 );

            double TOF1      = optimumInputs(nOptimum, 1 )/tudat::physical_constants::JULIAN_DAY;
            double TOFC      = optimumInputs(nOptimum, 2 )/tudat::physical_constants::JULIAN_DAY;
            double TOF3      = optimumInputs(nOptimum, 3 )/tudat::physical_constants::JULIAN_DAY;

            double radius    = optimumInputs(nOptimum, 4 );
            double azimuth   = optimumInputs(nOptimum, 5 );
            double elevation = optimumInputs(nOptimum, 6 );
            double Vr        = optimumInputs(nOptimum, 7 );
            double Vtheta    = optimumInputs(nOptimum, 8 );
            double Vphi      = optimumInputs(nOptimum, 9 );

            double TOFt = TOF1 + TOFC + TOF3;

            // Tempel-1 State
            Eigen::Vector6d keplerianState;
                    keplerianState(0) = 3.145692355253531*physical_constants::ASTRONOMICAL_UNIT;
                    keplerianState(1) = .5096307949367873;
                    keplerianState(2) = 10.47382730765828*mathematical_constants::PI/180.0;
                    keplerianState(3) = 179.203181927074	*mathematical_constants::PI/180.0;
                    keplerianState(4) = 68.74970463832437*mathematical_constants::PI/180.0;
            double meanMotion = .1766565332681526*mathematical_constants::PI/180.0/physical_constants::JULIAN_DAY;
            double initialEpoch = (2457533.5 - 2400000.5 - 51544.0)*tudat::physical_constants::JULIAN_DAY; // JD -> MJD -> MJD2000


            double initialMeanAnomaly = 347.7086591090512*mathematical_constants::PI/180.0;

            double tofEpoch         = TOFt*tudat::physical_constants::JULIAN_DAY;
            double finalMeanAnomaly = initialMeanAnomaly + meanMotion*(startDate + tofEpoch - initialEpoch);
            double eccentricAnomaly = orbital_element_conversions::convertMeanAnomalyToEccentricAnomaly( keplerianState(1), finalMeanAnomaly );
            double trueAnomaly      = orbital_element_conversions::convertEccentricAnomalyToTrueAnomaly( eccentricAnomaly, keplerianState(1) );
            keplerianState(5)       = trueAnomaly;

            Eigen::Vector6d cartesianState = orbital_element_conversions::convertKeplerianToCartesianElements( keplerianState, centralBodyGravitationalParameter );

            // States
            // Earth
            Eigen::Vector6d initialState = pointerToDepartureBodyEphemeris->getCartesianState( startDate );
            // Tempel 1
            Eigen::Vector6d finalState = cartesianState;



            /// ################################################################################################################ //
            /// ////////////////////        0. COASTING SET-UP          ///////////////////////////////////////////////////////////
            /// ################################################################################################################ //


            Eigen::Vector6d sphericalState;
            sphericalState(0) = radius*physical_constants::ASTRONOMICAL_UNIT; // r
            sphericalState(1) = azimuth; // theta
            sphericalState(2) = elevation; // phi
            sphericalState(3) = Vr; // Vr
            sphericalState(4) = Vtheta; // Vtheta
            sphericalState(5) = Vphi; // Vphi


            Eigen::Vector6d coastingState1 = coordinate_conversions::convertSphericalToCartesianState( sphericalState );


            fileTrajectories << Nruns << ' ' << 2 << ' ' << Nrev << ' ' << startDate/tudat::physical_constants::JULIAN_DAY << ' ' << TOFt << ' ' << TOF1 << ' ' << TOFC << ' ' << radius << ' ' << azimuth << ' ' << elevation << ' ' << Vr << ' ' << Vtheta << ' ' << Vphi << ' ' << 1 << ' ' << initialState[0] << ' ' << initialState[1] << ' ' << initialState[2] << ' ' << initialState[3] << ' ' << initialState[4] << ' '<< initialState[5] << ' ' << 0 << ' ' << 0 << '\n';


            /// ################################################################################################################ //
            /// ////////////////////        1. PHASE 1 - POWERED FLIGHT         ///////////////////////////////////////////////////
            /// ################################################################################################################ //
            bool sphericalShapingInfeasibleTOF1;
            bool sphericalShapingInfeasibleTOF3;

            bool invPolyShapingInfeasableTOF1;
            bool invPolyShapingInfeasableTOF3;

            try
            {

            // Define integrator settings.
                std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettingsSB1 =
                        std::make_shared< numerical_integrators::IntegratorSettings< double > > (
                            numerical_integrators::rungeKutta4, 0.0, TOF1 * tudat::physical_constants::JULIAN_DAY / 500.0 );


            if (mIPTrajectories)
            {
                 // Compute shaped trajectory. - invpolyshaping
                shape_based_methods::InvPolyShaping invPolyShaping1 = shape_based_methods::InvPolyShaping(
                            initialState, coastingState1, TOF1 * tudat::physical_constants::JULIAN_DAY,
                            Nrev, bodyMap, "Vehicle", "Sun", 0.928,
                            rootFinderSettings, -1.0e-3, 1.0, integratorSettingsSB1 );

                 deltaVinvPoly1   = invPolyShaping1.computeDeltaV();

                 // trajectory and thrust
                 currentTravelledAngle = invPolyShaping1.getTravelledAzimuthAngle();
                 stepSize = ( currentTravelledAngle ) / 100.0;
                 // Check that the trajectory is feasible, ie curved toward the central body.
                 peakAccelerationInvPoly1 = 0.0;
                 for ( int i = 0 ; i <= 100 ; i++ )
                     {
                         double currentThetaAngle = i * stepSize;

                         currentStateVector = invPolyShaping1.computeCurrentStateVector(currentThetaAngle);

                         if ( (std::abs(invPolyShaping1.computeCurrentThrustAccelerationMagnitude( currentThetaAngle ))  * physical_constants::ASTRONOMICAL_UNIT / std::pow(tudat::physical_constants::JULIAN_YEAR, 2.0) ) >  peakAccelerationInvPoly1 )
                         {
                             // thrust is non-dimensional so dimensionalisation is needed
                             peakAccelerationInvPoly1 = std::abs(invPolyShaping1.computeCurrentThrustAccelerationMagnitude( currentThetaAngle )) * physical_constants::ASTRONOMICAL_UNIT / std::pow(tudat::physical_constants::JULIAN_YEAR, 2.0);
                         }

                         fileTrajectories << Nruns << ' ' << 1 << ' ' << Nrev << ' ' << startDate/tudat::physical_constants::JULIAN_DAY << ' ' << TOFt << ' ' << TOF1 << ' ' << TOFC << ' ' << radius << ' ' << azimuth << ' ' << elevation << ' ' << Vr << ' ' << Vtheta << ' ' << Vphi << ' ' << 1 << ' ' << currentStateVector[0] << ' ' << currentStateVector[1] << ' ' << currentStateVector[2] << ' ' << currentStateVector[3] << ' ' << currentStateVector[4] << ' '<< currentStateVector[5] << ' ' << invPolyShaping1.convertIndependentVariableToTime(currentThetaAngle)*tudat::physical_constants::JULIAN_YEAR << ' ' << std::abs(invPolyShaping1.computeCurrentThrustAccelerationMagnitude( currentThetaAngle )) * physical_constants::ASTRONOMICAL_UNIT / std::pow(tudat::physical_constants::JULIAN_YEAR, 2.0) << '\n';
                     } // for loop stateVector

                 invPolyShapingInfeasableTOF1   = false;

            } // if IPTrajectories


//            if (fileTrajectories.is_open())
//            {

//            } // if file sopen
            } // try
            catch (const std::runtime_error& error)
            {
                if (mIPTrajectories)
                {
                    invPolyShapingInfeasableTOF1 = true;
                    deltaVinvPoly1 = 0.0;
                }

                // std::cout << "DeltaV - spherical " <<  deltaVspherical1 << std::endl;
            } // catch

            /// ################################################################################################################ //
            /// ////////////////////        2. PHASE 2 - COASTING               ///////////////////////////////////////////////////
            /// ################################################################################################################ //




            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            ///////////////////////             CREATE PROPAGATION SETTINGS            ////////////////////////////////////////////
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            // Set initial conditions for the Asterix satellite that will be propagated in this simulation.
            // The initial conditions are given in Keplerian elements and later on converted to Cartesian
            // elements.

            // Set simulation end epoch.
            const double simulationEndEpoch = TOFC*tudat::physical_constants::JULIAN_DAY;

            // Create propagator settings.
            std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
                    std::make_shared< TranslationalStatePropagatorSettings< double > >
                    ( centralBodies, accelerationModelMap, bodiesToPropagate, coastingState1, simulationEndEpoch );



            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            ///////////////////////             PROPAGATE ORBIT            ////////////////////////////////////////////////////////
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            // Create simulation object and propagate dynamics.
            SingleArcDynamicsSimulator< > dynamicsSimulator( bodyMap, integratorSettings, propagatorSettings );
            std::map< double, Eigen::VectorXd > integrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );



            Eigen::VectorXd coastingState2 = ( --integrationResult.end( ) )->second;



//            fileTrajectories << Nruns << ' ' << 2 << ' ' << Nrev << ' ' << startDate/tudat::physical_constants::JULIAN_DAY << ' ' << TOFt << ' ' << TOF1 << ' ' << TOFC << ' ' << radius << ' ' << azimuth << ' ' << elevation << ' ' << Vr << ' ' << Vtheta << ' ' << Vphi << ' ' << 2 << ' ' << coastingState1[0] << ' ' << coastingState1[1] << ' ' << coastingState1[2] << ' ' << coastingState1[3] << ' ' << coastingState1[4] << ' '<< coastingState1[5] << ' ' << 0 << ' ' << 0 << '\n';

            Eigen::Vector6d currentCoastingStateVector;
            for( std::map< double, Eigen::VectorXd >::const_iterator stateIterator = integrationResult.begin( );
                 stateIterator != integrationResult.end( ); stateIterator++ )
            {
                currentCoastingStateVector = stateIterator->second;
                fileTrajectories << Nruns << ' ' << 2 << ' ' << Nrev << ' ' << startDate/tudat::physical_constants::JULIAN_DAY << ' ' << TOFt << ' ' << TOF1 << ' ' << TOFC << ' ' << radius << ' ' << azimuth << ' ' << elevation << ' ' << Vr << ' ' << Vtheta << ' ' << Vphi << ' ' << 2 << ' ' << currentCoastingStateVector[0] << ' ' << currentCoastingStateVector[1] << ' ' << currentCoastingStateVector[2] << ' ' << currentCoastingStateVector[3] << ' ' << currentCoastingStateVector[4] << ' '<< currentCoastingStateVector[5] << ' ' << stateIterator->first/tudat::physical_constants::JULIAN_DAY << ' ' << 0 << '\n';
//                            std::cout << "currentCoastingState " << stateIterator->first/tudat::physical_constants::JULIAN_DAY << std::endl;

            }
//                        std::cout << "TOFC " << TOFC << std::endl;


//            std::cout << "initialState " << initialState << std::endl;
//            std::cout << "coastingState1 " << coastingState1 << std::endl;
//            std::cout << "coastingState2 " << coastingState2 << std::endl;
//            std::cout << "diff coastingState2 " << (coastingState1-coastingState2) << std::endl;
//            std::cout << "finalState " << finalState << std::endl;

            fileTrajectories << Nruns << ' ' << 2 << ' ' << Nrev << ' ' << startDate/tudat::physical_constants::JULIAN_DAY << ' ' << TOFt << ' ' << TOF1 << ' ' << TOFC << ' ' << radius << ' ' << azimuth << ' ' << elevation << ' ' << Vr << ' ' << Vtheta << ' ' << Vphi << ' ' << 2 << ' ' << coastingState2[0] << ' ' << coastingState2[1] << ' ' << coastingState2[2] << ' ' << coastingState2[3] << ' ' << coastingState2[4] << ' '<< coastingState2[5] << ' ' << TOFC*tudat::physical_constants::JULIAN_DAY << ' ' << 0 << '\n';


//            // Write satellite propagation history to file.
//            input_output::writeDataMapToTextFile( integrationResult,
//                                                  "shapeBasedCoastingPropagationHistory_run2.dat",
//                                                  "C:/tudatBundle/tudat/bin/unit_tests/",
//                                                  "",
//                                                  std::numeric_limits< double >::digits10,
//                                                  std::numeric_limits< double >::digits10,
//                                                  "," );




            /// ################################################################################################################ //
            /// ////////////////////        3. PHASE 3 - POWERED FLIGHT               /////////////////////////////////////////////
            /// ################################################################################################################ //

            peakAccelerationSpherical3   = 0.0;
            try
            {

            // Define integrator settings.
                    std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettingsSB3 =
                            std::make_shared< numerical_integrators::IntegratorSettings< double > > (
                                numerical_integrators::rungeKutta4, 0.0, TOF3 * tudat::physical_constants::JULIAN_DAY / 500.0 );

            if (mSSTrajectories)
            {
                // Compute shaped trajectory. - spherical shaping - constants influence not yet adjusted
                shape_based_methods::SphericalShaping sphericalShaping3 = shape_based_methods::SphericalShaping(
                            coastingState2, finalState, TOF3 * tudat::physical_constants::JULIAN_DAY,
                            Nrev, bodyMap, "Vehicle", "Sun", 0.0,
                            rootFinderSettings, -1.0e-1, 1.0e-1, integratorSettingsSB3 );

                sphericalShapingInfeasibleTOF3 = false;

                deltaVspherical3 = sphericalShaping3.computeDeltaV();

                // spherical
                currentTravelledAngle = sphericalShaping3.getFinalAzimuthAngle() - sphericalShaping3.getInitialAzimuthAngle();
                stepSize = ( currentTravelledAngle ) / 100.0;

                // Check that the trajectory is feasible, ie curved toward the central body.
                peakAccelerationSpherical3 = 0.0;
                for ( int i = 0 ; i <= 100 ; i++ )
                {
                        double currentThetaAngle = sphericalShaping3.getInitialAzimuthAngle() + i * stepSize;

                        currentStateVector = sphericalShaping3.computeCurrentStateVector(currentThetaAngle);

                        if ( sphericalShaping3.computeCurrentThrustAccelerationVector( currentThetaAngle ).norm() >  peakAccelerationSpherical3 )
                        {
                            peakAccelerationSpherical3 = sphericalShaping3.computeCurrentThrustAccelerationVector( currentThetaAngle ).norm();

                        }

                        fileTrajectories << Nruns << ' ' << 2 << ' ' << Nrev << ' ' << startDate/tudat::physical_constants::JULIAN_DAY << ' ' << TOFt << ' ' << TOF1 << ' ' << TOFC << ' ' << radius << ' ' << azimuth << ' ' << elevation << ' ' << Vr << ' ' << Vtheta << ' ' << Vphi << ' ' << 3 << ' ' << currentStateVector[0] << ' ' << currentStateVector[1] << ' ' << currentStateVector[2] << ' ' << currentStateVector[3] << ' ' << currentStateVector[4] << ' '<< currentStateVector[5] << ' ' << sphericalShaping3.convertIndependentVariableToTime(currentThetaAngle)*tudat::physical_constants::JULIAN_YEAR << ' ' << sphericalShaping3.computeCurrentThrustAccelerationVector( currentThetaAngle ).norm() << '\n';
                }



            } // SSTRAJECTORIES

//            if (fileTrajectories.is_open())
//            {
//            }
            }

            catch (const std::runtime_error& error)
            {
                if (mSSTrajectories)
                {
                    sphericalShapingInfeasibleTOF3 = true;
                    deltaVspherical3 = 0.0;
                }
            }

            if (mSSTrajectories)
            {
                file << Nruns << ' ' << 2 << ' ' << Nrev << ' ' << startDate/tudat::physical_constants::JULIAN_DAY << ' ' << TOFt << ' ' << TOF1 << ' ' << TOFC << ' ' << radius << ' ' << azimuth << ' ' << elevation << ' ' << Vr << ' ' << Vtheta << ' ' << Vphi << ' ' << deltaVinvPoly1 << ' ' << deltaVspherical3 << ' ' << peakAccelerationInvPoly1 << ' ' << peakAccelerationSpherical3 << ' ' << '\n';

            }


            // if deltaVspherical1


    }


    }


}






BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
