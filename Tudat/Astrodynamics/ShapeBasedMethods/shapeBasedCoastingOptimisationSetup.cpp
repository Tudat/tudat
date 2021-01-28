/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "shapeBasedCoastingOptimisationSetup.h"


using namespace tudat::ephemerides;
using namespace tudat::basic_astrodynamics;
using namespace tudat::orbital_element_conversions;
using namespace tudat::basic_mathematics;
using namespace tudat::input_output;
using namespace tudat::transfer_trajectories;

using namespace pagmo;

namespace tudat
{
namespace shape_based_methods
{



// Calculates the fitness
std::vector< double > AnvPolyShapingOptimisationProblem::fitness( const std::vector< double > &x ) const
{

    double departureTime = x.at( 0 );
    double timeOfFlight = x.at( 1 );
    double arrivalTime = departureTime + timeOfFlight;

    double centralBodyGravitationalParameter = bodyMap_["Sun"]->getGravityFieldModel()->getGravitationalParameter( );

    Eigen::Vector6d departureState = initialStateFunction_( departureTime );
    Eigen::Vector6d arrivalState = finalStateFunction_( arrivalTime, centralBodyGravitationalParameter );

    std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings =
            std::make_shared< numerical_integrators::IntegratorSettings< double > > (
                numerical_integrators::rungeKutta4, 0.0, timeOfFlight / 500.0 );


    // creating shaping objects
    InvPolyShaping invPolyShaping     = InvPolyShaping(departureState,arrivalState,timeOfFlight,numberOfRevolutions_,bodyMap_,bodyToPropagate_,centralBody_,initialValueFreeCoefficient_,rootFinderSettings_,lowerBoundFreeCoefficient_,upperBoundFreeCoefficient_,integratorSettings);


    // computing and comparing fitness values
    double deltaVinvpoly;
    bool infeasibleTOF = invPolyShaping.getInfeasibleTOF();
    if (infeasibleTOF)
    {
        deltaVinvpoly = 1e21;
    }
    else
    {
        deltaVinvpoly = invPolyShaping.computeDeltaV();
    }

    //storing fitnessvalues
    std::vector< double > fitnessVector;

    fitnessVector.push_back(deltaVinvpoly);



    return fitnessVector;
} // Invpoly Optimisation

// Calculates the fitness
std::vector< double > SphericalShapingOptimisationProblem::fitness( const std::vector< double > &x ) const
{

    double departureTime = x.at( 0 );
    double timeOfFlight = x.at( 1 );

    double arrivalTime = departureTime + timeOfFlight;

    double centralBodyGravitationalParameter = bodyMap_["Sun"]->getGravityFieldModel()->getGravitationalParameter( );

    Eigen::Vector6d departureState = initialStateFunction_( departureTime );
    Eigen::Vector6d arrivalState = finalStateFunction_( arrivalTime, centralBodyGravitationalParameter );

    std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings =
            std::make_shared< numerical_integrators::IntegratorSettings< double > > (
                numerical_integrators::rungeKutta4, 0.0, timeOfFlight / 500.0 );


    double deltaVspherical;
    try
    {
        SphericalShaping sphericalShaping = SphericalShaping(departureState,arrivalState,timeOfFlight,numberOfRevolutions_, bodyMap_,bodyToPropagate_,centralBody_,initialValueFreeCoefficient_,rootFinderSettings_,lowerBoundFreeCoefficient_,upperBoundFreeCoefficient_,integratorSettings );


        bool infeasibleTOF = sphericalShaping.getInfeasibleTOF();
        if (infeasibleTOF)
        {
            deltaVspherical = 1e22;
        }
        else
        {
            deltaVspherical = sphericalShaping.computeDeltaV();
        }

    }
    catch (const std::runtime_error& error)
    {
        deltaVspherical = 1e22;
    }

    //storing fitnessvalues
    std::vector< double > fitnessVector;

    fitnessVector.push_back(deltaVspherical);



    return fitnessVector;
} // Spherical Optimisation

// Calculates the fitness
std::vector< double > ShapeBasedCoastingOptimisationProblem::fitness( const std::vector< double > &x ) const
{
    using namespace tudat;
    using namespace tudat::simulation_setup;
    using namespace tudat::propagators;
    using namespace tudat::numerical_integrators;
    using namespace tudat::basic_astrodynamics;
    using namespace tudat::basic_mathematics;
    using namespace tudat::orbital_element_conversions;
    using namespace tudat::unit_conversions;
    using namespace tudat::input_output;


    // inputs
    double departureTime = x.at( 0 );

    double timeOfFlight1  = x.at( 1 );
    double timeOfFlightC  = x.at( 2 );
    double timeOfFlight3  = x.at( 3 );

    double radius            = x.at( 4 );
    double azimuth           = x.at( 5 );
    double elevation         = x.at( 6 );
    double radialVelocity    = x.at( 7 );
    double azimuthalVelocity = x.at( 8 );
    double phiVelocity       = x.at( 9 );

    double timeOfFlight = timeOfFlight1 + timeOfFlightC + timeOfFlight3;

    double arrivalTime = departureTime + timeOfFlight;

    double centralBodyGravitationalParameter = bodyMap_["Sun"]->getGravityFieldModel()->getGravitationalParameter( );

    Eigen::Vector6d departureState = initialStateFunction_( departureTime );
    Eigen::Vector6d arrivalState = finalStateFunction_( arrivalTime, centralBodyGravitationalParameter );

    std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings1 =
            std::make_shared< numerical_integrators::IntegratorSettings< double > > (
                numerical_integrators::rungeKutta4, 0.0, timeOfFlight1 / 500.0 );

    std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings3 =
            std::make_shared< numerical_integrators::IntegratorSettings< double > > (
                numerical_integrators::rungeKutta4, 0.0, timeOfFlight3 / 500.0 );



//    Eigen::Vector6d initialStateSphericalCoordinates = coordinate_conversions::convertCartesianToSphericalState( departureState );
//    Eigen::Vector6d finalStateSphericalCoordinates = coordinate_conversions::convertCartesianToSphericalState( arrivalState );

//     coastingState 1
    Eigen::Vector6d sphericalState1 ;
    sphericalState1(0) = radius*physical_constants::ASTRONOMICAL_UNIT; // r
    sphericalState1(1) = azimuth;                                      // theta
    sphericalState1(2) = elevation;                                    // elevation
    sphericalState1(3) = radialVelocity;//*initialStateSphericalCoordinates(3);                               // Vr
    sphericalState1(4) = azimuthalVelocity;//*initialStateSphericalCoordinates(4);                            // Vtheta
    sphericalState1(5) = phiVelocity;//*initialStateSphericalCoordinates(5);                                  // Vphi
    Eigen::Vector6d coastingState1 = coordinate_conversions::convertSphericalToCartesianState( sphericalState1 );
//    Eigen::Vector6d coastingState1 = arrivalState;



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
                bodyMap_, accelerationMap, bodiesToPropagate, centralBodies );


    double simulationStartEpoch = 0.0;
    const double fixedStepSize = 24*3600.0;
    std::shared_ptr< IntegratorSettings< > > integratorSettingsC =
            std::make_shared< IntegratorSettings< > >
            ( rungeKutta4, simulationStartEpoch, fixedStepSize );



    // Phase 1
    double deltaV1;
    try
    {
//        SphericalShaping shapeBasedMethod1 = SphericalShaping(departureState,coastingState1,timeOfFlight1,numberOfRevolutions_, bodyMap_,bodyToPropagate_,centralBody_,initialValueFreeCoefficient_,rootFinderSettings_,lowerBoundFreeCoefficient_,upperBoundFreeCoefficient_,integratorSettings1 );
        InvPolyShaping   shapeBasedMethod1 = InvPolyShaping(departureState,coastingState1,timeOfFlight1,numberOfRevolutions_,bodyMap_,bodyToPropagate_,centralBody_,initialValueFreeCoefficient_,rootFinderSettings_,lowerBoundFreeCoefficient_,upperBoundFreeCoefficient_,integratorSettings1);

        // spherical shaping
        bool infeasibleTOF = shapeBasedMethod1.getInfeasibleTOF();
        if (infeasibleTOF)
        {
            deltaV1 = 1e22;
        }
        else
        {
            deltaV1 = shapeBasedMethod1.computeDeltaV();
        }


    }
    catch (const std::runtime_error& error)
    {
        deltaV1 = 1e22;
    }


    double deltaV3;
    if (deltaV1!=1e22)
    {
        // Phase 2 - Coasting
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////             CREATE PROPAGATION SETTINGS            ////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // Set initial conditions for the Asterix satellite that will be propagated in this simulation.
        // The initial conditions are given in Keplerian elements and later on converted to Cartesian
        // elements.

        // Set simulation end epoch.
        const double simulationEndEpoch = timeOfFlightC;

        // Create propagator settings.
        std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
                std::make_shared< TranslationalStatePropagatorSettings< double > >
                ( centralBodies, accelerationModelMap, bodiesToPropagate, coastingState1, simulationEndEpoch );



        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////             PROPAGATE ORBIT            ////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // Create simulation object and propagate dynamics.
        SingleArcDynamicsSimulator< > dynamicsSimulator( bodyMap_, integratorSettingsC, propagatorSettings );
        std::map< double, Eigen::VectorXd > integrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );

        Eigen::Vector6d coastingState2 = ( --integrationResult.end( ) )->second;
    //    Eigen::Vector6d coastingState2 = departureState;


        // Phase 3 - powered

        try
        {
//            SphericalShaping shapeBasedMethod3 = SphericalShaping(coastingState2,arrivalState,timeOfFlight3,numberOfRevolutions_, bodyMap_,bodyToPropagate_,centralBody_,initialValueFreeCoefficient_,rootFinderSettings_,lowerBoundFreeCoefficient_,upperBoundFreeCoefficient_,integratorSettings3 );
            InvPolyShaping   shapeBasedMethod3 = InvPolyShaping(coastingState2,arrivalState,timeOfFlight3,numberOfRevolutions_,bodyMap_,bodyToPropagate_,centralBody_,initialValueFreeCoefficient_,rootFinderSettings_,lowerBoundFreeCoefficient_,upperBoundFreeCoefficient_,integratorSettings3 );


            bool infeasibleTOF = shapeBasedMethod3.getInfeasibleTOF();
            if (infeasibleTOF)
            {
                deltaV3 = 1e22;
            }
            else
            {
                deltaV3 = shapeBasedMethod3.computeDeltaV();
            }

        }
        catch (const std::runtime_error& error)
        {
            deltaV3 = 1e22;
        }
    } // if dv spherical !=1e22
    else
    {
        deltaV3 = 1e22;
    }

    double deltaV = deltaV1 + deltaV3;

//    double deltaVspherical = 20;
    //storing fitnessvalues
    std::vector< double > fitnessVector;

    fitnessVector.push_back(deltaV);

//    std::cout << "----------------------------------------------" << std::endl;
//    std::cout << "1 pop" << std::endl;
//    std::cout << "----------------------------------------------" << std::endl;


    return fitnessVector;
} // Optimisation

std::vector< double > ShapeBasedIPCoastingMultiOptimisationProblem::fitness( const std::vector< double > &x ) const
{
    using namespace tudat;
    using namespace tudat::simulation_setup;
    using namespace tudat::propagators;
    using namespace tudat::numerical_integrators;
    using namespace tudat::basic_astrodynamics;
    using namespace tudat::basic_mathematics;
    using namespace tudat::orbital_element_conversions;
    using namespace tudat::unit_conversions;
    using namespace tudat::input_output;


    // inputs
    double departureTime = x.at( 0 );

    double timeOfFlight1  = x.at( 1 );
    double timeOfFlightC  = x.at( 2 );
    double timeOfFlight3  = x.at( 3 );

    double radius            = x.at( 4 );
    double azimuth           = x.at( 5 );
    double elevation         = x.at( 6 );
    double radialVelocity    = x.at( 7 );
    double azimuthalVelocity = x.at( 8 );
    double phiVelocity       = x.at( 9 );

    double timeOfFlight = timeOfFlight1 + timeOfFlightC + timeOfFlight3;

    double arrivalTime = departureTime + timeOfFlight;

    double centralBodyGravitationalParameter = bodyMap_["Sun"]->getGravityFieldModel()->getGravitationalParameter( );

    Eigen::Vector6d departureState = initialStateFunction_( departureTime );
    Eigen::Vector6d arrivalState = finalStateFunction_( arrivalTime, centralBodyGravitationalParameter );

    std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings1 =
            std::make_shared< numerical_integrators::IntegratorSettings< double > > (
                numerical_integrators::rungeKutta4, 0.0, timeOfFlight1 / 500.0 );

    std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings3 =
            std::make_shared< numerical_integrators::IntegratorSettings< double > > (
                numerical_integrators::rungeKutta4, 0.0, timeOfFlight3 / 500.0 );



//    Eigen::Vector6d initialStateSphericalCoordinates = coordinate_conversions::convertCartesianToSphericalState( departureState );
//    Eigen::Vector6d finalStateSphericalCoordinates = coordinate_conversions::convertCartesianToSphericalState( arrivalState );

//     coastingState 1
    Eigen::Vector6d sphericalState1 ;
    sphericalState1(0) = radius*physical_constants::ASTRONOMICAL_UNIT; // r
    sphericalState1(1) = azimuth;                                      // theta
    sphericalState1(2) = elevation;                                    // elevation
    sphericalState1(3) = radialVelocity;//*initialStateSphericalCoordinates(3);                               // Vr
    sphericalState1(4) = azimuthalVelocity;//*initialStateSphericalCoordinates(4);                            // Vtheta
    sphericalState1(5) = phiVelocity;//*initialStateSphericalCoordinates(5);                                  // Vphi
    Eigen::Vector6d coastingState1 = coordinate_conversions::convertSphericalToCartesianState( sphericalState1 );
//    Eigen::Vector6d coastingState1 = arrivalState;



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
                bodyMap_, accelerationMap, bodiesToPropagate, centralBodies );


    double simulationStartEpoch = 0.0;
    const double fixedStepSize = 24*3600.0;
    std::shared_ptr< IntegratorSettings< > > integratorSettingsC =
            std::make_shared< IntegratorSettings< > >
            ( rungeKutta4, simulationStartEpoch, fixedStepSize );



    double currentTravelledAngle;
    double stepSize;

    // Phase 1
    double deltaV1;
    double peakAcceleration1;
    try
    {
//        SphericalShaping shapeBasedMethod1 = SphericalShaping(departureState,coastingState1,timeOfFlight1,numberOfRevolutions_, bodyMap_,bodyToPropagate_,centralBody_,initialValueFreeCoefficient_,rootFinderSettings_,lowerBoundFreeCoefficient_,upperBoundFreeCoefficient_,integratorSettings1 );
        InvPolyShaping   shapeBasedMethod1 = InvPolyShaping(departureState,coastingState1,timeOfFlight1,numberOfRevolutions_,bodyMap_,bodyToPropagate_,centralBody_,initialValueFreeCoefficient_,rootFinderSettings_,lowerBoundFreeCoefficient_,upperBoundFreeCoefficient_,integratorSettings1);

        // shaping
        bool infeasibleTOF = shapeBasedMethod1.getInfeasibleTOF();
        if (infeasibleTOF)
        {
            deltaV1           = 1e22;
            peakAcceleration1 = 1e22;
        }
        else
        {
            deltaV1 = shapeBasedMethod1.computeDeltaV();

            // Peak Acceleration
            // trajectory and thrust
            currentTravelledAngle = shapeBasedMethod1.getTravelledAzimuthAngle();
            stepSize = ( currentTravelledAngle ) / 100.0;
            // Check that the trajectory is feasible, ie curved toward the central body.
            peakAcceleration1 = 0.0;
            for ( int i = 0 ; i <= 100 ; i++ )
            {
                double currentThetaAngle = i * stepSize;
                if ( (std::abs(shapeBasedMethod1.computeCurrentThrustAccelerationMagnitude( currentThetaAngle ))  * physical_constants::ASTRONOMICAL_UNIT / std::pow(tudat::physical_constants::JULIAN_YEAR, 2.0) ) >  peakAcceleration1 )
                {
                    // thrust is non-dimensional so dimensionalisation is needed
                    peakAcceleration1 = std::abs(shapeBasedMethod1.computeCurrentThrustAccelerationMagnitude( currentThetaAngle )) * physical_constants::ASTRONOMICAL_UNIT / std::pow(tudat::physical_constants::JULIAN_YEAR, 2.0);
                }
            } // for loop stateVector
        }

    }
    catch (const std::runtime_error& error)
    {
        deltaV1 = 1e22;
    }


    double deltaV3;
    double peakAcceleration3;

    if (deltaV1!=1e22)
    {
        // Phase 2 - Coasting
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////             CREATE PROPAGATION SETTINGS            ////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // Set initial conditions for the Asterix satellite that will be propagated in this simulation.
        // The initial conditions are given in Keplerian elements and later on converted to Cartesian
        // elements.

        // Set simulation end epoch.
        const double simulationEndEpoch = timeOfFlightC;

        // Create propagator settings.
        std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
                std::make_shared< TranslationalStatePropagatorSettings< double > >
                ( centralBodies, accelerationModelMap, bodiesToPropagate, coastingState1, simulationEndEpoch );



        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////             PROPAGATE ORBIT            ////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // Create simulation object and propagate dynamics.
        SingleArcDynamicsSimulator< > dynamicsSimulator( bodyMap_, integratorSettingsC, propagatorSettings );
        std::map< double, Eigen::VectorXd > integrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );

        Eigen::Vector6d coastingState2 = ( --integrationResult.end( ) )->second;
    //    Eigen::Vector6d coastingState2 = departureState;


        // Phase 3 - powered

        try
        {
//            SphericalShaping shapeBasedMethod3 = SphericalShaping(coastingState2,arrivalState,timeOfFlight3,numberOfRevolutions_, bodyMap_,bodyToPropagate_,centralBody_,initialValueFreeCoefficient_,rootFinderSettings_,lowerBoundFreeCoefficient_,upperBoundFreeCoefficient_,integratorSettings3 );
            InvPolyShaping   shapeBasedMethod3 = InvPolyShaping(coastingState2,arrivalState,timeOfFlight3,numberOfRevolutions_,bodyMap_,bodyToPropagate_,centralBody_,initialValueFreeCoefficient_,rootFinderSettings_,lowerBoundFreeCoefficient_,upperBoundFreeCoefficient_,integratorSettings3 );


            bool infeasibleTOF = shapeBasedMethod3.getInfeasibleTOF();
            if (infeasibleTOF)
            {
                deltaV3           = 1e22;
                peakAcceleration3 = 1e22;
            }
            else
            {
                deltaV3 = shapeBasedMethod3.computeDeltaV();

                // Peak Acceleration
                // trajectory and thrust
                currentTravelledAngle = shapeBasedMethod3.getTravelledAzimuthAngle();
                stepSize = ( currentTravelledAngle ) / 100.0;
                // Check that the trajectory is feasible, ie curved toward the central body.
                peakAcceleration3 = 0.0;
                for ( int i = 0 ; i <= 100 ; i++ )
                {
                    double currentThetaAngle = i * stepSize;
                    if ( (std::abs(shapeBasedMethod3.computeCurrentThrustAccelerationMagnitude( currentThetaAngle ))  * physical_constants::ASTRONOMICAL_UNIT / std::pow(tudat::physical_constants::JULIAN_YEAR, 2.0) ) >  peakAcceleration3 )
                    {
                        // thrust is non-dimensional so dimensionalisation is needed
                        peakAcceleration3 = std::abs(shapeBasedMethod3.computeCurrentThrustAccelerationMagnitude( currentThetaAngle )) * physical_constants::ASTRONOMICAL_UNIT / std::pow(tudat::physical_constants::JULIAN_YEAR, 2.0);
                    }
                } // for loop stateVector
            }

        }
        catch (const std::runtime_error& error)
        {
            deltaV3           = 1e22;
            peakAcceleration3 = 1e22;

        }
    } // if dv spherical !=1e22
    else
    {
        deltaV3           = 1e22;
        peakAcceleration3 = 1e22;
    }


    double deltaV = deltaV1 + deltaV3;

    double peakAcc;
    if (peakAcceleration1 > peakAcceleration3)
    {
        peakAcc = peakAcceleration1;
    }
    else
    {
        peakAcc = peakAcceleration3;
    }

    //storing fitnessvalues
    std::vector< double > fitnessVector;

    fitnessVector.push_back(deltaV/1e3);

    fitnessVector.push_back(peakAcc*1e3);


//    std::cout << "----------------------------------------------" << std::endl;
//    std::cout << "1 pop" << std::endl;
//    std::cout << "----------------------------------------------" << std::endl;


    return fitnessVector;
} // IP Multi-Optimisation


std::vector< double > ShapeBasedIPSSCoastingMultiOptimisationProblem::fitness( const std::vector< double > &x ) const
{
    using namespace tudat;
    using namespace tudat::simulation_setup;
    using namespace tudat::propagators;
    using namespace tudat::numerical_integrators;
    using namespace tudat::basic_astrodynamics;
    using namespace tudat::basic_mathematics;
    using namespace tudat::orbital_element_conversions;
    using namespace tudat::unit_conversions;
    using namespace tudat::input_output;


    // inputs
    double departureTime = x.at( 0 );

    double timeOfFlight1  = x.at( 1 );
    double timeOfFlightC  = x.at( 2 );
    double timeOfFlight3  = x.at( 3 );

    double radius            = x.at( 4 );
    double azimuth           = x.at( 5 );
    double elevation         = x.at( 6 );
    double radialVelocity    = x.at( 7 );
    double azimuthalVelocity = x.at( 8 );
    double phiVelocity       = x.at( 9 );

    double timeOfFlight = timeOfFlight1 + timeOfFlightC + timeOfFlight3;

    double arrivalTime = departureTime + timeOfFlight;

    double centralBodyGravitationalParameter = bodyMap_["Sun"]->getGravityFieldModel()->getGravitationalParameter( );

    Eigen::Vector6d departureState = initialStateFunction_( departureTime );
    Eigen::Vector6d arrivalState = finalStateFunction_( arrivalTime, centralBodyGravitationalParameter );

    std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings1 =
            std::make_shared< numerical_integrators::IntegratorSettings< double > > (
                numerical_integrators::rungeKutta4, 0.0, timeOfFlight1 / 500.0 );

    std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings3 =
            std::make_shared< numerical_integrators::IntegratorSettings< double > > (
                numerical_integrators::rungeKutta4, 0.0, timeOfFlight3 / 500.0 );



//    Eigen::Vector6d initialStateSphericalCoordinates = coordinate_conversions::convertCartesianToSphericalState( departureState );
//    Eigen::Vector6d finalStateSphericalCoordinates = coordinate_conversions::convertCartesianToSphericalState( arrivalState );

//     coastingState 1
    Eigen::Vector6d sphericalState1 ;
    sphericalState1(0) = radius*physical_constants::ASTRONOMICAL_UNIT; // r
    sphericalState1(1) = azimuth;                                      // theta
    sphericalState1(2) = elevation;                                    // elevation
    sphericalState1(3) = radialVelocity;//*initialStateSphericalCoordinates(3);                               // Vr
    sphericalState1(4) = azimuthalVelocity;//*initialStateSphericalCoordinates(4);                            // Vtheta
    sphericalState1(5) = phiVelocity;//*initialStateSphericalCoordinates(5);                                  // Vphi
    Eigen::Vector6d coastingState1 = coordinate_conversions::convertSphericalToCartesianState( sphericalState1 );
//    Eigen::Vector6d coastingState1 = arrivalState;



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
                bodyMap_, accelerationMap, bodiesToPropagate, centralBodies );


    double simulationStartEpoch = 0.0;
    const double fixedStepSize = 24*3600.0;
    std::shared_ptr< IntegratorSettings< > > integratorSettingsC =
            std::make_shared< IntegratorSettings< > >
            ( rungeKutta4, simulationStartEpoch, fixedStepSize );



    double currentTravelledAngle;
    double stepSize;

    // Phase 1
    double deltaV1;
    double peakAcceleration1;
    try
    {
//        SphericalShaping shapeBasedMethod1 = SphericalShaping(departureState,coastingState1,timeOfFlight1,numberOfRevolutions_, bodyMap_,bodyToPropagate_,centralBody_,initialValueFreeCoefficient_,rootFinderSettings_,lowerBoundFreeCoefficient_,upperBoundFreeCoefficient_,integratorSettings1 );
        InvPolyShaping   shapeBasedMethod1 = InvPolyShaping(departureState,coastingState1,timeOfFlight1,numberOfRevolutions_,bodyMap_,bodyToPropagate_,centralBody_,initialValueFreeCoefficient_,rootFinderSettings_,lowerBoundFreeCoefficient_,upperBoundFreeCoefficient_,integratorSettings1);

        // shaping
        bool infeasibleTOF = shapeBasedMethod1.getInfeasibleTOF();
        if (infeasibleTOF)
        {
            deltaV1           = 1e22;
            peakAcceleration1 = 1e22;
        }
        else
        {
            deltaV1 = shapeBasedMethod1.computeDeltaV();

            // Peak Acceleration
            // trajectory and thrust
            currentTravelledAngle = shapeBasedMethod1.getTravelledAzimuthAngle();
            stepSize = ( currentTravelledAngle ) / 100.0;
            // Check that the trajectory is feasible, ie curved toward the central body.
            peakAcceleration1 = 0.0;
            for ( int i = 0 ; i <= 100 ; i++ )
            {
                double currentThetaAngle = i * stepSize;
                if ( (std::abs(shapeBasedMethod1.computeCurrentThrustAccelerationMagnitude( currentThetaAngle ))  * physical_constants::ASTRONOMICAL_UNIT / std::pow(tudat::physical_constants::JULIAN_YEAR, 2.0) ) >  peakAcceleration1 )
                {
                    // thrust is non-dimensional so dimensionalisation is needed
                    peakAcceleration1 = std::abs(shapeBasedMethod1.computeCurrentThrustAccelerationMagnitude( currentThetaAngle )) * physical_constants::ASTRONOMICAL_UNIT / std::pow(tudat::physical_constants::JULIAN_YEAR, 2.0);
                }
            } // for loop stateVector
        }

    }
    catch (const std::runtime_error& error)
    {
        deltaV1 = 1e22;
    }


    double deltaV3;
    double peakAcceleration3;

    if (deltaV1!=1e22)
    {
        // Phase 2 - Coasting
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////             CREATE PROPAGATION SETTINGS            ////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // Set initial conditions for the Asterix satellite that will be propagated in this simulation.
        // The initial conditions are given in Keplerian elements and later on converted to Cartesian
        // elements.

        // Set simulation end epoch.
        const double simulationEndEpoch = timeOfFlightC;

        // Create propagator settings.
        std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
                std::make_shared< TranslationalStatePropagatorSettings< double > >
                ( centralBodies, accelerationModelMap, bodiesToPropagate, coastingState1, simulationEndEpoch );



        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////             PROPAGATE ORBIT            ////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // Create simulation object and propagate dynamics.
        SingleArcDynamicsSimulator< > dynamicsSimulator( bodyMap_, integratorSettingsC, propagatorSettings );
        std::map< double, Eigen::VectorXd > integrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );

        Eigen::Vector6d coastingState2 = ( --integrationResult.end( ) )->second;
    //    Eigen::Vector6d coastingState2 = departureState;


        // Phase 3 - powered

        try
        {
            SphericalShaping shapeBasedMethod3 = SphericalShaping(coastingState2,arrivalState,timeOfFlight3,numberOfRevolutions_, bodyMap_,bodyToPropagate_,centralBody_,initialValueFreeCoefficient_,rootFinderSettings_,lowerBoundFreeCoefficient_,upperBoundFreeCoefficient_,integratorSettings3 );
//            InvPolyShaping   shapeBasedMethod3 = InvPolyShaping(coastingState2,arrivalState,timeOfFlight3,numberOfRevolutions_,bodyMap_,bodyToPropagate_,centralBody_,initialValueFreeCoefficient_,rootFinderSettings_,lowerBoundFreeCoefficient_,upperBoundFreeCoefficient_,integratorSettings3 );


            bool infeasibleTOF = shapeBasedMethod3.getInfeasibleTOF();
            if (infeasibleTOF)
            {
                deltaV3           = 1e22;
                peakAcceleration3 = 1e22;
            }
            else
            {
                deltaV3 = shapeBasedMethod3.computeDeltaV();

                // Peak Acceleration
                // trajectory and thrust
                currentTravelledAngle = shapeBasedMethod3.getFinalAzimuthAngle() - shapeBasedMethod3.getInitialAzimuthAngle();
                stepSize = ( currentTravelledAngle ) / 100.0;
                // Check that the trajectory is feasible, ie curved toward the central body.
                peakAcceleration3 = 0.0;
                for ( int i = 0 ; i <= 100 ; i++ )
                {
//                    double currentThetaAngle = i * stepSize;
                    double currentThetaAngle = shapeBasedMethod3.getInitialAzimuthAngle() + i * stepSize;

                    if ( (shapeBasedMethod3.computeCurrentThrustAccelerationVector( currentThetaAngle ).norm() ) >  peakAcceleration3 )
                    {
                        // thrust is non-dimensional so dimensionalisation is needed
                        peakAcceleration3 = shapeBasedMethod3.computeCurrentThrustAccelerationVector( currentThetaAngle ).norm();
                    }
                } // for loop stateVector
            }

        }
        catch (const std::runtime_error& error)
        {
            deltaV3           = 1e22;
            peakAcceleration3 = 1e22;

        }
    } // if dv spherical !=1e22
    else
    {
        deltaV3           = 1e22;
        peakAcceleration3 = 1e22;
    }


    double deltaV = deltaV1 + deltaV3;

    double peakAcc;
    if (peakAcceleration1 > peakAcceleration3)
    {
        peakAcc = peakAcceleration1;
    }
    else
    {
        peakAcc = peakAcceleration3;
    }

    //storing fitnessvalues
    std::vector< double > fitnessVector;

    fitnessVector.push_back(deltaV);

    fitnessVector.push_back(peakAcc);


//    std::cout << "----------------------------------------------" << std::endl;
//    std::cout << "1 pop" << std::endl;
//    std::cout << "----------------------------------------------" << std::endl;


    return fitnessVector;
} // IP-SS Multi-Optimisation


} // namespace shape_based_methods
} // namespace tudat



