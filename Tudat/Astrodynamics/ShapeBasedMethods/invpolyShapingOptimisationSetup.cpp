/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "invpolyShapingOptimisationSetup.h"

#include "Tudat/Astrodynamics/ShapeBasedMethods/compositeFunctionHodographicShaping.h"
#include "Tudat/Astrodynamics/ShapeBasedMethods/hodographicShaping.h"
#include "Tudat/Astrodynamics/ShapeBasedMethods/baseFunctionsHodographicShaping.h"
#include "Tudat/Astrodynamics/ShapeBasedMethods/createBaseFunctionHodographicShaping.h"



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
std::vector< double > InvPolyShapingOptimisationProblem::fitness( const std::vector< double > &x ) const
{

//    std::cout << x.at(0)  << " " << x.at(1) << std::endl;

    double windingParameter = 0.125;

    double departureTime = x.at( 0 );
    double timeOfFlight = x.at( 1 );
    double arrivalTime = departureTime + timeOfFlight;

    double centralBodyGravitationalParameter = bodyMap_["Sun"]->getGravityFieldModel()->getGravitationalParameter( );

    Eigen::Vector6d departureState = initialStateFunction_( departureTime );
    Eigen::Vector6d arrivalState = finalStateFunction_( arrivalTime, centralBodyGravitationalParameter );

    std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings =
            std::make_shared< numerical_integrators::IntegratorSettings< double > > (
                numerical_integrators::rungeKutta4, 0.0, timeOfFlight / 500.0 );


    // /////////////////////////// HODOGRAPHIC - START ////////////////////////////
    // Initialize free coefficients vector for radial velocity function.
    Eigen::VectorXd freeCoefficientsRadialVelocityFunction = Eigen::VectorXd::Zero( 0 );

    // Initialize free coefficients vector for normal velocity function.
    Eigen::VectorXd freeCoefficientsNormalVelocityFunction = Eigen::VectorXd::Zero( 0 );

    // Initialize free coefficients vector for axial velocity function.
    Eigen::VectorXd freeCoefficientsAxialVelocityFunction = Eigen::VectorXd::Zero( 0 );

    double frequency = 2.0 * mathematical_constants::PI / ( timeOfFlight );

    double scaleFactor = 1.0 / ( timeOfFlight );

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
    std::vector< std::shared_ptr< BaseFunctionHodographicShaping > > normalVelocityFunctionComponents;
    normalVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( constant, firstNormalVelocityBaseFunctionSettings ) );
    normalVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPower, secondNormalVelocityBaseFunctionSettings ) );
    normalVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( cosine, thirdNormalVelocityBaseFunctionSettings ) );


    // Create base function settings for the components of the axial velocity composite function.
    std::shared_ptr< BaseFunctionHodographicShapingSettings > firstAxialVelocityBaseFunctionSettings =
            std::make_shared< TrigonometricFunctionHodographicShapingSettings >( ( numberOfRevolutions_ + 0.5 ) * frequency );
    std::shared_ptr< BaseFunctionHodographicShapingSettings > secondAxialVelocityBaseFunctionSettings =
            std::make_shared< PowerTimesTrigonometricFunctionHodographicShapingSettings >
            ( 3.0, ( numberOfRevolutions_ + 0.5 ) * frequency, scaleFactor );
    std::shared_ptr< BaseFunctionHodographicShapingSettings > thirdAxialVelocityBaseFunctionSettings =
            std::make_shared< PowerTimesTrigonometricFunctionHodographicShapingSettings >(
                3.0, ( numberOfRevolutions_ + 0.5 ) * frequency, scaleFactor );

    // Set components for the axial velocity function.
    std::vector< std::shared_ptr< BaseFunctionHodographicShaping > > axialVelocityFunctionComponents;
    axialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( cosine, firstAxialVelocityBaseFunctionSettings ) );
    axialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPowerCosine, secondAxialVelocityBaseFunctionSettings ) );
    axialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPowerSine, thirdAxialVelocityBaseFunctionSettings ) );



    // /////////////////////////// HODOGRAPHIC ////////////////////////////



    // creating shaping objects
    InvPolyShaping invPolyShaping     = InvPolyShaping(departureState,arrivalState,timeOfFlight,numberOfRevolutions_,bodyMap_,bodyToPropagate_,centralBody_,initialValueFreeCoefficient_,rootFinderSettings_,lowerBoundFreeCoefficient_,upperBoundFreeCoefficient_,integratorSettings);

    ExposinsShaping exposinsShaping   = ExposinsShaping(departureState,arrivalState,timeOfFlight,numberOfRevolutions_,bodyMap_,bodyToPropagate_,centralBody_,windingParameter,rootFinderSettings_,integratorSettings);

    HodographicShaping hodographicShaping = HodographicShaping(departureState,arrivalState,timeOfFlight,numberOfRevolutions_,bodyMap_,bodyToPropagate_,centralBody_,radialVelocityFunctionComponents,normalVelocityFunctionComponents,axialVelocityFunctionComponents,freeCoefficientsRadialVelocityFunction,freeCoefficientsNormalVelocityFunction,freeCoefficientsAxialVelocityFunction,integratorSettings );

    double deltaVspherical;
    try
    {        
        SphericalShaping sphericalShaping = SphericalShaping(departureState,arrivalState,timeOfFlight,numberOfRevolutions_, bodyMap_,bodyToPropagate_,centralBody_,initialValueFreeCoefficient_,rootFinderSettings_,lowerBoundFreeCoefficient_,upperBoundFreeCoefficient_,integratorSettings );
        
        deltaVspherical = sphericalShaping.computeDeltaV();
    }
    catch (const std::runtime_error& error)
    {
        deltaVspherical = 1e22;
    }
    

    // computing and comparing fitness values
    double deltaVexposins;
    bool infeasibleTOF = exposinsShaping.getInfeasibleTOF();
    if (infeasibleTOF)
    {
        deltaVexposins = 1e20;
    }
    else
    {
        deltaVexposins = exposinsShaping.computeDeltaV() + exposinsShaping.computeDeltaVBoundaries();
    }

    double deltaVinvpoly;
    infeasibleTOF = invPolyShaping.getInfeasibleTOF();
    if (infeasibleTOF)
    {
        deltaVinvpoly = 1e21;
    }
    else
    {
        deltaVinvpoly = invPolyShaping.computeDeltaV( );
    }

    double deltaVhodographic = hodographicShaping.computeDeltaV();

    //storing fitnessvalues
    std::vector< double > fitnessVector;
    double deltaVbest = deltaVexposins;

    if ( deltaVinvpoly <= deltaVbest )
    {
        deltaVbest = deltaVinvpoly;
    }
    if ( deltaVspherical <= deltaVbest )
    {
        deltaVbest = deltaVspherical;
    }
    if ( deltaVhodographic <= deltaVbest )
    {
        deltaVbest = deltaVhodographic;
    }

    fitnessVector.push_back(deltaVbest);

//    std::cout << x.at(0)  << " " << x.at(1) << " " << invPolyShaping.computeDeltaV( ) << std::endl;


    return fitnessVector;
} // Full Optimisation

// Calculates the fitness
std::vector< double > ExposinsShapingOptimisationProblem::fitness( const std::vector< double > &x ) const
{
    double departureTime    = x.at( 0 );
    double timeOfFlight     = x.at( 1 );
    double windingParameter = x.at( 2 );

    double arrivalTime = departureTime + timeOfFlight;

    double centralBodyGravitationalParameter = bodyMap_["Sun"]->getGravityFieldModel()->getGravitationalParameter( );

    Eigen::Vector6d departureState = initialStateFunction_( departureTime );
    Eigen::Vector6d arrivalState = finalStateFunction_( arrivalTime, centralBodyGravitationalParameter );

    std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings =
            std::make_shared< numerical_integrators::IntegratorSettings< double > > (
                numerical_integrators::rungeKutta4, 0.0, timeOfFlight / 500.0 );


    // creating shaping objects
    ExposinsShaping exposinsShaping   = ExposinsShaping(departureState,arrivalState,timeOfFlight,numberOfRevolutions_,bodyMap_,bodyToPropagate_,centralBody_,windingParameter,rootFinderSettings_,integratorSettings);


    // computing and comparing fitness values
    double deltaVexposins;
    bool infeasibleTOF = exposinsShaping.getInfeasibleTOF();
    if (infeasibleTOF)
    {
        deltaVexposins = 1e20;
    }
    else
    {
        deltaVexposins = exposinsShaping.computeDeltaV() + exposinsShaping.computeDeltaVBoundaries();
    }

    //storing fitnessvalues
    std::vector< double > fitnessVector;

    fitnessVector.push_back(deltaVexposins);



    return fitnessVector;
} // Exposins Optimisation

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
} // Invpoly Optimisation

// Calculates the fitness
std::vector< double > HodographicShapingOptimisationProblem::fitness( const std::vector< double > &x ) const
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


    // /////////////////////////// HODOGRAPHIC - START ////////////////////////////
    // Initialize free coefficients vector for radial velocity function.
    Eigen::VectorXd freeCoefficientsRadialVelocityFunction = Eigen::VectorXd::Zero( 0 );

    // Initialize free coefficients vector for normal velocity function.
    Eigen::VectorXd freeCoefficientsNormalVelocityFunction = Eigen::VectorXd::Zero( 0 );

    // Initialize free coefficients vector for axial velocity function.
    Eigen::VectorXd freeCoefficientsAxialVelocityFunction = Eigen::VectorXd::Zero( 0 );

    double frequency = 2.0 * mathematical_constants::PI / ( timeOfFlight );

    double scaleFactor = 1.0 / ( timeOfFlight );

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
    std::vector< std::shared_ptr< BaseFunctionHodographicShaping > > normalVelocityFunctionComponents;
    normalVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( constant, firstNormalVelocityBaseFunctionSettings ) );
    normalVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPower, secondNormalVelocityBaseFunctionSettings ) );
    normalVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( cosine, thirdNormalVelocityBaseFunctionSettings ) );


    // Create base function settings for the components of the axial velocity composite function.
    std::shared_ptr< BaseFunctionHodographicShapingSettings > firstAxialVelocityBaseFunctionSettings =
            std::make_shared< TrigonometricFunctionHodographicShapingSettings >( ( numberOfRevolutions_ + 0.5 ) * frequency );
    std::shared_ptr< BaseFunctionHodographicShapingSettings > secondAxialVelocityBaseFunctionSettings =
            std::make_shared< PowerTimesTrigonometricFunctionHodographicShapingSettings >
            ( 3.0, ( numberOfRevolutions_ + 0.5 ) * frequency, scaleFactor );
    std::shared_ptr< BaseFunctionHodographicShapingSettings > thirdAxialVelocityBaseFunctionSettings =
            std::make_shared< PowerTimesTrigonometricFunctionHodographicShapingSettings >(
                3.0, ( numberOfRevolutions_ + 0.5 ) * frequency, scaleFactor );

    // Set components for the axial velocity function.
    std::vector< std::shared_ptr< BaseFunctionHodographicShaping > > axialVelocityFunctionComponents;
    axialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( cosine, firstAxialVelocityBaseFunctionSettings ) );
    axialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPowerCosine, secondAxialVelocityBaseFunctionSettings ) );
    axialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPowerSine, thirdAxialVelocityBaseFunctionSettings ) );



    // /////////////////////////// HODOGRAPHIC ////////////////////////////

    //storing fitnessvalues
    std::vector< double > fitnessVector;

    HodographicShaping hodographicShaping = HodographicShaping(departureState,arrivalState,timeOfFlight,numberOfRevolutions_,bodyMap_,bodyToPropagate_,centralBody_,radialVelocityFunctionComponents,normalVelocityFunctionComponents,axialVelocityFunctionComponents,freeCoefficientsRadialVelocityFunction,freeCoefficientsNormalVelocityFunction,freeCoefficientsAxialVelocityFunction,integratorSettings );

    double deltaVhodographic = hodographicShaping.computeDeltaV();

    fitnessVector.push_back(deltaVhodographic);



    return fitnessVector;
} // Invpoly Optimisation



} // namespace shape_based_methods
} // namespace tudat



