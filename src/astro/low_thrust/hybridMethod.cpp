/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */


#include "tudat/astro/LowThrustTrajectories/hybridMethod.h"
#include "tudat/astro/LowThrustTrajectories/hybridOptimisationSetup.h"
#include "pagmo/problems/unconstrain.hpp"
#include "pagmo/algorithms/compass_search.hpp"
#include "tudat/astro/propulsion/thrustMagnitudeWrapper.h"
#include "tudat/astro/propulsion/costateBasedThrustGuidance.h"


namespace tudat
{
namespace low_thrust_trajectories
{


//! Perform optimisation.
std::pair< std::vector< double >, std::vector< double > > HybridMethod::performOptimisation( )
{
    //Set seed for reproducible results
    pagmo::random_device::set_seed( 456 );

    // Create object to compute the problem fitness
    problem prob{ HybridMethodProblem( stateAtDeparture_, stateAtArrival_, maximumThrust_, specificImpulse_, timeOfFlight_,
                                       initialGuessThrustModel_, initialAndFinalMEEcostatesBounds_,
                                       optimisationSettings_->relativeToleranceConstraints_ )};

    std::vector< double > constraintsTolerance;
    for ( unsigned int i = 0 ; i < ( prob.get_nec() + prob.get_nic() ) ; i++ )
    {
        constraintsTolerance.push_back( 1.0e-3 );
    }
    prob.set_c_tol( constraintsTolerance );

    algorithm algo = optimisationSettings_->optimisationAlgorithm_;

    unsigned long long populationSize = optimisationSettings_->numberOfIndividualsPerPopulation_;

    island island{ algo, prob, populationSize };

    // Evolve for 10 generations
    for( int i = 0 ; i < optimisationSettings_->numberOfGenerations_ ; i++ )
    {
        island.evolve( );
        while( island.status( ) != pagmo::evolve_status::idle &&
               island.status( ) != pagmo::evolve_status::idle_error )
        {
            island.wait( );
        }
        island.wait_check( ); // Raises errors

    }

    championFitness_ = island.get_population().champion_f();
    championDesignVariables_ = island.get_population().champion_x();

    std::pair< std::vector< double >, std::vector< double > > output;
    output.first = championFitness_;
    output.second = championDesignVariables_;

    return output;

}


//! Compute current thrust vector.
Eigen::Vector3d HybridMethod::computeCurrentThrustForce(
        double time,
        std::function< double ( const double ) > specificImpulseFunction,
        std::shared_ptr<numerical_integrators::IntegratorSettings< double > > integratorSettings )
{

    std::shared_ptr< simulation_setup::ThrustMagnitudeSettings > thrustMagnitudeSettings = hybridMethodModel_->getMEEcostatesBasedThrustMagnitudeSettings( );

    std::function< Eigen::Vector3d( ) > bodyFixedThrustDirection = simulation_setup::getBodyFixedThrustDirection(
                thrustMagnitudeSettings, bodies_, bodyToPropagate_ );

    std::function< Eigen::Vector6d( ) > thrustingBodyStateFunction = [ = ] ( )
    {
        return computeCurrentStateVector( time );
    };

    std::function< Eigen::Vector6d( ) > centralBodyStateFunction = [ = ] ( )
    {
        return Eigen::Vector6d::Zero( );
    };

    std::function< double( ) > centralBodyGravitationalParameterFunction = [ = ]( )
    {
        return bodies_[ centralBody_ ]->getGravityFieldModel( )->getGravitationalParameter( );
    };

    std::function< double( ) > thrustingBodyMassFunction = std::bind( &simulation_setup::Body::getBodyMass, bodies_.at( bodyToPropagate_ ) );


    propulsion::MeeCostatesBangBangThrustMagnitudeWrapper thrustMagnitudeWrapper = propulsion::MeeCostatesBangBangThrustMagnitudeWrapper(
                thrustingBodyStateFunction, centralBodyStateFunction, centralBodyGravitationalParameterFunction,
                hybridMethodModel_->getCostatesFunction_( ), maximumThrust_, specificImpulseFunction, thrustingBodyMassFunction);

    propulsion::MeeCostateBasedThrustGuidance thrustGuidance = propulsion::MeeCostateBasedThrustGuidance(
                thrustingBodyStateFunction, centralBodyStateFunction, centralBodyGravitationalParameterFunction,
                hybridMethodModel_->getCostatesFunction_( ), bodyFixedThrustDirection );

    thrustGuidance.updateCalculator( time );
    thrustMagnitudeWrapper.update( time );

    return thrustMagnitudeWrapper.getCurrentThrustMagnitude( ) * thrustGuidance.getCurrentForceDirectionInPropagationFrame( );
}



//! Return thrust profile.
void HybridMethod::getThrustForceProfile(
        std::vector< double >& epochsVector,
        std::map< double, Eigen::VectorXd >& thrustProfile,
        std::function< double ( const double ) > specificImpulseFunction )
{
    thrustProfile.clear( );

    std::shared_ptr< simulation_setup::ThrustMagnitudeSettings > thrustMagnitudeSettings = hybridMethodModel_->getMEEcostatesBasedThrustMagnitudeSettings( );

    std::function< Eigen::Vector3d( ) > bodyFixedThrustDirection = simulation_setup::getBodyFixedThrustDirection(
                thrustMagnitudeSettings, bodies_, bodyToPropagate_ );

    std::map< double, Eigen::Vector6d > trajectory;
    getTrajectory( epochsVector, trajectory );

    for ( unsigned int i = 0 ; i < epochsVector.size() ; i++ )
    {
        if ( ( i > 0 ) && ( epochsVector[ i ] < epochsVector[ i - 1 ] ) )
        {
            throw std::runtime_error( "Error when retrieving the thrust profile of a hybrid method trajectory, "
                                      "epochs are not provided in increasing order." );
        }

        Eigen::Vector6d currentStateVector = trajectory[ epochsVector[ i ] ];

        std::function< Eigen::Vector6d( ) > thrustingBodyStateFunction = [ = ] ( )
        {
            return currentStateVector;
        };

        std::function< Eigen::Vector6d( ) > centralBodyStateFunction = [ = ] ( )
        {
            return Eigen::Vector6d::Zero( );
        };

        std::function< double( ) > centralBodyGravitationalParameterFunction = [ = ]( )
        {
            return bodies_[ centralBody_ ]->getGravityFieldModel( )->getGravitationalParameter( );
        };

        std::function< double( ) > thrustingBodyMassFunction = std::bind( &simulation_setup::Body::getBodyMass, bodies_.at( bodyToPropagate_ ) );


        propulsion::MeeCostatesBangBangThrustMagnitudeWrapper thrustMagnitudeWrapper = propulsion::MeeCostatesBangBangThrustMagnitudeWrapper(
                    thrustingBodyStateFunction, centralBodyStateFunction, centralBodyGravitationalParameterFunction,
                    hybridMethodModel_->getCostatesFunction_( ), maximumThrust_, specificImpulseFunction, thrustingBodyMassFunction);

        propulsion::MeeCostateBasedThrustGuidance thrustGuidance = propulsion::MeeCostateBasedThrustGuidance(
                    thrustingBodyStateFunction, centralBodyStateFunction, centralBodyGravitationalParameterFunction,
                    hybridMethodModel_->getCostatesFunction_( ), bodyFixedThrustDirection );

        thrustGuidance.updateCalculator( epochsVector[ i ] );
        thrustMagnitudeWrapper.update( epochsVector[ i ] );

        thrustProfile[ epochsVector[ i ] ] = thrustMagnitudeWrapper.getCurrentThrustMagnitude( ) * thrustGuidance.getCurrentForceDirectionInPropagationFrame( );

    }
}

//! Compute magnitude thrust acceleration.
double HybridMethod::computeCurrentThrustAccelerationMagnitude(
        double currentTime, std::function< double ( const double ) > specificImpulseFunction,
        std::shared_ptr<numerical_integrators::IntegratorSettings< double > > integratorSettings )
{

    double currentMass = computeCurrentMass( currentTime, specificImpulseFunction, integratorSettings );
    Eigen::Vector3d currentThrustVector = computeCurrentThrustForce( currentTime, specificImpulseFunction, integratorSettings_ );

    return currentThrustVector.norm( ) / currentMass;

}


//! Compute direction thrust acceleration in cartesian coordinates.
Eigen::Vector3d HybridMethod::computeCurrentThrustAccelerationDirection(
        double currentTime, std::function< double ( const double ) > specificImpulseFunction,
        std::shared_ptr<numerical_integrators::IntegratorSettings< double > > integratorSettings )
{

    Eigen::Vector3d currentThrustVector = computeCurrentThrustForce( currentTime, specificImpulseFunction, integratorSettings );

    Eigen::Vector3d thrustAcceleration = currentThrustVector.normalized( );

    return thrustAcceleration.normalized( );
}



//! Return thrust acceleration profile.
void HybridMethod::getThrustAccelerationProfile(
        std::vector< double >& epochsVector,
        std::map< double, Eigen::VectorXd >& thrustAccelerationProfile,
        std::function< double ( const double ) > specificImpulseFunction,
        std::shared_ptr<numerical_integrators::IntegratorSettings< double > > integratorSettings )
{
    thrustAccelerationProfile.clear();

    std::map< double, Eigen::VectorXd > thrustProfile;
    getThrustForceProfile( epochsVector, thrustProfile, specificImpulseFunction, integratorSettings );

    std::map< double, Eigen::VectorXd > massProfile;
    getMassProfile( epochsVector, massProfile, specificImpulseFunction, integratorSettings );

    for ( unsigned int i = 0 ; i < epochsVector.size() ; i++ )
    {
        if ( ( i > 0 ) && ( epochsVector[ i ] < epochsVector[ i - 1 ] ) )
        {
            throw std::runtime_error( "Error when retrieving the thrust profile of a hybrid method trajectory, "
                                      "epochs are not provided in increasing order." );
        }

        Eigen::Vector3d currentThrustVector = thrustProfile[ epochsVector[ i ] ];
        double currentMass = massProfile[ epochsVector[ i ] ][ 0 ];

        Eigen::Vector3d currentThrustAccelerationVector = currentThrustVector / currentMass;

        thrustAccelerationProfile[ epochsVector[ i ] ] = currentThrustAccelerationVector;

    }
}


//! Compute current cartesian state.
Eigen::Vector6d HybridMethod::computeCurrentStateVector( const double currentTime )
{
    Eigen::Vector6d stateVector;
    if ( currentTime == 0.0 )
    {
        stateVector = stateAtDeparture_;
    }
    else
    {
        stateVector = hybridMethodModel_->propagateTrajectory( 0.0, currentTime, stateAtDeparture_, initialSpacecraftMass_ );
    }

    return stateVector;

}


//! Retrieve acceleration map (thrust and central gravity accelerations).
basic_astrodynamics::AccelerationMap HybridMethod::retrieveLowThrustAccelerationMap(
        const simulation_setup::SystemOfBodies& bodies,
        const std::string& bodyToPropagate,
        const std::string& centralBody,
        const std::function< double ( const double ) > specificImpulseFunction,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings )
{
    basic_astrodynamics::AccelerationMap hybridMethodAccelerationMap = hybridMethodModel_->getLowThrustTrajectoryAccelerationMap( );
    return hybridMethodAccelerationMap;
}


//! Define appropriate translational state propagator settings for the full propagation.
std::pair< std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > >,
std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > > > HybridMethod::createLowThrustTranslationalStatePropagatorSettings(
        basic_astrodynamics::AccelerationMap accelerationModelMap,
        std::shared_ptr< propagators::DependentVariableSaveSettings > dependentVariablesToSave )
{

    // Create termination conditions settings.
    std::pair< std::shared_ptr< propagators::PropagationTerminationSettings >,
            std::shared_ptr< propagators::PropagationTerminationSettings > > terminationConditions;

    terminationConditions.first = std::make_shared< propagators::PropagationTimeTerminationSettings >( 0.0, true );
    terminationConditions.second = std::make_shared< propagators::PropagationTimeTerminationSettings >( timeOfFlight_, true );

    // Compute state vector at half of the time of flight.
    Eigen::Vector6d stateAtHalfOfTimeOfFlight = computeCurrentStateVector( timeOfFlight_ / 2.0 );

    // Define translational state propagator settings.
    std::pair< std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > >,
            std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > > > translationalStatePropagatorSettings;

    // Define backward translational state propagation settings.
    translationalStatePropagatorSettings.first = std::make_shared< propagators::TranslationalStatePropagatorSettings< double > >
            ( std::vector< std::string >{ centralBody_ }, accelerationModelMap,
              std::vector< std::string >{ bodyToPropagate_ }, stateAtHalfOfTimeOfFlight,
              terminationConditions.first, propagators::gauss_modified_equinoctial, dependentVariablesToSave );

    // Define forward translational state propagation settings.
    translationalStatePropagatorSettings.second = std::make_shared< propagators::TranslationalStatePropagatorSettings< double > >
            ( std::vector< std::string >{ centralBody_ }, accelerationModelMap,
              std::vector< std::string >{ bodyToPropagate_ }, stateAtHalfOfTimeOfFlight,
              terminationConditions.second, propagators::gauss_modified_equinoctial, dependentVariablesToSave );

    return translationalStatePropagatorSettings;
}


} // namespace low_thrust_trajectories
} // namespace tudat
