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

#ifndef HYBRIDMETHOD_H
#define HYBRIDMETHOD_H

#include <Tudat/SimulationSetup/tudatSimulationHeader.h>
#include <math.h>
#include <vector>
#include <Eigen/Dense>
#include <map>
#include "pagmo/algorithm.hpp"

namespace tudat
{
namespace low_thrust_direct_methods
{


//! Transform thrust model as a function of time into hybrid method thrust model.
Eigen::Matrix< double, 10 , 1 > convertToHybridMethodThrustModel( std::function< Eigen::Vector3d( const double ) > thrustModelWrtTime );

class HybridMethod
{
public:

    //! Constructor.
    HybridMethod(
            const Eigen::Vector6d& stateAtDeparture,
            const Eigen::Vector6d& stateAtArrival,
            const double maximumThrust,
            const double specificImpulse,
            const double timeOfFlight,
            simulation_setup::NamedBodyMap bodyMap,
            const std::string bodyToPropagate,
            const std::string centralBody,
            std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
            pagmo::algorithm optimisationAlgorithm,
            const int numberOfGenerations,
            const int numberOfIndividualsPerPopulation,
            const double relativeToleranceConstraints = 1.0e-6,
            std::pair< std::function< Eigen::Vector3d( const double ) >, double > initialGuessThrustModel = std::make_pair( nullptr, 0.0 ) ) :
        stateAtDeparture_( stateAtDeparture ),
        stateAtArrival_( stateAtArrival ),
        maximumThrust_( maximumThrust ),
        specificImpulse_( specificImpulse ),
        timeOfFlight_( timeOfFlight ),
        bodyMap_( bodyMap ),
        bodyToPropagate_( bodyToPropagate ),
        centralBody_( centralBody ),
        integratorSettings_( integratorSettings ),
        optimisationAlgorithm_( optimisationAlgorithm ),
        numberOfGenerations_( numberOfGenerations ),
        numberOfIndividualsPerPopulation_( numberOfIndividualsPerPopulation ),
        relativeToleranceConstraints_( relativeToleranceConstraints ) //,
//        initialGuessThrustModel_( initialGuessThrustModel )
    {

        // Store initial spacecraft mass.
        initialSpacecraftMass_ = bodyMap_[ bodyToPropagate_ ]->getBodyMass();

        // Convert the thrust model proposed as initial guess into simplified thrust model adapted to the hybrid method.
        if ( initialGuessThrustModel.first != nullptr )
        {
            initialGuessThrustModel_.first = convertToHybridMethodThrustModel( initialGuessThrustModel.first );
        }
        else
        {
            Eigen::VectorXd emptyVector;
            initialGuessThrustModel_.first = emptyVector; //Eigen::VectorXd::Zero( 10 ); // emptyVector;
        }
        initialGuessThrustModel_.second = initialGuessThrustModel.second;

        // Perform optimisation
        std::pair< std::vector< double >, std::vector< double > > bestIndividual = performOptimisation( );
        championFitness_ = bestIndividual.first;
        championDesignVariables_ = bestIndividual.second;

    }

    //! Default destructor.
    ~HybridMethod( ) { }

    //! Perform optimisation.
    std::pair< std::vector< double >, std::vector< double > > performOptimisation( );

    //! Compute DeltaV.
    double computeDeltaV( )
    {
        return championFitness_[ 0 ];
    }

    //! Return best individual.
    std::vector< double > getBestIndividual( )
    {
        return championDesignVariables_;
    }

    //! Return fitness of best individual.
    std::vector< double > getBestIndividualFitness( )
    {
        return championFitness_;
    }


    void computeHybridMethodTrajectoryAndFullPropagation(
         std::pair< std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > >,
            std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > > >& propagatorSettings,
         std::map< double, Eigen::VectorXd >& fullPropagationResults,
         std::map< double, Eigen::Vector6d >& hybridMethodResults,
         std::map< double, Eigen::VectorXd>& dependentVariablesHistory );


protected:

private:

    //! State vector of the vehicle at the leg departure.
    Eigen::Vector6d stateAtDeparture_;

    //! State vector of the vehicle at the leg arrival.
    Eigen::Vector6d stateAtArrival_;

    //! Maximum allowed thrust.
    double maximumThrust_;

    //! Specific impulse.
    double specificImpulse_;

    //! Time of flight for the leg.
    double timeOfFlight_;

    //! Body map object.
    simulation_setup::NamedBodyMap bodyMap_;

    //! Name of the body to be propagated.
    std::string bodyToPropagate_;

    //! Name of the central body.
    std::string centralBody_;

    //! Integrator settings.
    std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings_;

    //! Optimisation algorithm to be used to solve the Sims-Flanagan problem.
    pagmo::algorithm optimisationAlgorithm_;

    //! Number of generations for the optimisation algorithm.
    int numberOfGenerations_;

    //! Number of individuals per population for the optimisation algorithm.
    int numberOfIndividualsPerPopulation_;

    //! Relative tolerance for optimisation constraints.
    double relativeToleranceConstraints_;

    //! Initial guess for the optimisation.
    //! The first element contains the thrust throttles corresponding to the initial guess for the thrust model.
    //! The second element defines the bounds around the initial time (in percentage).
    std::pair< Eigen::VectorXd, double > initialGuessThrustModel_;

    //! Fitness vector of the optimisation best individual.
    std::vector< double > championFitness_;

    //! Design variables vector corresponding to the optimisation best individual.
    std::vector< double > championDesignVariables_;

    //! Initial mass of the spacecraft.
    double initialSpacecraftMass_;

};


} // namespace low_thrust_direct_methods
} // namespace tudat

#endif // HYBRIDMETHOD_H