/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_PROPAGATION_CR3BP_FULL
#define TUDAT_PROPAGATION_CR3BP_FULL

#include "tudat/simulation/simulation.h"


namespace tudat
{

namespace propagators
{

simulation_setup::BodyListSettings setupBodySettingsCR3BP(
        const double distancePrimarySecondary,
        const std::string& namePrimaryBody,
        const std::string& nameSecondaryBody,
        const std::string& frameOrientation = "ECLIPJ2000",
        const double primaryGravitationalParameter = TUDAT_NAN,
        const double secondaryGravitationalParameter = TUDAT_NAN );

//! Setup CR3BP system of bodies.
/*!
 * Setup CR3BP system of bodies. The two primaries, as well as the third, smaller body to be propagated are defined in the system of bodies.
 * The two primaries are in circular orbit about their barycenter, orbiting it with the same mean motion, so that they stay
 * aligned during propagation.
 * \param distancePrimarySecondary Distance between primaries    [m].
 * \param namePrimaryBody Name of the primary body.
 * \param nameSecondaryBody Name of the secondary body.
 * \param nameBodyToPropagate Name of the third, smaller body to be propagated.
 * \param frameOrientation Orientation of frame in which to propagate
 * \param primaryGravitationalParameter Gravitational parameter of primary
 * \param secondaryGravitationalParameter Gravitational parameter of secondary
 * \return Body Map modelling the CR3BP.
 */
simulation_setup::SystemOfBodies setupBodyMapCR3BP(
        const double distancePrimarySecondary,
        const std::string& namePrimaryBody,
        const std::string& nameSecondaryBody,
        const std::string& nameBodyToPropagate,
        const std::string& frameOrientation = "ECLIPJ2000",
        const double primaryGravitationalParameter = TUDAT_NAN,
        const double secondaryGravitationalParameter = TUDAT_NAN );

//! Setup CR3BP acceleration map.
/*!
 * Setup CR3BP acceleration map. Define the acceleration map for the CR3BP, from the corresponding system of bodies.
 * The only accelerations acting on the system are point-mass gravity from the two primaries.
 * \param namePrimaryBody Name of the primary body.
 * \param nameSecondaryBody Name of the secondary body.
 * \param nameBodyToPropagate Name of the third, smaller body to be propagated.
 * \param centralBody Central bodys for the propagation.
 * \param bodies CR3BP system of bodies.
 * \return Acceleration map for the CR3BP.
 */
basic_astrodynamics::AccelerationMap setupAccelerationMapCR3BP(
        const std::string& namePrimaryBody,
        const std::string& nameSecondaryBody,
        const std::string& nameBodyToPropagate,
        const std::string& centralBody,
        const simulation_setup::SystemOfBodies& bodies );

//! Propagate CR3BP from CR3BP environment
/*!
 * Propagate CR3BP from CR3BP environment.
 * \param initialTime Initial time for the propagation  [s].
 * \param finalPropagationTime Final time at which the propagation ends [s].
 * \param initialState Initial state of the third, smaller body to be propagated (Cartesian position in velocity in m and m/s)
 * \param integratorSettings Integrator settings for the propagation (the initial time and time-step have to be defined as
 * dimensional times, their conversion to adimensional times is implemented within the function).
 * \param bodies Body map for the CR3BP.
 * \param bodiesCR3BP Name of the two primaries defining the CR3BP.
 * \param stateHistory Propagated states of the CR3BP (returned by reference)
 * \param outputInNormalizedCoordinates Boolean denoting whether output is to be in dimensionless quantities (if true)
 */
void propagateCR3BPFromEnvironment(
        const double initialTime,
        const double finalPropagationTime,
        const Eigen::Vector6d& initialState,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
        const simulation_setup::SystemOfBodies& bodies,
        const std::vector < std::string >& bodiesCR3BP,
        std::map< double, Eigen::Vector6d >& stateHistory,
        const bool outputInNormalizedCoordinates = false );

void propagateCR3BPAndFullDynamicsProblem(
        const double initialTime,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
        const std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > > propagatorSettings,
        const simulation_setup::SystemOfBodies& bodies,
        const std::vector < std::string >& bodiesCR3BP,
        std::map< double, Eigen::Vector6d >& directPropagationResult,
        std::map< double, Eigen::Vector6d >& cr3bpPropagationResult,
        std::map< double, Eigen::VectorXd >& dependentVariableValues );

//! Propagate the CR3BP and the full dynamics problem
/*!
 * Propagate the CR3BP and the full dynamics problem
 * \param initialTime Initial time for the propagation [s].
 * \param finalTime Final time at which the propagation ends [s].
 * \param initialState Initial state of the third, smaller body to be propagated (Cartesian position in velocity in m and m/s)
 * \param integratorSettings Integrator settings for the propagation (the initial time and time-step have to be defined as
 * dimensional times, their conversion to adimensional times is implemented within the function).
 * \param accelerationModelMap Acceleration map for the CR3BP.
 * \param bodiesToPropagate Bodies to be propagated.
 * \param centralBodies Central bodies for the propagation.
 * \param bodies Body Map for the CR3BP.
 * \param bodiesCR3BP Name of the two primaries defining the CR3BP.
 * \param directPropagation Propagated states of the full dynamics problem (returned by reference)
 * \param cr3bpPropagation Propagated states of the CR3BP, converted to dimensional coordinates (returned by reference)
 */
void propagateCR3BPAndFullDynamicsProblem(
        const double initialTime,
        const double finalTime,
        const Eigen::Vector6d& initialState,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
        const basic_astrodynamics::AccelerationMap& accelerationModelMap,
        const std::vector< std::string >& bodiesToPropagate,
        const std::vector< std::string >& centralBodies,
        const simulation_setup::SystemOfBodies& bodies,
        const std::vector < std::string >& bodiesCR3BP,
        std::map< double, Eigen::Vector6d >& directPropagation,
        std::map< double, Eigen::Vector6d >& cr3bpPropagation );

//! Propagate the CR3BP and the full dynamics problem and compute the state difference at the end of the propagation.
/*!
 * Propagate the CR3BP and the full dynamics problem and compute the state difference at the end of the propagation.
 * \param initialTime Initial time for the propagation [s].
 * \param finalTime Final time at which the propagation ends [s].
 * \param initialState Initial state of the third, smaller body to be propagated (Cartesian position in velocity in m and m/s)
 * \param integratorSettings Integrator settings for the propagation (the initial time and time-step have to be defined as
 * dimensional times, their conversion to adimensional times is implemented within the function).
 * \param accelerationModelMap Acceleration map for the CR3BP.
 * \param bodiesToPropagate Bodies to be propagated.
 * \param centralBodies Central bodies for the propagation.
 * \param bodies Body Map for the CR3BP.
 * \param bodiesCR3BP Name of the two primaries defining the CR3BP.
 * \return State difference between the full dynamics problem and the CR3BP at the final propagation time
 * (expressed in inertial cartesian coordinates).
 */
Eigen::Vector6d getFinalStateDifferenceFullPropagationWrtCR3BP(
        const double initialTime,
        const double finalTime,
        const Eigen::Vector6d& initialState,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
        const basic_astrodynamics::AccelerationMap& accelerationModelMap,
        const std::vector< std::string >& bodiesToPropagate,
        const std::vector< std::string >& centralBodies,
        const simulation_setup::SystemOfBodies& bodies,
        const std::vector < std::string >& bodiesCR3BP );

}

}


#endif // TUDAT_PROPAGATION_CR3BP_FULL
