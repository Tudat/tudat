.. _walkthroughsFullPropagationCircularRestrictedThreeBodyProblem:

Full propagation of the circular restricted three body problem
==============================================================
The example described on this page focuses on the comparison of the results obtained with the circular restricted three body problem with those provided by the propagation of the full dynamics problem. The code for this tutorial is available on GitHub, and is located in the Tudat bundle at::

   tudatBundle/tudatExampleApplications/satellitePropagatorExamples/SatellitePropagatorExamples/fullPropagationCircularRestrictedThreeBodyProblem.cpp

As a reminder, the circular restricted three body problem relies on the following assumptions: 
    - The two main bodies are in co-rotating circular orbits about their barycentre. They are referred to as the primary and secondary in the following.
    - The mass of the third body is neglected with respect to the masses of the two main bodies.
    - The only perturbations acting on the massless body are the point-mass gravitational accelerations exerted by the two main bodies.

This example deals with a horseshoe orbit around the L3 point of the Sun-Jupiter system.

Ideal case corresponding to the assumptions of the circular restricted three-body problem
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

First, the initial state of the spacecraft to be propagated is defined, in inertial cartesian coordinates. It corresponds to the initial conditions of a horseshoe orbit around the L3 point of the Sun-Jupiter system (ADD REFERENCE)

.. code-block:: cpp
   
   // Global characteristics of the problem
   double distanceSunJupiter = 778.0e9;

   // Initialise the spacecraft state
   Eigen::Vector6d initialState = Eigen::Vector6d::Zero();
   initialState[0] = - 7.992e11;
   initialState[4] =  -1.29e4;


The integrator settings to be used for the propagation of the full dynamics are then specified:

.. code-block:: cpp

   // Create integrator settings.
   double initialTime = 0.0;
   const double fixedStepSize = 100000.0;
   std::shared_ptr< numerical_integrators::IntegratorSettings< > > integratorSettings =
            std::make_shared < numerical_integrators::IntegratorSettings < > >
            ( numerical_integrators::rungeKutta4, initialTime, fixedStepSize );


Before propagating the full dynamics problem and calculating the circular restricted three body problem solution, the user needs to set up the environment. To this end, the function :literal:`setupBodyMapCR3BP` which automatically defines the body map corresponding to the circular restricted three body problem can be used. It takes as inputs the distance between the primary and secondary bodies and the names of the three bodies involved in the problem. Their ephemerides and gravity fields are directly defined in such a way that they respect the CR3BP assumptions. 

.. code-block:: cpp

   // Create body map.
   std::vector < std::string > bodiesCR3BP;
   bodiesCR3BP.push_back( "Sun" );
   bodiesCR3BP.push_back( "Jupiter" );

   simulation_setup::NamedBodyMap bodyMap = propagators::setupBodyMapCR3BP(
                distanceSunJupiter, "Sun", "Jupiter", "Spacecraft" );

In this problem, the Sun (primary body) and the Jupiter (secondary body) are the co-rotating, main bodies of the CR3BP while the spacecraft is the third, massless body.

The accelerations which are exerted on the spacecraft have also to be defined. In a similar way to what was done for the body map, the function :literal:`setupAccelerationMapCR3BP` directly returns a :literal:`accelerationMap` object with the accelerations corresponding to the CR3BP assumptions (point-mass gravitational attraction from the main bodies only).

.. code-block:: cpp

   // Define propagator settings variables.
    std::vector< std::string > bodiesToPropagate;
    std::vector< std::string > centralBodies;
    bodiesToPropagate.push_back( "Spacecraft" );
    centralBodies.push_back( "SSB" );

    // Create acceleration map.
    basic_astrodynamics::AccelerationMap accelerationModelMap = propagators::setupAccelerationMapCR3BP(
                "Sun", "Jupiter", bodiesToPropagate.at( 0 ), centralBodies.at( 0 ), bodyMap );


The results obtained with the CR3BP analytical solution and those provided by the propagation of the full dynamics problem are calculated by calling the function :literal:`propagateCR3BPAndFullDynamicsProblem`. This function fills the maps which are provided as inputs with the state history of the spacecraft calculated from the CR3BP solution and from the full problem propagation, respectively.

.. code-block:: cpp

    // Calculate the difference between CR3BP and full problem.
    std::map< double, Eigen::Vector6d> fullPropagation;
    std::map< double, Eigen::Vector6d> cr3bpPropagation;

    propagators::propagateCR3BPAndFullDynamicsProblem(initialTime, finalTime, initialState, integratorSettings, accelerationModelMap,
                                                      bodiesToPropagate, centralBodies, bodyMap, bodiesCR3BP, fullPropagation,
                                                      cr3bpPropagation);   

To directly retrieve the state difference between the CR3BP solution and the result of the full problem propagation, the function :literal:`getFinalStateDifferenceFullPropagationWrtCR3BP` can be used and returns the difference in cartesian state at the required final time between the two computational methods.  

.. code-block:: cpp

   Eigen::Vector6d stateDifference = propagators::getFinalStateDifferenceFullPropagationWrtCR3BP(
                initialTime, finalTime, initialState, integratorSettings, accelerationModelMap,
                bodiesToPropagate, centralBodies,bodyMap, bodiesCR3BP );

As the body map and acceleration map have here been defined in such a way that they actually fullfil the CR3BP assumptions, no significant state differences are expected between the CR3BP and the full propagation results.

.. figure:: images/horseshoeOrbit.png
.. figure:: images/differenceCR3BPfullProblem.png


Perturbed case
~~~~~~~~~~~~~~

The previous example has been developed in the ideal case in which the full dynamics problem actually corresponds to the CR3BP and respects its assumptions. However, for real-world applications, such a simple dynamical model is rather unrealistic and the CR3BP solution is actually an approximate solution for which the results of the full problem propagation can significantly differ. In the following example, a more complex and realistic model is considered. 

First of all, the orbits of the two main orbits are neither perfectly circular nor their orbital periods about their barycentre are equal. Instead of this simplified model for their orbits, use can be made of the default settings, which include more realistic ephemerides and gravity fields. The ephemerides are initialised in such a way that the Sun and Jupiter are initially aligned.

.. code-block:: cpp

   // Create body map.
   simulation_setup::NamedBodyMap bodyMapPerturbedCase;

   bodyMapPerturbedCase[ "Jupiter" ] = std::make_shared< simulation_setup::Body >( );
   bodyMapPerturbedCase[ "Jupiter" ]->setEphemeris( std::make_shared< ephemerides::ApproximatePlanetPositions>(
                        ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::jupiter) );
   bodyMapPerturbedCase[ "Jupiter" ]->setGravityFieldModel( simulation_setup::createGravityFieldModel(
                  simulation_setup::getDefaultGravityFieldSettings("Jupiter", TUDAT_NAN, TUDAT_NAN), "Jupiter"));


   // Ensure that the Sun and Jupiter are initially aligned.
   double distanceJupiterBarycentrePerturbedCase = bodyMapPerturbedCase[ "Jupiter" ]->getEphemeris()->getCartesianState(initialTime)
            .segment(0,3).norm();
   double distanceBarycenterSunPerturbedCase = distanceSunJupiter - distanceJupiterBarycentrePerturbedCase;

   double xSun = -( distanceBarycenterSunPerturbedCase / distanceJupiterBarycentrePerturbedCase ) * bodyMapPerturbedCase[ "Jupiter" ]
            ->getEphemeris()->getCartesianState(initialTime)[0];
   double ySun = -( distanceBarycenterSunPerturbedCase / distanceJupiterBarycentrePerturbedCase ) * bodyMapPerturbedCase[ "Jupiter" ]
            ->getEphemeris()->getCartesianState(initialTime)[1];
   double zSun = bodyMapPerturbedCase[ "Jupiter" ]->getEphemeris()->getCartesianState(initialTime)[2];


   bodyMapPerturbedCase[ "Sun" ] = std::make_shared< simulation_setup::Body >( );
   bodyMapPerturbedCase[ "Sun" ]->setEphemeris( std::make_shared< ephemerides::ConstantEphemeris> (
                                                     ( Eigen::Vector6d( ) << xSun, ySun, zSun, 0.0, 0.0, 0.0 ).finished( ),
                                                     frameOrigin, frameOrientation ) );
   bodyMapPerturbedCase[ "Sun" ]->setGravityFieldModel( simulation_setup::createGravityFieldModel(
                  simulation_setup::getDefaultGravityFieldSettings("Sun", TUDAT_NAN, TUDAT_NAN), "Sun"));

Additionally, the accelerations experienced by the spacecraft usually do not restrict themselves to point-mass gravitational attractions from the two main bodies. A more complete set of accelerations can be defined for the spacecraft, as it is done below.

.. code-block:: cpp

   // Set of accelerations experienced by the spacecraft.
   std::map< std::string, std::vector< std::shared_ptr< simulation_setup::AccelerationSettings > > > bodyToPropagateAccelerations;
   bodyToPropagateAccelerations["Sun"].push_back(std::make_shared< simulation_setup::AccelerationSettings >(
                                                                basic_astrodynamics::central_gravity ) );
   bodyToPropagateAccelerations["Jupiter"].push_back(std::make_shared< simulation_setup::AccelerationSettings >(
                                                                basic_astrodynamics::central_gravity ) );

The integrator settings and spacecraft initial time remain the same with respect to those of the ideal case. The calculation of the CR3BP solution and of the results of the full problem propagation are obtained similarly to what was done in the previous example. The difference in cartesian state between the simplified CR3BP solution and the propagation results as a function of time are plotted below.

.. figure:: images/horseshoeOrbitPerturbedCase.png
.. figure:: images/differenceCR3BPfullProblemPerturbedCase.png



