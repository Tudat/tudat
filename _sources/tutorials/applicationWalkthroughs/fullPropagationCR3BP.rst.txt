.. _walkthroughsFullPropagationCR3BP:

Full propagation of the CR3BP
=============================
The example described on this page focuses on the comparison of the results obtained with the circular restricted three body problem (CR3BP) with those provided by the propagation of the full dynamics problem. The code for this tutorial is available on GitHub, and is located in the Tudat bundle at::

   tudatBundle/tudatExampleApplications/satellitePropagatorExamples/SatellitePropagatorExamples/fullPropagationCircularRestrictedThreeBodyProblem.cpp

As a reminder, the CR3BP relies on the following assumptions: 
    - The two main bodies are in co-rotating circular orbits about their barycentre. They are referred to as the primary and secondary in the following.
    - The mass of the third body is neglected with respect to the masses of the two main bodies.
    - The only perturbations acting on the massless body are the point-mass gravitational accelerations exerted by the two main bodies.
However, this model never exactly holds in reality. It is therefore necessary to investigate how much the full dynamics problem differs from this simplified formation.
This tutorial takes the propagation of a horseshoe orbit around the L3 point of the Sun-Jupiter system as an example for comparing the CR3BP and full dynamics problem.

Ideal case corresponding to the assumptions of the CR3BP
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

First, the initial state of the spacecraft to be propagated is defined, in inertial cartesian coordinates. It corresponds to the initial conditions of a horseshoe orbit around the L3 point of the Sun-Jupiter system (B. Taylor, D. (1981). Horseshoe periodic orbits in the restricted problem of three bodies for a sun-Jupiter mass ratio. Astronomy and Astrophysics. 103. 288-294.).

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


Before propagating the full dynamics problem and calculating the CR3BP solution, the environment also  needs to be set up. To this end, the function :literal:`setupBodyMapCR3BP` which automatically defines the body map corresponding to the CR3BP assumptions can be used. It takes as inputs the distance between the primary and secondary bodies and the names of the three bodies involved in the problem. Their ephemerides and gravity fields are directly defined in such a way that they respect the CR3BP assumptions. 

.. code-block:: cpp

   // Create body map.
   std::vector < std::string > bodiesCR3BP;
   bodiesCR3BP.push_back( "Sun" );
   bodiesCR3BP.push_back( "Jupiter" );

   simulation_setup::NamedBodyMap bodyMap = propagators::setupBodyMapCR3BP(
                distanceSunJupiter, "Sun", "Jupiter", "Spacecraft" );

In this problem, the Sun (primary body) and Jupiter (secondary body) are the co-rotating, main bodies of the CR3BP while the spacecraft is the third, massless body.

The accelerations which are exerted on the spacecraft also have to be defined. In a similar way to what was done for the body map, the function :literal:`setupAccelerationMapCR3BP` directly returns a :literal:`accelerationMap` object which contains the accelerations corresponding to the CR3BP assumptions (point-mass gravitational attraction from the main bodies only).

.. code-block:: cpp

   // Define propagator settings variables.
    std::vector< std::string > bodiesToPropagate;
    std::vector< std::string > centralBodies;
    bodiesToPropagate.push_back( "Spacecraft" );
    centralBodies.push_back( "SSB" );

    // Create acceleration map.
    basic_astrodynamics::AccelerationMap accelerationModelMap = propagators::setupAccelerationMapCR3BP(
                "Sun", "Jupiter", bodiesToPropagate.at( 0 ), centralBodies.at( 0 ), bodyMap );


The results obtained from the CR3BP solution and those provided by the propagation of the full dynamics problem are calculated by calling the function :literal:`propagateCR3BPAndFullDynamicsProblem`. This function fills the maps :literal:`directPropagationResult` and :literal:`cr3bpPropagationResult` provided as inputs with the state history of the spacecraft calculated from the CR3BP solution and from the full problem propagation respectively.

.. code-block:: cpp

    // Calculate the difference between CR3BP and full problem.
    std::map< double, Eigen::Vector6d> fullPropagation;
    std::map< double, Eigen::Vector6d> cr3bpPropagation;

    propagators::propagateCR3BPAndFullDynamicsProblem(initialTime, finalTime, initialState, integratorSettings, accelerationModelMap,
                                                      bodiesToPropagate, centralBodies, bodyMap, bodiesCR3BP, fullPropagation,
                                                      cr3bpPropagation);   

More precisely, the function :literal:`propagateCR3BPAndFullDynamicsProblem` takes the following variables as input parameters:
   - :literal:`initialTime` Initial time of the propagation, given in seconds.
   - :literal:`finalTime` Time at which the propagation of the full problem and calculation of the CR3BP stop, provided in seconds.
   - :literal:`initialState` Initial state of the body to be propagated within the CR3BP and full dynamics problem (given in cartesian coordinates).
   - :literal:`integratorSettings` Integrator settings to be used for the propagation.
   - :literal:`accelerationModelMap` Map of the accelerations to be taken into account in the propagation of the full dynamics problem.
   - :literal:`bodyToPropagate` Name of the body to be propagated.
   - :literal:`centralBody` Name of the central body of the CR3BP (ie barycenter of the primary and secondary bodies).
   - :literal:`bodyMap` Body map defining the environment of the full dynamics problem. 
   - :literal:`bodiesCR3BP` Vector containing the names of the co-rotating primary and secondary bodies which define the CR3BP.
    
To directly retrieve the state difference between the CR3BP solution and the result of the full problem propagation, the function :literal:`getFinalStateDifferenceFullPropagationWrtCR3BP` can be used and returns the difference in cartesian state at the required final time between the two computational methods.  

.. code-block:: cpp

   Eigen::Vector6d stateDifference = propagators::getFinalStateDifferenceFullPropagationWrtCR3BP(
                initialTime, finalTime, initialState, integratorSettings, accelerationModelMap,
                bodiesToPropagate, centralBodies,bodyMap, bodiesCR3BP );

As the body map and acceleration map have here been defined in such a way that they actually fulfill the CR3BP assumptions, no significant state differences are expected between the CR3BP and the full problem propagation results. Below are first provided the spacecraft trajectories obtained from both the CR3BP solution and from the full problem propagation. The differences between the two orbits are actually so small compared to the Sun-Jupiter distance that the two orbits are not distinguishable from this plot. A second plot gives the position difference between the CR3BP and full propagation as a function of time.

.. figure:: images/horseshoeOrbit.png
.. figure:: images/differenceCR3BPfullProblem.png

The difference between the CR3BP and numerical solution remains small (200 m after 350 years of propagation), with the difference most likely due to different integration errors between the two cases.

Perturbed case
~~~~~~~~~~~~~~

The previous example has been developed in the ideal case in which the full dynamics problem actually corresponds to the CR3BP and respects its assumptions. However, for real-world applications, such a simple dynamical model is rather unrealistic. The CR3BP solution is an approximate solution from which the results of the full problem propagation can significantly differ. In the following example, a more complex and realistic model is considered. 

First of all, the orbits of the two main bodies are neither perfectly circular nor their orbital periods about their barycentre are equal. Instead of this simplified model for their orbits, use can be made of the default settings, which include more realistic ephemerides and gravity fields. In addition, extra perturbations due to thrid-bodies and no-gravitational forces, will perturb the solution away from the idealized CR3BP solution. In this example, we will only consider the latter modification, by adding the third-body perturbations due to Earth, Mars, Venus and Saturn.

.. code-block:: cpp

	NamedBodyMap perturbedBodyMap;
        std::vector< std::string > additionalBodies = { "Earth", "Mars", "Venus", "Saturn" };
        for( unsigned int i = 0; i < additionalBodies.size( ); i++ )
        {

            perturbedBodyMap[ additionalBodies.at( i ) ] = std::make_shared< Body >( );
            perturbedBodyMap[ additionalBodies.at( i ) ]->setEphemeris(
                        std::make_shared< ephemerides::ApproximatePlanetPositions>(
                                additionalBodies.at( i ) ) );
            perturbedBodyMap[ additionalBodies.at( i ) ]->setGravityFieldModel(
                        createGravityFieldModel(
                            std::make_shared< CentralGravityFieldSettings >(
                                spice_interface::getBodyGravitationalParameter(
                                    additionalBodies.at( i ) ) ), additionalBodies.at( i ) ) );
        }

Where we have created body objects for each of these perturbers, with a simplified ephemeris, and a point-mass gravity model. The ephemeris and gravity field of the Sun and Jupiter are copied from the ideal case:

.. code-block:: cpp

        perturbedBodyMap[ "Sun" ] = idealBodyMap[ "Sun" ];
        perturbedBodyMap[ "Jupiter" ] = idealBodyMap[ "Jupiter" ]; 

The calculation of the CR3BP solution and of the results of the full problem propagation are obtained similarly to what was done in the previous example. Taking into account the extra perturbations acting on spacecraft leads the results of the full propagation to differ from those of the CR3BP solution. The trajectory of the spacecraft in the two cases as well, as the difference in cartesian state between the simplified CR3BP solution and the propagation results as a function of time, are plotted below.

.. figure:: images/horseshoeOrbitPerturbedCase.png
.. figure:: images/differenceCR3BPfullProblemPerturbedCase.png

Although the additional perturbations introduce a substantial difference in position, on the order of about 50 million km, it is clear that the spacecraft still follows the horshoe orbit relatively well. Below, the normalized in-plane coordinates of the perturbed case are shown:

.. figure:: images/horseshoeOrbitPerturbedCaseTime.png


Finally, we show the numerical solution of the CR3BP, as well as the perturbed full numerical solution, in unnormalized coordinates below. Note the different scale of the z-axis compared to the in-plane axes.

.. figure:: images/horseshoeUnnormalized.png

