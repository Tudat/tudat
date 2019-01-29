.. _walkthroughsFullPropagationCircularRestrictedThreeBodyProblem:

Full propagation of the circular restricted three body problem
==============================================================
The example described on this page focuses on the comparison of the results obtained with the circular restricted three body problem with those provided by the propagation of the full dynamics problem. The code for this tutorial is available on GitHub, and is located in the Tudat bundle at::

   tudatBundle/tudatExampleApplications/satellitePropagatorExamples/SatellitePropagatorExamples/fullPropagationCircularRestrictedThreeBodyProblem.cpp

As a reminder, the circular restricted three body problem relies on the following assumptions: 
    - The two main bodies are in co-rotating circular orbits about their barycentre. They are referred to as the primary and secondary in the following.
    - The mass of the third body is neglected with respect to the masses of the two main bodies.
    - The only perturbations acting on the massless body are the point-mass gravitational accelerations exerted by the two main bodies.


Ideal case corresponding to the assumptions of the circular restricted three-body problem
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Before propagating the full dynamics problem and calculating the circular restricted three body problem solution, the user needs to set up the environment. To this end, the function :literal:`setupBodyMapCR3BP` which automatically defines the body map corresponding to the circular restricted three body problem can be used. It takes as inputs the distance between the primary and secondary bodies and the names of the three bodies involved in the problem. Their ephemerides and gravity fields are directly defined in such a way that they respect the CR3BP assumptions. 

.. code-block:: cpp

   // Create body map for CR3BP.
   simulation_setup::NamedBodyMap bodyMap = setupBodyMapCR3BP(
                physical_constants::ASTRONOMICAL_UNIT, "Sun", "Earth", "Spacecraft" );

In this problem, the Sun (primary body) and the Earth (secondary body) are the co-rotating, main bodies of the CR3BP while the spacecraft is the third, massless body.

The accelerations which are exerted on the spacecraft have also to be defined. In a similar way to what was done for the body map, the function :literal:`setupAccelerationMapCR3BP` directly returns a :literal:`accelerationMap` object with the accelerations corresponding to the CR3BP assumptions (point-mass gravitational attraction from the main bodies only).

.. code-block:: cpp

   // Define propagator settings variables.
   std::vector< std::string > bodiesToPropagate;
   bodiesToPropagate.push_back( "Spacecraft" );
   std::vector< std::string > centralBodies;
   centralBodies.push_back( "SSB" );

   basic_astrodynamics::AccelerationMap accelerationModelMap = setupAccelerationMapCR3BP(
                "Sun", "Earth", bodiesToPropagate.at( 0 ), centralBodies.at( 0 ), bodyMap );


The integrator settings to be used for the propagation of the full dynamics are then to be specified:

.. code-block:: cpp

   // Define integrator settings
   double initialTime = 0.0;
   double finalTime = 120000000.0;
   double fixedStepSize = 1000.0;
   std::shared_ptr< numerical_integrators::IntegratorSettings< > > integratorSettings =
            std::make_shared < numerical_integrators::IntegratorSettings < > >
            ( numerical_integrators::rungeKutta4, initialTime, fixedStepSize ); 

Finally, the initial state of the spacecraft to be propagated is defined, in inertial cartesian coordinates:

.. code-block:: cpp
   
   // Define initial state of the spacecraft.
   Eigen::Vector6d initialState;
   initialState[0] = 2.991957413820000e+10;
   initialState[1] = 1.295555563704656e+11;
   initialState[2] = 0.0;
   initialState[3] = -2.579433850734350e+04;
   initialState[4] = 5.956947312313238e+03;
   initialState[5] = 0.0;

The results obtained with the CR3BP analytical solution and those provided by the propagation of the full dynamics problem are calculated by calling the function :literal:`propagateCR3BPAndFullDynamicsProblem`. This function fills the maps which are provided as inputs with the state history of the spacecraft calculated from the CR3BP solution and from the full problem propagation, respectively.

.. code-block:: cpp

   // Calculate the CR3BP solution and propagate the full dynamics problem simultaneously
   

To directly retrieve the state difference between the CR3BP solution and the result of the full problem propagation, the function :literal:`getFinalStateDifferenceFullPropagationWrtCR3BP` can be used and returns the difference in cartesian state at the required final time between the two computational methods.  

.. code-block:: cpp

   // Calculate the difference between CR3BP and full problem.
    Eigen::Vector6d stateDifference = getFinalStateDifferenceFullPropagationWrtCR3BP(
                initialTime, finalTime, initialState, integratorSettings, accelerationModelMap,
                bodiesToPropagate, centralBodies,bodyMap, bodiesCR3BP );

As the body map and acceleration map have here been defined in such a way that they actually fullfil the CR3BP assumptions, no significant state differences are expected between the CR3BP and the full propagation results.


Perturbed case
~~~~~~~~~~~~~~

The previous example has been developed in the ideal case in which the full dynamics problem actually corresponds to the CR3BP and respects its assumptions. However, for real-world applications, such a simple dynamical model is rather unrealistic and the CR3BP solution is actually an approximate solution for which the results of the full problem propagation can significantly differ. In the following example, a more complex and realistic model is considered. 

First of all, the orbits of the two main orbits are neither perfectly circular nor their orbital periods about their barycentre are equal. Instead of this simplified model for their orbits, use can be made of the default settings, which include more realistic ephemerides and gravity fields.

.. code-block:: cpp

   // Define body settings and create the body map.


Additionally, the accelerations experienced by the spacecraft usually do not restrict themselves to point-mass gravitational attractions from the two main bodies. A more complete set of accelerations can be defined for the spacecraft, as it is done below.

.. code-block:: cpp

   // Define the acceleration model for the spacecraft to be propagated.


The integrator settings and spacecraft initial time remain the same with respect to those of the ideal case. The calculation of the CR3BP solution and of the results of the full problem propagation are obtained similarly to what was done in the previous example. The difference in cartesian state between the simplified CR3BP solution and the propagation results as a function of time are plotted below.





