.. _satelliteRendezVous:

Mission Linking and Optimization
================================

In this tutorial we introduce the features of the Optimization package, including the MissionLinker class and how to set up an optimization problem in Tudat.

 The code for this tutorial is given on GitHub and is also located in the tudatBundle at:

 .. code-block:: cpp

   tudatBundle/tudatExampleApplications/satellitePropagatorExamples/SatellitePropagatorExamples/satelliteRendezVousExample.cpp

The code follows the following example: two satellites are approximately on the same orbit around Earth. The two are separated by a small angle, with the Chaser
falling behind the Target. We need to optimize the initial burn time for the chaser in order to rendez-vous with the target.

These are the initial orbital values for the two satellites:

 .. tabularColumns:: l R R

 +-------------------------------------------+------------+------------+
 | **Parameter**                             | **Target** | **Chaser** |
 +-------------------------------------------+------------+------------+
 | Semi-major axis :math:`a` [m]             | 6782650.0  | 6782650.0  |
 +-------------------------------------------+------------+------------+
 | Eccentricity :math:`e`                    | 0.00050    | 0.00025    |
 +-------------------------------------------+------------+------------+
 | Inclination :math:`i` [°]                 | 51.640     | 51.635     |
 +-------------------------------------------+------------+------------+
 | Argument of perigee :math:`\omega` [°]    | 235.700    | 84.100     |
 +-------------------------------------------+------------+------------+
 | RAAN :math:`\Omega` [°]                   | 23.400     | 23.400     |
 +-------------------------------------------+------------+------------+
 | Initial true anomaly :math:`\theta_0` [°] | 139.800    | 290.500    |
 +-------------------------------------------+------------+------------+
  

From the table one can discern that, since the two spacecraft have approximately the same inclination and exactly the same RAAN, the angular separation between the two is
 
 .. math::

   \Delta u_0 = \omega_T + (\theta_0)_T - (\omega_C + (\theta_0)_C) = 0.9 [^\circ]

which makes up for an arc distance of
 
 .. math::

   \Delta y_0 = a \cdot \Delta \theta_0 \approx 106.604 [\text{km}]

In order to rendez-vous, the Chaser needs to gain velocity w.r.t. the Target. Because of the geometry of the problem, the needed impulse and velocity direction can be approximated
using the Clohessy Whiltshire equations. The approximation sees the Chaser on the line of the velocity vector, in delay w.r.t. the Target:

 .. math::

   y_0 = -\Delta y_0

We decide to perform a burn parallel to the direction of motion, i.e., all the parameters except for :math:`y_0` and :math:`\dot y_0` are set to 0:

 .. math::

   \dot y_0 = \dfrac{y_0 n}{6 \pi}

where :math:`n = \sqrt{\mu / a^3 }` with :math:`\mu` being the orbital parameter of Earth. The fact that :math:`y_0 < 0` means that **the required burn is retrograde**.

The total impulse required is

 .. math::
   
   \Delta V = | \dot y_0 | \approx 6.39207 [\text{m/s}]

The propulsion system of the Chaser is modelled on the Orbital Maneouvering System of the STS. The values are given in the following table:

 .. tabularColumns:: l
 
 +---------------------------------------------+
 | Initial mass :math:`M_0 = 90116.4` [kg]     |
 +---------------------------------------------+
 | Thrust magnitude :math:`T = 26700.0` [N]    |
 +---------------------------------------------+
 | Specific impulse :math:`I_{sp} = 316.0` [s] |
 +---------------------------------------------+

The initial mass represents an empty space-shuttle with enough fuel for orbital maneouver, de-orbiting and a small margin.

Using the Tsiolkowsky equation:

 .. math::
  
   \Delta V = I_{sp} g_0 ln( \Lambda )

and the relation with a constant thrust:

 .. math::

   t_b =  \dfrac {I_{sp} g_0 M_0}{T} \left( 1 - \dfrac{1}{\Lambda} \right)

The total burn time can be estimated to be :math:`t_b \approx 21.552` [s].

The objetive of this example is to optimize this value to minimize the final separation between the spacecraft. In order to do so we have to set-up two simulations in Tudat: one with a thrust acceleration
model for the chaser and one with only the gravitational accelerations due to Earth, to represent the coasting phase. We will use the MissionLinker class
to link the two simulations and define the decision variable and the objective function.

Creating the first mission segment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For the first simulation we need to define the acceleration model in order to include a thrust in the opposite direction of motion for the Chaser.
After defining the accelerations of the target, which only include the spherical armonics of Earth up to degree and order 4, the accelerations of
the Chaser are defined in the following piece of code:

 .. code-block:: cpp

    // Space shuttle OMS thrust (AJ10-190 engine)
    double thrustMagnitude = 26700.0; //N
    double specificImpulse = 316.0; //s

    // Thrust opposite to direction of motion
    boost::shared_ptr< ThrustDirectionGuidanceSettings > thrustDirectionGuidanceSettings =
            boost::make_shared< ThrustDirectionFromStateGuidanceSettings >(
                "Earth", true, true);
    boost::shared_ptr< ThrustEngineSettings > thrustMagnitudeSettings =
            boost::make_shared< ConstantThrustEngineSettings >(
                thrustMagnitude, specificImpulse );

    // Define chaser acceleration model settings.
    std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > accelerationsOfChaser_1;
    accelerationsOfChaser_1[ "Earth" ].push_back( boost::make_shared< SphericalHarmonicAccelerationSettings >( 4, 0 ) );
    accelerationsOfChaser_1[ "Chaser" ].push_back(
                boost::make_shared< ThrustAccelerationSettings >( thrustDirectionGuidanceSettings, thrustMagnitudeSettings) );
    accelerationMap_1[  "Chaser" ] = accelerationsOfChaser_1;
    bodiesToPropagate_1.push_back( "Chaser" );
    centralBodies_1.push_back( "Earth" );

    // Create acceleration models and propagation settings.
    basic_astrodynamics::AccelerationMap accelerationModelMap_1 = createAccelerationModelsMap(
                bodyMap_1, accelerationMap_1, bodiesToPropagate_1, centralBodies_1 );


The initial Cartesian states are translated from the initial orbital parameters in the table.

A mass propagator for the Chaser is included in the propagator settings:

 .. code-block:: cpp

    // Preliminary time termination settings for first simulation
    boost::shared_ptr< PropagationTimeTerminationSettings > timeTerminationSettings_1 =
            boost::make_shared< propagators::PropagationTimeTerminationSettings >( estimatedBurnTime );

    // State propagator for first simulation
    boost::shared_ptr< TranslationalStatePropagatorSettings< double > > translationalPropagatorSettings =
            boost::make_shared< TranslationalStatePropagatorSettings< double > >
            ( centralBodies_1, accelerationModelMap_1, bodiesToPropagate_1, systemInitialState, timeTerminationSettings_1 );

    // Create mass rate model for the chaser
    boost::shared_ptr< MassRateModelSettings > massRateModelSettings =
            boost::make_shared< FromThrustMassModelSettings >( 1 );
    std::map< std::string, boost::shared_ptr< basic_astrodynamics::MassRateModel > > massRateModels;
    massRateModels[ "Chaser" ] = createMassRateModel(
            "Chaser", massRateModelSettings, bodyMap_1, accelerationModelMap_1 );

    // Create mass propagator settings for the chaser
    std::vector< std::string > bodiesWithMassToPropagate;
    bodiesWithMassToPropagate.push_back( "Chaser" );

    Eigen::VectorXd initialBodyMasses = Eigen::VectorXd( 1 );
    initialBodyMasses( 0 ) = chaserMass;

    boost::shared_ptr< MassPropagatorSettings< double > > massPropagatorSettings =
            boost::make_shared< MassPropagatorSettings< double > >(
                    bodiesWithMassToPropagate, massRateModels, initialBodyMasses, timeTerminationSettings_1 );

    // Create multi-type propagator settings (state + chaser mass propagator)
    std::vector< boost::shared_ptr< PropagatorSettings< double > > > propagatorSettingsVector;
                    propagatorSettingsVector.push_back( translationalPropagatorSettings);
                    propagatorSettingsVector.push_back( massPropagatorSettings );

    boost::shared_ptr< MultiTypePropagatorSettings< double > > propagatorSettings_1 =
            boost::make_shared< MultiTypePropagatorSettings< double > >(
                                    propagatorSettingsVector, timeTerminationSettings_1 );

The time termination settings are just preliminary, as we want to define the time of this part of the simulation as the decision variable. Nevertheless,
they need to be included, as the optimizer does not create one by itself.

Afterwards, the integrator settings and the dynamics simulator are created. The chosen fixed step size is very small, so to augment the focus of the time decision variable:

 .. code-block:: cpp

    const double fixedStepSize = 0.001; //s

    // Create integrator settings for first simulation
    boost::shared_ptr< IntegratorSettings< > > integratorSettings_1 =
            boost::make_shared< IntegratorSettings< > >
            ( rungeKutta4, simulationStartEpoch, fixedStepSize );

    // Create dynamics simulator for first simulation
    boost::shared_ptr< SingleArcDynamicsSimulator< > > dynamicsSimulator_1
            = boost::make_shared< SingleArcDynamicsSimulator< > >(
                bodyMap_1, integratorSettings_1, propagatorSettings_1, false, false, false );


We now proceed to create the decision variable. In this case we want to optimize the burn time, i.e. the simulation time.
The boundaries are set to search for the optimum time in a 10 seconds interval :math:`t_b - 5 \text{ s} < t < t_b + 5 \text{ s}` :

 .. code-block::

    // Set the simulation time as decision variable, with boundaries +-5 seconds from
    // the estimated burn time.
    boost::shared_ptr< optimization::SingleDecisionVariableSettings > decisionVariable_1 =
            boost::make_shared< optimization::SingleDecisionVariableSettings >(
                    optimization::simulation_time_decision_variable,
                            estimatedBurnTime - 5.0, estimatedBurnTime + 5.0 );

Since the start epoch is set to 0, we coud also use the option ``optimization::SingleDecisionVariableFromTerminationSettings``, which includes the time
termination settings. It would be set as following:

 .. code-block:: cpp

    boost::shared_ptr< optimization::SingleDecisionVariableFromTerminationSettings > decisionVariable_1 =
            boost::make_shared< optimization::SingleDecisionVariableFromTerminationSettings >(
                    dynamicsSimulator_1, estimatedBurnTime - 5.0, estimatedBurnTime + 5.0 );

In order to link the dynamics simulator and the decision variable, we use the class MissionSegment:

 .. code-block:: cpp

    // Create first mission segment with the dynamics simulator and the decision variable
    boost::shared_ptr< optimization::MissionSegmentSettings > missionSegment_1 =
            boost::make_shared< optimization::MissionSegmentSettings >( dynamicsSimulator_1, 
                    decisionVariable_1 );

Creating the second mission segment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Now we proceed to create the second simulation. 

Both the Target and the Chaser have as acceleration the spherical armonics gravitational field of Earth up to order and degree 4.

Since we want to minimize the separation between the two spacecraft we are going to need to save the distance between the two bodies as
a dependent variable:


 .. code-block:: cpp

    // Create dependent variable save settings to retrieve
    // the distance between the two spacecraft
    std::vector< boost::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariableSaveSettingsVector;
    
    dependentVariableSaveSettingsVector.push_back( boost::make_shared< SingleDependentVariableSaveSettings > (
            relative_distance_dependent_variable, "Chaser", "Target" ) );
    
    boost::shared_ptr< DependentVariableSaveSettings > dependentVariableSaveSettings =
            boost::make_shared< DependentVariableSaveSettings >( dependentVariableSaveSettingsVector, false );


We also want to propagate the simulation until the separation of the two spacecraft is at least 1000 m. We create, therefore, hibrid termination settings, in order to
include a stop at the minimum required separation or at an adeguate amount of time, here defined as the estimated rendez-vous time plus a small margin:

 .. code-block:: cpp

    // Create hybrid termination settings to simulate until estimated
    // rendez vous + margin or until minimum separation

    const double minimumSeparation = 1000; //m

    std::vector< boost::shared_ptr< propagators::PropagationTerminationSettings > > multiTerminationSettings;

    multiTerminationSettings.push_back( boost::make_shared< propagators::PropagationTimeTerminationSettings >(
                                            estimatedBurnTime + estimatedRendezVousTime + 10*fixedStepSize) );

    multiTerminationSettings.push_back( boost::make_shared<
            propagators::PropagationDependentVariableTerminationSettings >(
                    dependentVariableSaveSettingsVector[0], minimumSeparation, true ) );

    boost::shared_ptr< propagators::PropagationHybridTerminationSettings > hybridTerminationSettings =
        boost::make_shared< propagators::PropagationHybridTerminationSettings >( multiTerminationSettings, true );

The estimated rendez-vous time is calculated using the Clohessy-Wiltshire equations as: :math:`t_r = 2 \pi / n`, i.e. the orbital period of the target.

We then proceed to finalize the simulation with the propagator settings, integrator settings and the dynamics simulator:

 .. code-block:: cpp

    // Create propagator settings for second simulation
    boost::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings_2 =
            boost::make_shared< TranslationalStatePropagatorSettings< double > >
                    ( centralBodies_2, accelerationModelMap_2, bodiesToPropagate_2, Eigen::VectorXd::Zero(12),
                            hybridTerminationSettings, cowell, dependentVariableSaveSettings );

    fixedStepSize = 5.0; //s

    // Create integrator settings for second simulation
    boost::shared_ptr< IntegratorSettings< > > integratorSettings_2 =
            boost::make_shared< IntegratorSettings< > >( rungeKutta4, 0.0, fixedStepSize );

    // Create dynamics simulator for second simulation
    boost::shared_ptr< SingleArcDynamicsSimulator< > > dynamicsSimulator_2 =
            boost::make_shared< SingleArcDynamicsSimulator< > >( bodyMap_2, integratorSettings_2,
                    propagatorSettings_2, false );

The initial state is defined as a vector of zeros (``Eigen::VectorXd::Zero(12)``) and the initial time is set to 0, since, during the optimization process, they will be automatically set to the (interpolated) final values of the
previous mission segment. For the second simulation we are allowed to use a larger step size of 5 s, since the simulation is less sensitive (there is no thrust this time).

It is convenient to set the **objective function** as the minimum separation between the two spacecraft, since the defined simulation time is arbitrary it is unlikely that the distance between spacecraft
will reach its minimum value at the final entry. We use the :ref:`ObjectiveFunction` class as following:

 .. code-block:: cpp

    // Use the minimum separation between spacecraft as objective function 

    const int maxNumberOfEvolutions = 30;
    const double objectiveValue = 0.0;

    boost::shared_ptr< optimization::ObjectiveFunctionFromMinOrMaxDependentVariableSettings > objectiveFunction =
            boost::make_shared< optimization::ObjectiveFunctionFromMinOrMaxDependentVariableSettings >(
                    dependentVariableSaveSettingsVector[0], objectiveValue, 0, true,
                            minimumSeparation, maxNumberOfEvolutions );

We used the previously defined ``dependentVariableSaveSettingsVector[0]`` to get the dependent variable save settings of the distance between spacecraft.
The objective value is 0, so that the objective function will trend to converge to that value, but we set the tolerance to ``minimumSeparation``, so that
the optimization process stops when the minimum distance between spacecraft is inferior to 1000 m. Finally, we want to have a maximum of 30 evolutions.

We then define the :ref:`MissionSegment` object of this second simulation to link the dynamics simulator with the objective function:

 .. code-block:: cpp

    // Create mission segment for second simulation containing the second
    // dynamics simulator and the objective function
    boost::shared_ptr< optimization::MissionSegmentSettings > missionSegment_2 =
            boost::make_shared< optimization::MissionSegmentSettings >( dynamicsSimulator_2,
                    objectiveFunction


Linkage and optimization
~~~~~~~~~~~~~~~~~~~~~~~~

Finally, we link together the two mission segments in a chronological order with the :ref:`MissionLinker` class:


 .. code-block:: cpp

    // Create object to link the two mission segments and optimize the mission constraint

    std::vector< boost::shared_ptr< optimization::MissionSegmentSettings > > missionSegments;

    missionSegments.push_back( missionSegment_1 );
    missionSegments.push_back( missionSegment_2 );

    optimization::MissionLinker missionLinker(missionSegments,
            boost::make_shared< optimization::OptimizationSettings >( optimization::global_optimization,
                    32, true, 1 ), false  );

As :ref:`OptimizationSettings`, we are going to use a global optimization with a population of 32 individuals, set the flag that 
stops the simulation at the maximum number of evolutions to ``true`` and the verbosity to 1, as to visualize the results of each
iteration.

You could either set the last flag to ``true`` or use the statement

 .. code-block:: cpp

    // Start optimization
    missionLinker.optimize();

to optimize the configuration.

Since the decision variable is the simulation time of the first mission segment, in order to retrieve the optimum,
which is set automatically at the end of each simulation, you can directly scout the termination settings of the first dynamics simulator:

 .. code-block:: cpp

    std::cout << "\n Calculated burn time = " << timeTerminationSettings_1->terminationTime_ << "s\n";

The same is done with the optimized results of the two simulations:

 .. code-block:: cpp

    std::map< double, Eigen::VectorXd > integrationResult_1 = dynamicsSimulator_1->getEquationsOfMotionNumericalSolution();

    std::map< double, Eigen::VectorXd > integrationResult_2 = dynamicsSimulator_2->getEquationsOfMotionNumericalSolution();


Results
~~~~~~~

The output for this configuration is:

 .. code-block:: cpp

    Initial spacecraft seaparation = 106.604km

    Required impulse = 6.39207m/s

    Estimated burn time = 21.5519s

    Evolution n: 1 Objective Function: 2806.44
    Evolution n: 2 Objective Function: 2323.05
    [...]
    Evolution n: 6 Objective Function: 2316.05
    Evolution n: 7 Objective Function: 2315.48
    [...]
    Evolution n: 11 Objective Function: 2315.48
    Evolution n: 12 Objective Function: 2315.3
    [...]
    Evolution n: 30 Objective Function: 2315.3

     Calculated burn time = 21.3993s

The calculated burn time is the optimized value. As you can see it is very much similar to the estimated one. Moreover, due to the geometry of the
problem the optimizer fails to reach a value of 1000 m of separation, stopping at about 2315.3 m.

An image of the simulation shows the problem with our assumptions:

 .. figure:: images/RendezVousXY.svg
   :scale: 100 %
   :alt: x-y visualization of rendez-vous in target frame.

   Path of chaser in the target-frame after retrograde burn, planar visualization


 .. figure:: images/RendezVousXZ.svg
   :scale: 100 %
   :alt: x-z visualization of rendez-vous in target frame.
   
   Path of chaser in the target-frame after retrograde burn, x-z visualization.

The picures above show the path of the Chaser in the Target frame, where the x-direction is the radial direction from the centre of Earth,
the z-direction is the direction of the angular momentum of the target w.r.t. the centre of Earth and the y-direction completes the
right-hand rule of thumb frame.

One can see that even if the two orbits are similar, the two spacecraft are separated initially in the radial direction
by some 60 km. The burn, opposite to the direction of motion of the Chaser, has also a radial component. Nevertheless, the calculated burn
brings the orbit at just 2315.3 meters from the target at its minimum separation.

See also: `"Clohessy-Wiltshire equations" (PDF). University of Texas at Austin. Retrieved 23 October 2017. <http://www.ae.utexas.edu/courses/ase366k/cw_equations.pdf>`_



