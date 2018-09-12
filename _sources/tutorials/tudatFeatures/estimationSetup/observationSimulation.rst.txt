.. _observationSimulation:

Simulating Observations
=======================

.. _creatingObservationSimulators:

Creating the Observation Simulator
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In Tudat, a set of observations are simulated using the :class:`ObservationSimulatorBase` class. An object of this class is used to simulate observations of a single type, for any number of :class:`LinkEnds`. The :class:`ObservationSimulatorBase` has, like many of the classes in Tudat a number of template arguments, one for the scalar type of the observable and one for the time type. Note that the state scalar used in numerical propagation should be equal to the scalar type of the observable to make full use of the functionality. 

There are two ways in which to obtain :class:`ObservationSimulatorBase` objects:

* When creating a :class:`OrbitDeterminationManager` class, a set of :class:`ObservationSimulatorBase` objects are automatically created. Retrieving these is done simply by:

   .. code-block:: cpp

      OrbitDeterminationManager< ObservationScalarType, TimeType > orbitDeterminationManager = ..... //OrbitDeterminationManager object created here.
      
      std::map< ObservableType, boost::shared_ptr< ObservationSimulatorBase< ObservationScalarType, TimeType > > > observationSimulators =
         orbitDeterminationManager.getObservationSimulators( );
         
   This returns the :literal:`observationSimulators` list, which contains an observation simulator for each :literal:`ObservableType` defined in your orbit determination manager.
   
* You can also create an :class:`ObservationSimulatorBase` directly, without using an :literal:`OrbitDeterminationManager` object.    
   
   .. code-block:: cpp

      // Define type of observable
      ObservableType observableType = .... 
       
      // Define observation settings for each require LinkEnds 
      std::map< LinkEnds, boost::shared_ptr< ObservationSettings  > > settingsPerLinkEnds = .... 
      
      // Define environment
      NamedBodyMap bodyMap = ....
      
      // Create observation simulator      
      boost::shared_ptr< ObservationSimulatorBase< ObservationScalarType, TimeType > > observationSimulator = 
         createObservationSimulator( observableType, settingsPerLinkEnds, bodyMap );

In either case, this provides you with an object of type :literal:`ObservationSimulatorBase`.  

Observation Simulation Settings
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The times at which observations are simulated may be defined directly by the user, or they may depend on some scheduling algorithm, which is then used to determine the observation times. The final observation times are determined by a combination of :class:`ObservationSimulationTimeSettings` objects and :literal:`ObservationViabilitySettings` objects. The former allows you to define observation times directly, or an observation schedule algorithm. The latter defines constraints that must be met for an observation to be possible. The viability settings are discussed on the page :ref:`observationViability`.

Each type of observation settings is defined by a dedicated derived class of :class:`ObservationSimulationTimeSettings`. This class has a :literal:`TimeType` argument, which is discussed in more detail :ref:`tudatTemplatingStateTime`. The following classes are presently available:

.. class:: TabulatedObservationSimulationTimeSettings

   The :literal:`TabulatedObservationSimulationTimeSettings` class is used to define a simple list of times at which observations are simulated.

   .. code-block:: cpp

      boost::shared_ptr< TabulatedObservationSimulationTimeSettings< TimeType > > observationSettings =
            boost::make_shared< TabulatedObservationSimulationTimeSettings< TimeType > >( 
                linkEndType, simulationTimes );

   The input is:

   - :literal:`linkEndType`

      A :literal:`LinkEndType` variable denoting the reference link end type for the observation times.
      
      
   - :literal:`simulationTimes`

      A :literal:`std::vector< TimeType >` variable, containing the list of times at which observations are to be simulated.

.. _observationViability:

Observation Viability Setttings
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In many cases, you will not have the list of observation times *a priori*. Instead, the observation times could be a function of the states of the link ends, and depend on a number of constraints that must be satisfied for an observation to be possible. The constraints defined in Tudat are listed in the :literal:`ObservationViabilityType` enum, which can take the following values: 

* :literal:`minimum_elevation_angle`: Minimum elevation angle at a ground station: target must be at least a certain elevation above the horizon.
* :literal:`body_avoidance_angle`: Body avoidance angle: the line-of-sight vector from a link end to a given third body must have an angle w.r.t. the line-of-sight between link ends that is sufficiently large. This constraint is typically used to prevent the Sun from being too close to the field-of-view of the telescope(s).
* :literal:`body_occultation`: Body occultation: the link must not be obscured by a given third body.  For instance: the Moon occulting a link between Earth and Mars.

In Tudat, such constraints are defined by objects of the :literal:`ObservationViabilitySettings` class.

.. class:: ObservationViabilitySettings

   The :literal:`ObservationViabilitySettings` class is used to define a simple list of times at which observations are simulated.

   .. code-block:: cpp

      boost::shared_ptr< ObservationViabilitySettings > observationViabilitySettings =
            boost::make_shared< ObservationViabilitySettings >( 
                observationViabilityType, associatedLinkEnd, stringParameter, doubleParameter );

   The input is:

   - :literal:`observationViabilityType`

      A :literal:`ObservationViabilityType` variable denoting the type of constraint that is to be created
      
      
   - :literal:`associatedLinkEnd`

      A :literal:`std::pair< std::string, std::string >` variable, denoting the link end for which the constraint is to be applied

      .. note:: 
         When leaving the second entry of the :literal:`associatedLinkEnd` empty (for instance :literal:`std::make_pair( "Earth", "" )`, the constraint will be applied for all ground stations on that body.
         
   - :literal:`stringParameter` 
    
      An :literal:`std::string` input parameter defining a property of the constraint. Its meaning is different for different constraint types:
    
       * :literal:`minimum_elevation_angle`: None (stringParameter must be :literal:`""`)
       * :literal:`body_avoidance_angle`: Name of body to which viewing angle should be larger than value defined by :literal:`doubleParameter`
       * :literal:`body_occultation`: Name of body for which occultation is to be taken into account
        
   - :literal:`doubleParameter` 
    
      A :literal:`double` input parameter defining a property of the constraint. Its meaning is different for different constraint types:
    
       * :literal:`minimum_elevation_angle`: Minimum value of elevation angle (in radians) at ground station 
       * :literal:`body_avoidance_angle`: Minimum value of body viewing angle (in radians) of body that is to be avoided.
       * :literal:`body_occultation`: None (doubleParameter must be :literal:`TUDAT_NAN`)

As is the case for many other Tudat functionalities, the actual objects that perform the viability calculcations (of the :literal:`ObservationViabilityCalculator` class) are created from the settings objects as follows:


   .. code-block:: cpp

      // Define environment
      NamedBodyMap bodyMap = .... ;

      // Define link ends for each observable
      std::map< ObservableType, std::vector< LinkEnds > > linkEndsList = .... ;
      
      // Define observation viability settings
      std::vector< boost::shared_ptr< ObservationViabilitySettings > > observationViabilitySettings = .... ;
      
      //  Create viability calculators
      PerObservableObservationViabilityCalculatorList viabilityCalculators = createObservationViabilityCalculators(
                bodyMap, testLinkEndsList, observationViabilitySettings );
      

Where :literal:`PerObservableObservationViabilityCalculatorList` is a typedef for :literal:`std::map< ObservableType, std::map< LinkEnds, std::vector< boost::shared_ptr< ObservationViabilityCalculator > > > >`, which is a list of viability calculators for each set of link ends and observable type.      

.. _observationNoise:

Observation Noise
~~~~~~~~~~~~~~~~~

In addition to the observation biases (see :ref:`observationBiases`), which are part of the observation model and typically deterministic, stochastic noise may be added to the observations when simulating them. 

The interface for observation noise is made general, allowing both time-correlated and time-uncorrelated noise to be added: a function of type :literal:`boost::function< double( const double ) >` must be created. Here, the function input is the current time, and the output the noise value. You are free to define this function in any way you like. Refer to the documentation of :literal:`boost::function` and :literal:`boost::bind` (see :ref:`externalBoost`).

In typical basic simulation studies, time-uncorrelated white noise is used. To easily add this type of noise, you can make use of the Tudat interface to boost probability distributions/random number generation (see :ref:`tudatFeaturesProbabilityDistributions`). As an example, the following will generate a function which generates which noise with a mean of 0.005 and a standard deviationof 0.003.

.. code-block:: cpp

   // Define (arbitrary) noise properties
   double meanValue = 5.0E-3
   double standardDeviation = 3.0E-3;

   // Create noise function
   boost::function< double( ) > inputFreeNoiseFunction = createBoostContinuousRandomVariableGeneratorFunction(
       normal_boost_distribution, boost::assign::list_of( meanValue )( standardDeviation ), 0.0 );
   boost::function< double( const double ) > noiseFunction =
       boost::bind( &utilities::evaluateFunctionWithoutInputArgumentDependency< double, const double >,
          inputFreeNoiseFunction, _1 );
          
You may use a similar approach to use any of the boost distrbutions for noise. Note that the second step, in which the :literal:`evaluateFunctionWithoutInputArgumentDependency` is called, is needed for consistency with the observation noise interface.

.. _generatingObservations:

Generating the observations
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Before discussing in detail how to generate simulated observations, we need to define the manner in which these observations are return. Presently, they are stored in the following complicated data type:

:literal:`std::map< ObservableType, std::map< LinkEnds, std::pair< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >, std::pair< std::vector< TimeType >, LinkEndType > > > >`

The first part of this type :literal:`std::map< ObservableType, std::map< LinkEnds, ... ` denotes that a separate set of observations is generated for each requested observable type and set of link ends. For each of these, the simulated data is stored in the following data type:

:literal:`std::pair< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >, std::pair< std::vector< TimeType >, LinkEndType > >`

This pair contains:

 * A vector with the values of the observables, as a :literal:`Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >` (equal to :literal:`Eigen::VectorXd` when :literal:`ObservationScalarType = double`).
 * Another pair, this time: :literal:`std::pair< std::vector< TimeType >, LinkEndType >`, which first contains the observation times, and second the reference link end type of these observations (e.g. is the time valiud at reception or transmission of the signal).

For observations of size 1, the :literal:`Eigen::Vector` of observations and :literal:`std::vector` of times are the same length. For observations with size larger than 1, however, they are not, with the vector of observations being N times the size of the vector of times (with N the size of a single observable). For instance, for an angular position observable (N=2), entry 0 of the time vector gives the observation time of entry 0 and 1 of the observation vector, entry 1 of the time vector gives the time of entry 2 and 3 of the observation vector, etc.

Using the above, you can create all the required input to generate observations. Note that while the :class:`ObservationSimulatorBase` and :class:`ObservationSimulationTimeSettings` are required for this, the noise function and viability calculators need not be provided (no noise and no observation constraints are then used). The simplest way to generate observations, without noise or viability checks, is by using the following:

.. code-block:: cpp

   // Define times at which to simulate the observations
   std::map< ObservableType, std::map< LinkEnds, boost::shared_ptr< ObservationSimulationTimeSettings< TimeType > > > > observationTimeSettings = .... ;
   
   // Define observation simulator objects
   std::map< ObservableType, boost::shared_ptr< ObservationSimulatorBase< ObservationScalarType, TimeType > > > observationSimulators = .... ;
   
   // Define (arbitrary) noise properties
   FullSimulatedObservationSet = simulateObservations( observationsToSimulate, observationSimulators );

When including checks on the viability of the observations, this must be extended to:

.. code-block:: cpp

   // Define times at which to simulate the observations
   std::map< ObservableType, std::map< LinkEnds, boost::shared_ptr< ObservationSimulationTimeSettings< TimeType > > > > observationTimeSettings = .... ;
   
   // Define observation simulator objects
   std::map< ObservableType, boost::shared_ptr< ObservationSimulatorBase< ObservationScalarType, TimeType > > > observationSimulators = .... ;
   
   // Define viability calculators for observations
   PerObservableObservationViabilityCalculatorList viabilityCalculatorList = .... ; 
   
   // Define (arbitrary) noise properties
   FullSimulatedObservationSet = simulateObservations( observationsToSimulate, observationSimulators, viabilityCalculatorList );

Which will limit the simulated observation set to those that comply with the conditions defined by the :literal:`viabilityCalculatorList`, see THIS PAGE for more details.

Finally, when including noise on the simulated observations, we provide a number of interfaces of varying levels of generality. The interface that provides the greatest degree of freedom is the following:

.. code-block:: cpp

   // Define times at which to simulate the observations
   std::map< ObservableType, std::map< LinkEnds, boost::shared_ptr< ObservationSimulationTimeSettings< TimeType > > > > observationTimeSettings = .... ;
   
   // Define observation simulator objects
   std::map< ObservableType, boost::shared_ptr< ObservationSimulatorBase< ObservationScalarType, TimeType > > > observationSimulators = .... ;
   
   // Define viability calculators for observations
   PerObservableObservationViabilityCalculatorList viabilityCalculatorList = .... ; 
   
   // Define observation noise functions
   std::map< ObservableType, std::map< LinkEnds, boost::function< Eigen::VectorXd( const double ) > > > noiseFunctions = .... ;
   
   // Define (arbitrary) noise properties
   FullSimulatedObservationSet = simulateObservationsWithNoise( observationsToSimulate, observationSimulators, noiseFunctions, viabilityCalculatorList );

Which requires a noise function defined as a :literal:`Eigen::VectorXd` as a function of time (:literal:`const double`), where we use a vector representation of the observation noise to allow noise models to be applied to multi-valued observables (e.g. angular position). However, The :literal:`noiseFunctions` may also be of one of the following:

* :literal:`std::map< ObservableType, std::map< LinkEnds, boost::function< double( const double ) > > >` Here the noise is defined as a single output. If the observable is multi-valued, the same function is called to generate the noise for each of the entries of the observable. Note that the function is called separately for each entry.
* :literal:`std::map< ObservableType, boost::function< Eigen::VectorXd( const double ) > >` Here, the noise is not defineed separately for each set of :literal:`LinkEnds`, only per :literal:`ObservableType`, the same function is used for each set link ends of a given type of observable.
* :literal:`std::map< ObservableType, boost::function< double( const double ) > >` A combination of the previous two input types.
* :literal:`boost::function< double( const double ) >` The same noise function is used for each observable, link ends, and observable entry (for multi-valued observables)