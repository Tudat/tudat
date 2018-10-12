.. _tudatFeaturesPropagatorSettingsTermination:

Propagator Settings: Termination Settings
=========================================
The :class:`PropagationTerminationSettings` are a key parameter in the propagation of a body's orbit, since these will determine the computational time and the size of the output file. Depending on the application, the user may want to end the body propagation according to different criteria. Currently, the following :class:`PropagationTerminationSettings` are offered in Tudat:

- Termination once a certain time has passed.
- Termination once a certain CPU time is reached.
- Termination once a dependent variable meets a certain criterion.
- Termination once a user-defined function returns ``true``. 
- Termination once multiple criteria are met.

The different types of :class:`PropagationTerminationSettings` are implemented by means of derived classes and custom constructors, as detailed in this page. 

.. class:: PropagationTerminationSettings

   Base class for the termination settings of the propagation. Different options are implemented in the derived classes below.

.. class:: PropagationTimeTerminationSettings

   As the name suggests, these settings will cause the propagation to terminate after a certain time has passed. Please note that the simulator will finish the final time-step, which may cause the termination time to be slightly overpassed. The constructor is:

   .. code-block:: cpp

      PropagationTimeTerminationSettings( terminationTime, terminateExactlyOnFinalCondition )

   where:

   - :literal:`terminationTime`

      :literal:`double` that dictates the time that must pass before the simulation is terminated.

   - :literal:`terminateExactlyOnFinalCondition`

      :literal:`bool` that determines if the integrator will either stop exactly on the final condition (true), or if the integrator will terminate on the first step where it is violated (false), which is the default value.

.. class:: PropagationCPUTimeTerminationSettings

   You may want to make sure that the propagation does not exceed a certain computation time. In this case, you can easily set the termination settings as follows:

   .. code-block:: cpp

      PropagationCPUTimeTerminationSettings( cpuTerminationTime )

   where :literal:`cpuTerminationTime` is the maximum allowed CPU time, after which propagation should be stopped.

.. class:: PropagationDependentVariableTerminationSettings

   A more powerful method is the use of a :class:`DependentVariable` to terminate a propagation. Using this class will terminate the propagation once the chosen dependent variable meets a certain criterion. To do so, the user needs to first create the dependent variable that will dictate the termination of the propagation:

   .. code-block:: cpp

      std::shared_ptr< SingleDependentVariableSaveSettings > terminationDependentVariable =
            std::make_shared< SingleDependentVariableSaveSettings >(
                altitude_dependent_variable, "Vehicle", "Earth" );

   where in this example the variable created is an :literal:`altitude_dependent_variable` which provides the altitude of :literal:`Vehicle` with respect to :literal:`Earth`.

   The versatility in the use of :class:`PropagationDependentVariableTerminationSettings` to terminate a simulation lies in the large number of dependent variable types available in Tudat. The user is referred to :ref:`tudatFeaturesPropagatorSettingsDependentVariables` for a full list of the available dependent variables as well as their use.

   The constructor for this derived class is:

   .. code-block:: cpp

       PropagationDependentVariableTerminationSettings( terminationDependentVariable,
      							limitValue,
      							useAsLowerLimit,
                                                        terminateExactlyOnFinalCondition )

   where:

   - :literal:`terminationDependentVariable`

      :literal:`std::shared_ptr< SingleDependentVariableSaveSettings >` contains the dependent variable used for termination, in this example: :literal:`altitude_dependent_variable` as declared above.

   - :literal:`limitValue`

      :literal:`double` that provides the lower-limit (or higher-limit) that the :literal:`terminationDependentVariable` can take before terminating the propagation.

   - :literal:`useAsLowerLimit`

      :literal:`bool` that dictates whether the :literal:`limitValue` is used as a lower-limit (true) or as a higher-limit (false).

   - :literal:`terminateExactlyOnFinalCondition`

      :literal:`bool` that determines if the integrator will either stop exactly on the final condition (true), or if the integrator will terminate on the first step where it is violated (false), which is the default value.

.. class:: PropagationCustomTerminationSettings

   With this class, you can set a custom function that based on some internal calculations will return whether to stop propagation. The function should take the current time as input (hence a :literal:`double`) and output a boolean (i.e., :literal:`bool`). This boolean should be :literal:`true` when the propagation has to be stopped and :literal:`false` otherwise. The constructor looks like this:

   .. code-block:: cpp

      PropagationCustomTerminationSettings( checkStopCondition )

   where :literal:`checkStopCondition` is the only input and is of type :literal:`std::function< bool( const double ) >`. 

   .. tip::
      In case your custom function requires more inputs (e.g., it may depend on the position of the spacecraft or other variables that are not the current time), you can use :literal:`std::bind` to add more inputs.

   As an example, the case where the state of the spacecraft is added as an input is shown below: 

   .. code-block:: cpp

      std::function< Eigen::Vector6d( ) > spacecraftStateFunction =
              std::bind( &Body::getState, bodyMap.at( "Satellite" ) );
      std::shared_ptr< PropagationTerminationSettings > terminationSettings =
              std::make_shared< PropagationCustomTerminationSettings >(
                  std::bind( &customTerminationFunction, std::placeholders::_1, spacecraftStateFunction ) );

.. class:: PropagationHybridTerminationSettings

   It may be possible that the user desires to terminate a propagation according several criteria, where such criteria may or may not be fulfilled simulataneously. The constructor for this derived class is:

   .. code-block:: cpp

      PropagationHybridTerminationSettings( terminationSettingsList, 
      					    fulFillSingleCondition,
                                            terminateExactlyOnFinalCondition )

   where:

   - :literal:`terminationSettingsList`

      :literal:`std::vector< std::shared_ptr< PropagationTerminationSettings > >` where each of its elements contains derived classes of :class:`PropagationTerminationSettings`. The desired :class:`PropagationTerminationSettings` can be added by using the :literal:`push_back` method of :literal:`std::vector`, where the pushed elements are objects of the classes discussed above.

   - :literal:`fulFillSingleCondition`

      :literal:`bool` that determines whether the propagation terminates once a single condition is met (true) or whether all conditions must be met (false).

   - :literal:`terminateExactlyOnFinalCondition`

      :literal:`bool` that determines if the integrator will either stop exactly on the final condition (true), or if the integrator will terminate on the first step where it is violated (false), which is the default value.

   .. tip::  It is possible to combine both :class:`PropagationTimeTerminationSettings` and :class:`PropagationDependentVariableTerminationSettings`. 

.. note:: 
   
   For both :class:`PropagationCPUTimeTerminationSettings` and :class:`PropagationCustomTerminationSettings` the termination cannot be set to occur exactly on the final condition.
