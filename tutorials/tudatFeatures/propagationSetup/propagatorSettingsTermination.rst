.. _tudatFeaturesPropagatorSettingsTermination:

Propagator Settings: Termination Settings
=========================================
The :class:`PropagationTerminationSettings` are a key parameter in the propagation of a body's orbit, since these will determine the computational time and the size of the output file. Depending on the application, the user may want to end the body propagation according to different criteria. Currently, the following :class:`PropagationTerminationSettings` are offered in Tudat:

- Termination once a certain time has passed.
- Termination once a dependent variable meets a certain criterion.
- Termination once multiple criteria are met.

The different types of :class:`PropagationTerminationSettings` are implemented by means of derived classes and custom constructors, as detailed in this page. 

.. class:: PropagationTimeTerminationSettings

   As the name suggests, these settings will cause the propagation to terminate after a certain time has passed. Please note that the simulator will finish the final time-step, which may cause the termination time to be slightly overpassed. The constructor is:

   .. code-block:: cpp

      PropagationTerminationSettings( terminationTime )

   where:

   - :literal:`terminationTime`

      :literal:`double` that dictates the time that must pass before the simulation is terminated.

.. class:: PropagationDependentVariableTerminationSettings

   A more powerful method is the use of a :class:`DependentVariable` to terminate a propagation. Using this class will terminate the propagation once the chosen dependent variable meets a certain criterion. To do so, the user needs to first create the dependent variable that will dictate the termination of the propagation:

   .. code-block:: cpp

      boost::shared_ptr< SingleDependentVariableSaveSettings > terminationDependentVariable =
            boost::make_shared< SingleDependentVariableSaveSettings >(
                altitude_dependent_variable, "Vehicle", "Earth" );

   where in this example the variable created is an :literal:`altitude_dependent_variable` which provides the altitude of :literal:`Vehicle` with respect to :literal:`Earth`.

   The versatility in the use of :class:`PropagationDependentVariableTerminationSettings` to terminate a simulation lies in the large number of dependent variable types available in Tudat. The user is referred to :ref:`tudatFeaturesPropagatorSettingsDependentVariables` for a full list of the available dependent variables as well as their use.

   The constructor for this derived class is:

   .. code-block:: cpp

       PropagationDependentVariableTerminationSettings( terminationDependentVariable,
      							limitValue,
      							useAsLowerLimit )

   where:

   - :literal:`terminationDependentVariable`

      :literal:`boost::shared_ptr< SingleDependentVariableSaveSettings >` contains the dependent variable used for termination, in this example: :literal:`altitude_dependent_variable` as declared above.

   - :literal:`limitValue`

      :literal:`double` that provides the lower-limit (or higher-limit) that the :literal:`terminationDependentVariable` can take before terminating the propagation.

   - :literal:`useAsLowerLimit`

      :literal:`bool` that dictates whether the :literal:`limitValue` is used as a lower-limit (true) or as a higher-limit (false).

.. class:: PropagationHybridTerminationSettings

   It may be possible that the user desires to terminate a propagation according several criteria, where such criteria may or may not be fulfilled simulataneously. The constructor for this derived class is:

   .. code-block:: cpp

      PropagationHybridTerminationSettings( terminationSettingsList, 
      					    fulFillSingleCondition )

   where:

   - :literal:`terminationSettingsList`

      :literal:`std::vector< boost::shared_ptr< PropagationTerminationSettings > >` where each of its elements contains derived classes of :class:`PropagationTerminationSettings`. The desired :class:`PropagationTerminationSettings` can be added by using the :literal:`push_back` method of :literal:`std::vector`, where the pushed elements are objects of the classes discussed above.

   - :literal:`fulFillSingleCondition`

      :literal:`bool` that determines whether the propagation terminates once a single condition is met (true) or whether all conditions must be met (false).


   .. tip::  It is possible to combine both :class:`PropagationTimeTerminationSettings` and :class:`PropagationDependentVariableTerminationSettings`. 
