.. _tudatFeaturesIntegratorSettings:

Integrator Settings
===================
As the name suggests, the integrator settings tell the dynamics simulator how to integrate numerically the equations of motion that govern the orbital mechanics to simulate. The :class:`IntegratorSettings` are defined using four derived classes, depending on whether the integrator to be used is a fixed step-size integrator or a variable step-size integrator.


.. class:: IntegratorSettings 
   
   This class is used to define the settings for fixed step-size integration. The constructor for this base class is:

   .. code-block:: cpp

      IntegratorSettings< TimeType >( integratorType,
      			              simulationStartEpoch,
      			              fixedStepSize )

   where:

   - :literal:`TimeType`
   
      Template argument used to set the precision of the time, in general :literal:`double` is used. For some application where a high precision is required this can be changed to e.g. :literal:`long double`. 

   - :literal:`integratorType`

      :class:`AvailableIntegrators` which defines the fixed step-size integrator type to be used. Currently the only options available are :literal:`euler` and :literal:`rungeKutta4`.

   - :literal:`simulationStartEpoch`

      :literal:`TimeType` that defines the simulation's start epoch. 

   - :literal:`fixedStepSize`

      :literal:`TimeType` that defines the fixed step-size to be used either by the :literal:`euler` or the :literal:`rungeKutta4` numerical integrator. 

.. class:: RungeKuttaVariableStepSizeSettings
   
   This class is used to define the settings for the Runge-Kutta variable step-size integration methods, where the error tolerances are defined as a scalar (i.e., each state element has the same tolerance). The constructor for this derived class is:

   .. code-block:: cpp
   
      RungeKuttaVariableStepSizeSettings< TimeType >(
            simulationStartEpoch,
            initialTimeStep,
            coefficientSet,
            minimumStepSize,
            maximumStepSize,
            relativeErrorTolerance,
            absoluteErrorTolerance,
            saveFrequency,
            assessPropagationTerminationConditionDuringIntegrationSubsteps,
            safetyFactorForNextStepSize,
            maximumFactorIncreaseForNextStepSize,
            minimumFactorDecreaseForNextStepSize )

   where:

    - :literal:`TimeType`
   
      Template argument used to set the precision of the time, in general :literal:`double` is used. For some application where a high precision is required this can be changed to e.g. :literal:`long double`.

   - :literal:`simulationStartEpoch`

      :literal:`TimeType` that defines the simulation's initial time. It must be a :literal:`double` variable-type.
   
   - :literal:`initialTimeStep`

      :literal:`TimeType` that defines the initial step-size to be used either by the :literal:`rungeKuttaVariableStepSize` numerical integrator. It must be a :literal:`double` variable-type. 

   - :literal:`coefficientSet`

      :class:`RungeKuttaCoefficients::CoefficientSets` that defines the coefficient set to be used by the :literal:`rungeKuttaVariableStepSize` numerical integrator. The list of available coefficient sets is given in :ref:`tudatFeaturesIntegrators`.

   - :literal:`minimumStepSize`

      :literal:`TimeType` that defines the minimum step-size that the :literal:`rungeKuttaVariableStepSize` numerical integrator can take. 

   - :literal:`maximumStepSize`

      :literal:`TimeType` that defines the maximum step-size that the :literal:`rungeKuttaVariableStepSize` numerical integrator can take.

   - :literal:`relativeErrorTolerance`

      :literal:`TimeType` that defines the relative error tolerance for step size control of the :literal:`rungeKuttaVariableStepSize` numerical integrator.

   - :literal:`absoluteErrorTolerance`

      :literal:`TimeType` that defines the absolute error tolerance for step size control of the :literal:`rungeKuttaVariableStepSize` numerical integrator.

   - :literal:`saveFrequency`

      Frequency at which to save the numerical integrated states. For instance, you may want to save one every 15 time steps, to give an output that is less demanding in terms of storage (in this case 15 would be the :literal:`saveFrequency`). The default value is 1.

   - :literal:`assessPropagationTerminationConditionDuringIntegrationSubsteps`

      Whether the propagation termination conditions should be evaluated during the intermediate sub-steps of the integrator (``true``) or only at the end of each integration step (``false``). The default value is ``false``.

   - :literal:`safetyFactorForNextStepSize`

      Safety factor for step size control. The default value is 0.8.

   - :literal:`maximumFactorIncreaseForNextStepSize`

      Maximum increase factor in time step in subsequent iterations. The default value is 4.0.

   - :literal:`minimumFactorDecreaseForNextStepSize`

      Minimum decrease factor in time step in subsequent iterations. The default value is 0.1.

.. class:: RungeKuttaVariableStepSizeSettingsScalarTolerances

.. note:: The :class:`RungeKuttaVariableStepSizeSettings` class is actually a shorthand for the alias :class:`RungeKuttaVariableStepSizeSettingsScalarTolerances`, for compatibility with the previous definition of the Runge-Kutta variable step-size integrator.

.. class:: RungeKuttaVariableStepSizeSettingsVectorTolerances
   
   This class is used to define the settings for the Runge-Kutta variable step-size integration methods, where the error tolerances are defined as a vector (i.e., you could set a different absolute tolerance for position and velocity, if the propagated state is expressed in Cartesian elements). The constructor for this derived class is:

   .. code-block:: cpp
   
      RungeKuttaVariableStepSizeSettingsVectorTolerances< TimeType, StateType >(
            simulationStartEpoch,
            initialTimeStep,
            coefficientSet,
            minimumStepSize,
            maximumStepSize,
            relativeErrorTolerance,
            absoluteErrorTolerance,
            saveFrequency,
            assessPropagationTerminationConditionDuringIntegrationSubsteps,
            safetyFactorForNextStepSize,
            maximumFactorIncreaseForNextStepSize,
            minimumFactorDecreaseForNextStepSize )

   where most of the input variables are the same as for the previous constructor, except for the following:

    - :literal:`StateType`
   
      Template argument used to set the format of the state, in general :literal:`Eigen::VectorXd` is used. For applications where covariance propagation is also performed, this may be :literal:`Eigen::MatrixXd`. One can also change the precision of the state scalar, such as in :literal:`Eigen::VectorXld`, where :literal:`long double` is used instead of :literal:`double`.

   - :literal:`relativeErrorTolerance`

      :literal:`StateType` that defines the relative error tolerance for each state entry, for step size control of the :literal:`rungeKuttaVariableStepSize` numerical integrator.

   - :literal:`absoluteErrorTolerance`

      :literal:`StateType` that defines the absolute error tolerance for each state entry, for step size control of the :literal:`rungeKuttaVariableStepSize` numerical integrator.

.. class:: BulirschStoerIntegratorSettings
   
   This class is used to define the settings for variable step-size integration using the Bulirsch-Stoer method. The constructor for this derived class is:

   .. code-block:: cpp
   
      BulirschStoerIntegratorSettings< TimeType >( initialTime,
                            		           initialTimeStep,
                            		           extrapolationSequence,
                            		           maximumNumberOfSteps,
                            		           minimumStepSize,
					           maximumStepSize,
                            		           relativeErrorTolerence,
                            		           absoluteErrorTolerence )

   where:

    - :literal:`TimeType`
   
      Template argument used to set the precision of the time, in general :literal:`double` is used. For some application where a high precision is required this can be changed to e.g. :literal:`long double`. 


   - :literal:`initialTime`

      :literal:`TimeType` that defines the simulation's initial time. It must be a :literal:`double` variable-type.
   
   - :literal:`initialTimeStep`

      :literal:`TimeType` that defines the initial step-size to be used either by the :literal:`BulirschStoerIntegrator` numerical integrator. It must be a :literal:`double` variable-type. 

   - :literal:`extrapolationSequence`
      	
      :literal:`ExtrapolationMethodStepSequences` that defines the extrapolation sequence that is used for the :literal:`BulirschStoerIntegrator` numerical integrator.

   - :literal:`maximumNumberOfSteps`

      Number of integrations that are used for a single extrapolation. It must be a  :literal:`int` variable-type.
  
   - :literal:`minimumStepSize`

      :literal:`TimeType` that defines the minimum step-size that the :literal:`BulirschStoerIntegrator` numerical integrator can take. 

   - :literal:`maximumStepSize`

      :literal:`TimeType` that defines the maximum step-size that the :literal:`BulirschStoerIntegrator` numerical integrator can take.

   - :literal:`relativeErrorTolerance`

      :literal:`TimeType` that defines the relative error tolerance for step size control of the :literal:`BulirschStoerIntegrator` numerical integrator.

   - :literal:`absoluteErrorTolerance`

      :literal:`TimeType` that defines the absolute error tolerance for step size control of the :literal:`BulirschStoerIntegrator` numerical integrator.

.. class:: AdamsBashforthMoultonSettings
   
   This class is used to define the settings for variable step-size integration using the Adams-Bashfort-Moulton method. The constructor for this derived class is:

   .. code-block:: cpp
   
      AdamsBashforthMoultonSettings< TimeType >( initialTime,
                            		         initialTimeStep,
                            		         minimumStepSize,
					         maximumStepSize,
                            		         relativeErrorTolerence,
                            		         absoluteErrorTolerence,
					         minimumOrder,
					         maximumOrder )

   where:

    - :literal:`TimeType`
   
      Template argument used to set the precision of the time, in general :literal:`double` is used. For some application where a high precision is required this can be changed to e.g. :literal:`long double`. 


   - :literal:`initialTime`

      :literal:`TimeType` that defines the simulation's initial time. It must be a :literal:`double` variable-type.
   
   - :literal:`initialTimeStep`

      :literal:`TimeType` that defines the initial step-size to be used either by the :literal:`AdamsBashforthMoultonIntegrator` numerical integrator. It must be a :literal:`double` variable-type. 
  
   - :literal:`minimumStepSize`

      :literal:`TimeType` that defines the minimum step-size that the :literal:`AdamsBashforthMoultonIntegrator` numerical integrator can take. 

   - :literal:`maximumStepSize`

      :literal:`TimeType` that defines the maximum step-size that the :literal:`AdamsBashforthMoultonIntegrator` numerical integrator can take.

   - :literal:`relativeErrorTolerance`

      :literal:`TimeType` that defines the relative error tolerance for step size control of the :literal:`AdamsBashforthMoultonIntegratorr` numerical integrator.

   - :literal:`absoluteErrorTolerance`

      :literal:`TimeType` that defines the absolute error tolerance for step size control of the :literal:`AdamsBashforthMoultonIntegrator` numerical integrator.

   - :literal:`minimumOrder`

      The minimum order of the integrator, the default value is 6. It must be a :literal:`int` variable-type.  

   - :literal:`maximumOrder`

      The maximum order of the integrator, the default value is 11. It must be a :literal:`int` variable-type.  


.. note:: Aside from the arguments listed in this page, the :class:`IntegratorSettings` class and derived classes described here offer a number of optional arguments. The reader is advised to examine the Doxygen documentation included in the code for further details.

.. warning:: Make sure that a compatible :literal:`integratorType` is selected, otherwise a runtime exception will be thrown.
