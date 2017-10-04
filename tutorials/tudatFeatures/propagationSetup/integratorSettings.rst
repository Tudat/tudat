.. _tudatFeaturesIntegratorSettings:

Integrator Settings
===================
As the name suggests, the integrator settings tell the dynamics simulator how to integrate numerically the equations of motion that govern the orbital mechanics to simulate. The :class:`IntegratorSettings` are defined using two derived classes, depending on whether the integrator to be used is a fixed step-size integrator or a variable step-size integrator.


.. class:: IntegratorSettings 
   
   This class is used to define the settings for fixed step-size integration. The constructor for this base class is:

   .. code-block:: cpp

      IntegratorSettings( integratorType,
      			  simulationStartEpoch,
      			  fixedStepSize )

   where:

   - :literal:`integratorType`

      :class:`AvailableIntegrators` which defines the fixed step-size integrator type to be used. Currently the only options available are :literal:`euler` and :literal:`rungeKutta4`.

   - :literal:`simulationStartEpoch`

      :literal:`double` that defines the simulation's start epoch. 

   - :literal:`fixedStepSize`

      :literal:`double` that defines the fixed step-size to be used either by the :literal:`euler` or the :literal:`rungeKutta4` numerical integrator. 

.. class:: RungeKuttaVariableStepSizeSettings
   
   This class is used to define the settings for variable step-size integration. The constructor for this derived class is:

   .. code-block:: cpp
   
      RungeKuttaVariableStepSizeSettings( integratorType,
                            		  initialTime,
                            		  initialTimeStep,
                            		  coefficientSet,
                            		  minimumStepSize,
                            		  maximumStepSize )

   where:

   - :literal:`integratorType`

      :class:`AvailableIntegrators` which defines the fixed step-size integrator type to be used. The only option available is :literal:`rungeKuttaVariableStepSize`.

   - :literal:`initialTime`

      :literal:`double` that defines the simulation's initial time. It must be a :literal:`double` variable-type.
   
   - :literal:`initialTimeStep`

      :literal:`double` that defines the initial step-size to be used either by the :literal:`rungeKuttaVariableStepSize` numerical integrator. It must be a :literal:`double` variable-type. 

   - :literal:`coefficientSet`

      :class:`RungeKuttaCoefficients::CoefficientSets` that defines the coefficient set to be used by the :literal:`rungeKuttaVariableStepSize` numerical integrator. The list of available coefficient sets is given in :ref:`tudatFeaturesIntegrators`.

   - :literal:`minimumStepSize`

      :literal:`double` that defines the minimum step-size that the :literal:`rungeKuttaVariableStepSize` numerical integrator can take. 

   - :literal:`maximumStepSize`

      :literal:`double` that defines the maximum step-size that the :literal:`rungeKuttaVariableStepSize` numerical integrator can take.

.. note:: Aside from the arguments listed in this page, the :class:`IntegratorSettings` class and derived classes described here offer a number of optional arguments. The reader is advised to examine the Doxygen documentation included in the code for further details.

.. warning:: Make sure that a compatible :literal:`integratorType` is selected, otherwise a runtime exception will be thrown.
