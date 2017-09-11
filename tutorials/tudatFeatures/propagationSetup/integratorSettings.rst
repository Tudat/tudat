.. _tudatFeaturesIntegratorSettings:

Integrator Settings
===================
As the name suggests, the integrator settings tell the dynamics simulator how to integrate numerically the equations of motion that govern the orbital mechanics to simulate. The integrator settings are defined using two derived classes, depending on whether the integrator to be used is a fixed step-size integrator or a variable step-size integrator.

Propagation with a fixed step-size integrator
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
In order to use a fixed step-size integrator, the following constructor needs to be created in your :literal:`main` function. The constructor is created as a :literal:`shared_ptr` where the arguments to the constructor are passed by using :literal:`make_shared`:

.. code-block:: cpp

    boost::shared_ptr< IntegratorSettings< > > integratorSettings =
                boost::make_shared< IntegratorSettings< > >
                ( integratorType, simulationStartEpoch, fixedStepSize );

where the following arguments need to be passed:

- :literal:`integratorType`
    Defines the fixed step-size integrator type to be used. Currently the only options available are :literal:`euler` and :literal:`rungeKutta4`.
- :literal:`simulationStartEpoch`
    Defines the simulation's start epoch. It must be a :literal:`double` variable-type.
- :literal:`fixedStepSize`
    Defines the fixed step-size to be used either by the :literal:`euler` or the :literal:`rungeKutta4` numerical integrator. It must be a :literal:`double` variable-type.

Propagation with a variable step-size integrator
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
In order to use a variable step-size integrator, the following constructor needs to be created in your :literal:`main` function. Similarly to the fixed-step size integrator settings, this is done by means of a boost pointer libraries:

.. code-block:: cpp

    boost::shared_ptr< IntegratorSettings< > > integratorSettings =
                boost::make_shared< RungeKuttaVariableStepSizeSettings< > >
                            integratorType,
                            initialTime,
                            initialTimeStep,
                            coefficientSet,
                            minimumStepSize
                            maximumStepSize );

where the following arguments need to be passed:

- :literal:`integratorType`
    Defines the fixed step-size integrator type to be used. The only option available is :literal:`rungeKuttaVariableStepSize`.
- :literal:`initialTime`
    Defines the simulation's initial time. It must be a :literal:`double` variable-type.
- :literal:`initialTimeStep`
    Defines the initial step-size to be used either by the :literal:`rungeKuttaVariableStepSize` numerical integrator. It must be a :literal:`double` variable-type.
- :literal:`coefficientSet`
    Defines the coefficient set to be used by the :literal:`rungeKuttaVariableStepSize` numerical integrator. The list of available coefficient sets is given in :ref:`tudatFeaturesIntegrators`.
- :literal:`minimumStepSize`
    Defines the minimum step-size that the :literal:`rungeKuttaVariableStepSize` numerical integrator can take. It must be a :literal:`double` variable-type.
- :literal:`maximumStepSize`
    Defines the maximum step-size that the :literal:`rungeKuttaVariableStepSize` numerical integrator can take. It must be a :literal:`double` variable-type.

.. note:: Aside from the arguments listed in this page, the :class:`IntegratorSettings` derived classes offer a number of optional arguments. The reader is advised to examine the Doxygen documentation included in the code for further details.

.. warning:: Make sure that a compatible :literal:`integratorType` is selected, otherwise a runtime exception will be thrown.
