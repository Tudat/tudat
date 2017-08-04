.. _tudatFeaturesSimulatorCreation:

Create a Dynamics Simulator Object
==================================
At the moment, the following dynamics simulators are available or under development in Tudat:

- Single-arc dynamics simulator.
- Multi-arc dynamics simulator (under development).
- Hybrid dynamics simulator (under development).

The :literal:`bodyMap`, :literal:`integratorSettigns` and the :literal:`propagatorSettings` are binded together when the :literal:`dynamicsSimulator` is created as shown below:

.. code-block:: cpp

    // Create simulation object and propagate dynamics.
    SingleArcDynamicsSimulator< > dynamicsSimulator(
                bodyMap, integratorSettings, propagatorSettings );

By default, the equations of motion are integrated once the object is created. This can be changed by adding additional arguments, as shown below:

.. code-block:: cpp

    // Create simulation object and propagate dynamics.
    SingleArcDynamicsSimulator< > dynamicsSimulator(
                bodyMap, integratorSettings, propagatorSettings , areEquationsOfMotionToBeIntegrated , clearNumericalSolutions , setIntegratedResult );

where:

- :literal:`areEquationsOfMotionToBeIntegrated`
    Boolean to denote whether equations of motion should be integrated immediately at the end of the contructor or not (default true).
- :literal:`clearNumericalSolutions`
    Boolean to determine whether to clear the raw numerical solution member variables after propagation and resetting ephemerides (default true).
- :literal:`setIntegratedResult`
    Boolean to determine whether to automatically use the integrated results to set ephemerides (default true).

.. warning:: It is important to ensure that the propagator settings are compatible with the dynamics simulator type selected. During otherwise will result in exception being thrown during run-time.
