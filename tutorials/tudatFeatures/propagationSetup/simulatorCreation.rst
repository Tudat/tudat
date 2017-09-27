.. _tudatFeaturesSimulatorCreation:

Create a Dynamics Simulator Object
==================================
At the moment, the following :class:`DynamicSimulator` options are available or under development in Tudat:

- Single-arc dynamics simulator.
- Multi-arc dynamics simulator.
- Hybrid dynamics simulator (under development).

These are implemented in derived classes and are discussed below. 

.. class:: SingleArcDynamicsSimulator
   
   This derived class simulates single arc dynamics and its constructor is:

   .. code-block:: cpp

      SingleArcDynamicsSimulator( bodyMap,
      				  integratorSettings, 
      				  propagatorSettings );

   where:

   - :literal:`bodyMap`

      :class:`NamedBodyMap` the map containing al the objects of type :class:`Body` used in the simulation.

   - :literal:`integratorSettings`

      :class:`IntergratorSettings` contains the settings of the integrator used, as discussed in :ref:`tudatFeaturesIntegratorSettings`.

   - :literal:`propagatorSettings`

      :class:`PropagatorSettings` contains the settings that defines how the orbit is propagated, as described in :ref:`tudatFeaturesPropagatorSettings`.

.. class:: MultiArcDynamicsSimulator
   
   This derived class allows the numerical propagation to be performed in an arc-wise manner. It is constructed using:

   .. code-block:: cpp
   
    MultiArcDynamicsSimulator(
            bodyMap,
            integratorSettings,
            propagatorSettings,
            arcStartTimes )

   where:

   - :literal:`arcStartTimes`

      :literal:`std::vector< double >` contains the times at which the separate arcs start.
      

By default, the equations of motion are integrated once the object is created. This can be changed by adding additional arguments to the cosntructors of the :class:`DynamicsSimulator`, as shown below for the :class:`SingleArcDynamicsSimulator`:

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

.. warning:: It is important to ensure that the propagator settings are compatible with the dynamics simulator type selected. Otherwise it will result in an exception being thrown during run-time.



Retrieving the propagation history
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Once the :class:`DynamicsSimulator` object has been created and the equations of motion have been integrated, the propagation history of the selected bodies is stored within the :class:`DynamicsSimulator`. To make use of it using software, such history needs to be retrieved and saved to a file.

If the state propagation history needs to be saved, the following code needs to be placed after the :class:`DynamicsSimulator` object creation:

.. code-block:: cpp

    // Write body propagation history to file.
    writeDataMapToTextFile( dynamicsSimulator.getEquationsOfMotionNumericalSolution( ),
                            "bodyPropagationHistory.dat",
                            outputPath,
                            "",
                            std::numeric_limits< double >::digits10,
                            std::numeric_limits< double >::digits10,
                            "," );

If the dependent variable history needs to be saved, the following code needs to be placed after the :class:`DynamicsSimulator` object creation:

.. code-block:: cpp

    // Write body dependent variable history to file.
    writeDataMapToTextFile( dynamicsSimulator.getDependentVariableHistory( ),
                            "bodyDependentVariableHistory.dat",
                            outputPath,
                            "",
                            std::numeric_limits< double >::digits10,
                            std::numeric_limits< double >::digits10,
                            "," );

These two code snippets will save two :literal:`.dat` files in the folder specified by the :literal:`outputPath`. You can make use of the :literal:`tudat_applications::getOutputPath( )` function to get a folder name relative to the project folder.
