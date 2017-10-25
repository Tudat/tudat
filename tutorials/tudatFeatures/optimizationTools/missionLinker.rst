 .. _missionLinker:

Mission linker
==============

The **mission linker** is the tool that allows to link together and optimize a chain of **mission segments**.
Here it is represented by the class ``MissionLinker``

MissionLinker
~~~~~~~~~~~~~

**Constructor**
 .. code-block:: cpp

    MissionLinker(  
            std::vector< boost::shared_ptr< MissionSegmentSettings > > missionSegmentSettings,
            boost::shared_ptr< OptimizationSettings > optimizationSettings = boost::make_shared< OptimizationSettings >(),
            const bool optimizeConfiguration = true )

 * The ``missionSegmentSettings`` is the vector of ``MissionSegmentSettings`` objects defining the chain of mission segments to optimize.
 * The ``optimizationSettings`` is a user defined object which allows some degrees of freedom over the optimization process.
 * The ``optimizeConfiguration`` flag can be set to ``true`` to start the optimization right after creating the
   MissionLinker object.
   

**Methods**
 .. code-block:: cpp

    void optimize( void );

This is the method to call if you want to optimize the configuration. All the optimum values as well as the results
of the optimal propagation are subsequently stored in the dynamics simulators of each mission segment.

OptimizationSettings
~~~~~~~~~~~~~~~~~~~~

This class allows the user some degrees of freedom over the optimization process.

**Constructor**
 .. code-block:: cpp

    OptimizationSettings( OptimizationTypes optimizationType = global_optimization,
                          const int populationSize = 32, 
                          const bool stopAtMaxNumberOfEvolutions = true,
                          const int verbosity = 0 )

 * The ``optimizationType`` can be set to either of the following enumerations:

   * ``global_optimization`` will use an evolutionary algorithm to find the optimum value (generally slower, but reliable),
   * ``local_optimization`` will use a local derivative algorithm to find the optimum value (faster, less accurate);

 * The ``populationSize`` determines the number of individual in a population. The larger the number, the larger the possibilities;
   to converge in less iterations to the target value, although the optimization takes a proportional amount of time;
 * The ``stopAtMaxNumberOfEvolutions`` flag can be set to ``false`` if you want your optimization to only stop if convergence is reached (not advisable);
 * The ``verbosity`` can be set to:

   * 0: no message will be shown during the optimization,
   * n > 0: it will show the value of the objective function in a message every n iterations.

Example
~~~~~~~

See an example that uses the ``MissionLinker`` at the page :ref:`satelliteRendezVous`.
