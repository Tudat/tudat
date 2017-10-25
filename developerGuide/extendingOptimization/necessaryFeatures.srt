.. _optimizationNecessaryFeatures:

Ideas for future development
============================

At the moment several features are missing from the optimization package, which should be included to
provide the best Tudat experience for the users.

Linkage between different reference frames
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The linking of two mission segments in the :ref:`extendMissionLinker` does not take care to recognize the two reference frames are
centered at the same body. The program should recognize whether the two reference frames are the same, and, if not,
perform the due translations, maybe aided by the ``ReferenceFrameManager`` class. This part should be implemented in the 
``optimize()`` method of the ``MissionLinker`` class. See file ``Tudat/Optimization/missionLinker.cpp``.
Beware that at the moment **Tudat do not support different frame orientations**, therefore only the centre of propagation can be
modified.

Centre of frame for orbital element decision variables
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The classes ``SingleCartesianComponentDecisionVariableSettings``, ``SingleKeplerElementDecisionVariableSettings`` and
``SingleSphericalOrbitalElementDecisionVariableSettings``, as well as ``SingleDecisionVariableSettings`` using the flags
``initial_cartesian_state_decision_variable``, ``initial_cartesian_velocity_decision_variable`` and 
``initial_cartesian_position_decision_variable`` allow to use the initial orbital states of a body as decision variables.

At the moment the centre of the reference frame around which these orbital states are defined is the current ``centralBodies`` vector
associated index in the propagator settings, although the user should be able to define its own reference frame centre.

These modifications should be implemented both in the structures of the above mentioned classes, in the files 
``Tudat/Optimization/decisionVariableSettings.h`` and ``Tudat/Optimization/decisionVariableSettings.cpp`` and in the methods
``MissionLinker::getInitialOrbitalState( )``, ``MissionLinker::resetOrbitalState( )``, ``MissionLinker::modifyCartesianComponent( )``,
``MissionLinker::modifyKeplerElement( )`` and ``MissionLinker::modifySphericalComponent( )``, found in the file ``Tudat/Optimization/missionSegmentSettings.cpp``.
Also here the challenges are mainly about the frames in which each body is centered. The developer can use the ``ReferenceFrameManager`` class for this purpose.

Pareto front multi-objective optimization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Despite PaGMO allows some methods to create a multi-objective optimization problem, at the moment this option is not yet introduced in
the optimization namespace. The two options allowed in the ``OptimizationSettings`` class are ``local_optimization`` and ``global_optimization``,
each referencing to only one optimization alghorithm provided by PaGMO.


Other optimization settings
~~~~~~~~~~~~~~~~~~~~~~~~~~~

At the moment the ``OptimizationSettings`` class, defined in the file ``Tudat/Optimization/missionLinker.h`` allows for a sparse number of optimization
settings. Some more options should be introduced to allow the user a more flexible optimization capability.


Miscellaneous
~~~~~~~~~~~~~

 * An exact Lambert targeter, starting at a certain orbit around a departure Planet and ending at a certain orbit around a target Planet. This can only be implemented with an optimization process.
 * Other?
