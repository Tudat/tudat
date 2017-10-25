 .. _missionSegmentSettings:

Mission segment settings
========================

A **mission segment** in this framework is a piece of simulation placed in a chain with other pieces of simulation.
Each simulation from the secon on in the chain starts at the end point of the previous simulation.

The simulation is defined by a ``SingleArcDynamicsSimulator< double >`` object, and the mission segment can also have a **decision variable**
and an **objective function**, which allow to optimize the chain using a :ref:`missionLinker`.

The mission segment is here defined by the object ``MissionSegmentSettings``.

MissionSegmentSettings
~~~~~~~~~~~~~~~~~~~~~~

The ``MissionSegmentSettings`` class allows to link a decision variable or an objective function to a dynamics simulator.

**Constructor:**
 .. code-block:: cpp

    MissionSegmentSettings( 
             boost::shared_ptr< propagators::SingleArcDynamicsSimulator< > > dynamicsSimulator,
             boost::shared_ptr< DecisionVariableSettings > decisionVariableSettings = NULL,
             boost::shared_ptr< ObjectiveFunctionSettings > objectiveFunctionSettings = NULL )

   
    MissionSegmentSettings( 
            boost::shared_ptr< propagators::SingleArcDynamicsSimulator< > > dynamicsSimulator,
            boost::shared_ptr< SingleDecisionVariableSettings > singleDecisionVariableSettings,
            boost::shared_ptr< ObjectiveFunctionSettings > objectiveFunctionSettings = NULL )

    
    MissionSegmentSettings( 
            boost::shared_ptr< propagators::SingleArcDynamicsSimulator< > > dynamicsSimulator,
            boost::shared_ptr< ObjectiveFunctionSettings > objectiveFunctionSettings )

 * The ``dynamicsSimulator`` is a ``SingleArcDynamicsSimulator< double >`` object containing the pre-assembled features
   of the simulation;
 * The ``decisionVariableSettings`` is a ``DecisionVariableSettings`` object containing multiple decision variables retrievable from
   the ``dynamicsSimulator``;
 * The ``singleDecisionVariableSettings`` can be used alternatively is only one decision variable is needed;
 * The ``objectiveFunctionSettings`` is the objective function. A mission segment containing this object must
   be placed at the end of the chain.

Following the provided constructors one can create mission segments with the following combinations:

 * Mission segment with **only a dynamics simulator**: if placed at the beginning or at the end of the chain this
   mission segment is not involved in the optimization. It will anyway be linked to the other mission segments and 
   propagated in the optimizaton process if placed between a mission segment with a decision variable and a mission segment
   with an objective function;
 * Mission segment with a ``singleDecisionVariableSettings``: this mission segment cannot be placed at the end of the chain.
   It contains only one decision variable;
 * Mission segment with a ``singleDecisionVariableSettings`` and a ``objectiveFunctionSettings``: this mission segment can be
   placed anywhere in the chain and it is sufficient for the optimization process. It contains only one decision variable;
 * Mission segment with a ``decisionVariableSettings``: this mission segment cannot be placed at the end of the chain.
   It can have multiple decision variables;
 * Mission segment with a ``decisionVariableSettings`` and a ``objectiveFunctionSettings``: this mission segment can be
   placed anywhere in the chain and it is sufficient for the optimization process. It can have multiple decision variables;
 * Mission segment with an ``objectiveFunctionSettings``: this mission segment must be preceded by at least one
   mission segment with a ``singleDecisionVariableSettings`` or a ``decisionVariableSettings`` and without 
   an ``objectiveFunctionSettings``.





