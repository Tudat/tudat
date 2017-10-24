.. _decisionVariableSettings:

Decision variable settings
==========================

The decision variable settings classes allow the user to define some of the variables of a 
``SingleArcDynamicsSimulator`` and its members as decision variables for an optimization process.

SingleDecisionVariableSettings
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The class ``SingleDecisionVariableSettings`` is the base class which allows to define a single variable as
decision variable.

**Constructor:**

 .. code-block:: cpp

    SingleDecisionVariableSettings( DecisionVariables decisionVariable,
                                    double lowerBoundary, double upperBoundary,
                                    std::string associatedBody = "" );

    SingleDecisionVariableSettings( DecisionVariables decisionVariable,
                                    Eigen::VectorXd& lowerBoundary, Eigen::VectorXd& upperBoundary,
                                    std::string associatedBody = "" );

 * The ``decisionVariable`` can be either of the following enumerations:

   * ``simulation_time_decision_variable``: the simulation time (from start epoch to end epoch) is tuned to obtain the optimal solution. It is a scalar quantity, therefore the first constructor can be used;
   * ``initial_cartesian_state_decision_variable``: the initial Cartesian state is tuned to obtain the optimal value. This is a 6x1 matrix, therefore only the second constructor, with ``Eigen::Vector6d`` boundaries can be used;
   * ``initial_cartesian_velocity_decision_variable``: the initial Cartesian velocity is tuned to obtain the optimal value. This is a 3x1 matrix, therefore only the second constructor, with ``Eigen::Vector3d`` boundaries can be used;
   * ``initial_cartesian_position_decision_variable``: the initial Cartesian position is tuned to obtain the optimal value. This is a 3x1 matrix, therefore only the second constructor, with ``Eigen::Vector3d`` boundaries can be used.

 * The ``lowerBoundary`` and ``upperBoundary`` represent the interval in which the optimum is localized. They can be set either as scalar, with the first constructor, or as vectors,
   with the second constructor, depending on the dimension of the decision variable.

 * The ``associatedBody`` is the name of the body, as determined in the ``bodyMap`` of the dynamics simulator, to which the decision variable applies.
 

SingleDecisionVariableFromTerminationSettings
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This class allows to use one of the propagation termination settings as decision variables. It can also be used to modulate the time of simulation,
if a ``PropagationTimeTerminationSettings`` object is present in the propagator settings.

**Constructor:**
 .. code-block:: cpp
    
    SingleDecisionVariableFromTerminationSettings(
            boost::shared_ptr< propagators::SingleArcDynamicsSimulator< double > > dynamicsSimulator,
            const double lowerBoundary, 
            const double upperBoundary,
            int positionInVectorOfTerminationSettings = 0 )


 * The ``dynamicsSimulator`` is the dynamics simulator object in which the propagator settings object with the termination
   settings is found
 * The ``positionInVectorOfTerminationSettings`` is the index of the termination settings if a ``PropagatorHybridTerminationSettings`` is used.


SingleCartesianComponentDecisionVariableSettings
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This class allows to use one of the initial Cartesian components as decision variable.

**Constructor:**
 .. code-block:: cpp
   
    SingleCartesianComponentDecisionVariableSettings( 
            orbital_elements::CartesianElements cartesianComponent,
            double lowerBoundary, double upperBoundary, 
            std::string associatedBody )


 * The ``cartesianComponent`` can be chosen as either of the following enumerations:
   
   * ``xCartesianPosition``
   * ``yCartesianPosition``
   * ``zCartesianPosition``
   * ``xCartesianVelocity``
   * ``yCartesianVelocity``
   * ``zCartesianVelocity``


SingleKeplerElementDecisionVariableSettings
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This class allows to define one of the initial Kepler elements as the decision variable.


**Constructor:**
 .. code-block:: cpp

    SingleKeplerElementDecisionVariableSettings( 
            orbital_elements::KeplerianElements keplerElement,
            std::string centralBody, double lowerBoundary,
            double upperBoundary, std::string associatedBody ) :

 * The ``keplerElement`` can be chosen as either of the following enumerations:

   * ``semiMajorAxis`` 
   * ``eccentricity``
   * ``inclination``
   * ``argumentOfPeriapsis``
   * ``longitudeOfAscendingNode``
   * ``trueAnomaly``

 * The ``centralBody`` is the name of the massive body containing the orbital parameter used to calculated the
   Kepler elements. Beware that at the moment there is no option to use a different reference frame than the one used in
   the simulation, therefore the ``centralBody`` should be the name of the associated central body of the
   ``associatedBody`` in the propagator settings.



SingleSphericalOrbitalElementDecisionVariableSettings
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This class allows to define one of the initial spherical orbital elements as the decision variable.

**Constructor:**
 .. code-block:: cpp

    SingleSphericalOrbitalElementDecisionVariableSettings( 
            orbital_elements::SphericalOrbitalStateElements sphericalOrbitalElement,
            double lowerBoundary, double upperBoundary, std::string associatedBody )

 * The ``sphericalOrbitalElement`` can be chosen as either of the following enumerations:
   
   * ``radius``
   * ``latitude``
   * ``longitude``
   * ``speed``
   * ``flightPath``
   * ``headingAngle``


DecisionVariableSettings
~~~~~~~~~~~~~~~~~~~~~~~~

You can use this class if you want to set multiple decision variables for one mission segment. 

**Constructor:**
 .. code-block:: cpp

    DecisionVariableSettings( 
        std::vector< boost::shared_ptr< SingleDecisionVariableSettings > > multiDecisionVariableSettings )

 * ``multiDecisionVariableSettings`` is a vector of ``boost::shared_ptr< SingleDecisionVariableSettings >`` objects or any
   of the above mentioned classes.



