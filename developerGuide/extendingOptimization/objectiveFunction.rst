.. _extendObjectiveFunction:

Objective Function Settings
===========================

How does it work
~~~~~~~~~~~~~~~~

The class ``ObjectiveFunctionSettings`` and its derived are used to retrieve some value from the the last
mission segment after propagation.

The structure of the base class for the objective function segments is the following:

    .. code-block:: cpp
        :caption: :class:`Tudat/Optimization/objectiveFunction.h`
        :name: objectiveFunction-h

        struct ObjectiveFunctionSettings{
        
            ObjectiveFunctionSettings( ObjectiveFunctionType objectiveVariable,
                                       const double objectiveValue, const double tolerance = 1e-3,
                                       const unsigned int maxNumberOfEvolutions = 100 ) :
                objectiveVariable_( objectiveVariable ), objectiveValue_( objectiveValue ), tolerance_( tolerance ),
                maxNumberOfEvolutions_( maxNumberOfEvolutions ), objectiveValueIsSet_( true ){ }

            ObjectiveFunctionSettings( ObjectiveFunctionType objectiveVariable, const bool minimizeOrMaximize = true,
                                       const unsigned int maxNumberOfEvolutions = 100 ) :
                objectiveVariable_(objectiveVariable), minimizeOrMaximize_( minimizeOrMaximize ),
                maxNumberOfEvolutions_( maxNumberOfEvolutions ), objectiveValueIsSet_( false ){ }

            virtual ~ObjectiveFunctionSettings( ){ }

            ObjectiveFunctionType objectiveVariable_;
            double objectiveValue_;
            double tolerance_;
            bool minimizeOrMaximize_;
            unsigned int maxNumberOfEvolutions_;
            bool objectiveValueIsSet_;

        };

The two constructor are used in two different cases: either an objective value is set and in that case the objective function will
converge to that value, or the objective value will just be minimized or maximized accordingly to the ``minimizeOrMaximize`` variable
(``true`` for minimize and ``false`` for maximize).

Types of objective functions currently available
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The current types of objective function are reflected in the enumerator ``ObjectiveFunctionType``.

    .. code-block:: cpp
        :caption: :class:`Tudat/Optimization/objectiveFunction.h`
        :name: objectiveFunction-h

        enum ObjectiveFunctionType
        {
            final_cartesian_component_objective_function = 0,
            final_kepler_orbital_element_objective_function = 1,
            final_spherical_orbital_element_objective_function = 2,
            final_dependent_variable_objective_function = 3,
            min_max_dependent_variable_objective_function = 4,
            user_defined_objective_function = 5
        };

These enumerations are used in the class ``ObjectiveFunctionSettings`` and its derived to define the member
``ObjectiveFunctionSettings::objectiveVariable_``.
This member is not directly defined by the Tudat user, as each enumerator belong to different derived classes with different settings, instead, when creating a 
the derived class the member ``objectiveVariable_`` is set automatically. Here's a
definition of the above defined enumerations:

    * ``final_cartesian_component_objective_function`` is the value of ``objectiveVariable_`` when an object of class ``FinalCartesianComponentObjectiveFunctionSettings`` is created;
    * ``final_kepler_orbital_element_objective_function`` is the value of ``objectiveVariable_`` when an object of class ``FinalKeplerOrbitalElementObjectiveFunction`` is created;
    * ``final_spherical_orbital_element_objective_function`` is the value of ``objectiveVariable_`` when an object of class ``FinalSphericalObitalElementObjectiveFunction`` is created;
    * ``final_dependent_variable_objective_function`` is the value of ``objectiveVariable_`` when an object of class ``FinalDependentVariableObjectiveFunction``;
    * ``min_max_dependent_variable_objective_function`` is the value of ``objectiveVariable_`` when an object of class ``MinMaxDependentVariableObjectiveFunction``;
    * ``user_defined_objective_function`` is the value of ``objectiveVariable_`` when an object of class ``UserDefinedObjectiveFunction`` is created.

You are invited to explore the file ``objectiveFunction.h`` to understand the differences of the various derived classes

Retrieving the objective function value
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See the section :ref:`getObjectiveFunctionValue` for an understading on how the objective function value is retrieved from the ``dynamicsSimulator_`` member of the
``MissionSegmentSettings`` class.




