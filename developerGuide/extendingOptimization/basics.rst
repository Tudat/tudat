.. _optimization_basics:

Basics
======

.. _optimization-linearInterpolation:

Linear interpolation
~~~~~~~~~~~~~~~~~~~~

The interpolation used by the linking methods is linear. Assume we have two available epochs
:math:`t_0` and :math:`t_1` and the respective function values :math:`f(t_0)` and :math:`f(t_1)`.
If we want to retrieve the value of :math:`f` at :math:`t_e`, with

    .. math::
         
        t_0 < t_e < t_1

then, assuming that :math:`e = t_1 - t_0` is small enough, the value can be retrieved as following:

    .. math::
        
        f(t_e) = f(t_0) + (f(t_1) - f(t_0))(t_e - t_0)/(t_1 - t_0)

.. _optimization-dynamicPointerCast:

Dynamic pointer cast
~~~~~~~~~~~~~~~~~~~~

Let us consider the class ``DecisionVariableSettings``:

    .. code-block:: cpp
        :caption: :class:`Tudat/Optimization/decisionVariableSettings.h`
        :name: decisionVariableSettings-h

        struct DecisionVariableSettings{

            DecisionVariableSettings( std::vector< boost::shared_ptr< SingleDecisionVariableSettings > >
                    multiDecisionVariableSettings ) : decisionVariableSettings_( multiDecisionVariableSettings )
            { }

            DecisionVariableSettings( boost::shared_ptr< SingleDecisionVariableSettings >
                    singleDecisionVariableSettings )
            {
                decisionVariableSettings_.push_back( singleDecisionVariableSettings );
            }

            std::vector< boost::shared_ptr< SingleDecisionVariableSettings > > decisionVariableSettings_;

        };

The member ``decisionVariableSettings_`` is a vector of ``boost::shared_ptr< SingleDecisionVariableSettings >``,
but any index can be as well initialized as ``boost::shared_ptr< >`` of any of the 
``SingleDecisionVariableSettings`` derived classes.

For example, if the first member of the constructor parameter ``multiDecisionVariableSettings`` is of the type
``boost::shared_ptr< SingleKeplerElementDecisionVariableSettings >`` the value is accepted and stored as
``boost::shared_ptr< SingleDecisionVariableSettings >``.

``SingleKeplerElementDecisionVariableSettings`` contains methods and members that do not belong to 
``SingleDecisionVariableSettings``. In order to access these methods one has to perform a **dynamic pointer cast**.

Let us stick to the example above. That means that ``decisionVariableSettings_[0]`` is stored as 
``boost::shared_ptr< SingleDecisionVariableSettings >`` but contains information of ``SingleKeplerElementDecisionVariableSettings``.

This can be accertained with the flag ``SingleDecisionVariableSettings::decisionVariable_``, which in this case is set
automatically by the constructor of ``SingleKeplerElementDecisionVariableSettings`` to be:

    .. code-block:: cpp

        decisionVariableSettings_[0]->decisionVariable_ == single_kepler_element_decision_variable;

The dynamic cast to retrieve the information of the derived class is then performed as following:


    .. code-block:: cpp

        boost::shared_ptr< SingleKeplerElementDecisionVariableSettings > dynamicCastDecisionVariable =
                boost::dynamic_pointer_cast< SingleKeplerElementDecisionVariableSettings >( decisionVariableSettings_[0] );


Now all the ``SingleKeplerElementDecisionVariableSettings`` are available inside ``dynamicCastDecisionVariable``.

