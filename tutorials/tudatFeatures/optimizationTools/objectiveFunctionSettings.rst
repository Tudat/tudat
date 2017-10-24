.. _objectiveFunctionSettings:

Objective function settings
===========================

The objective function settings allows to define what is the objective function of the optimization process. Several classes can be
used, each focused on a certain type of objective function.

 .. warning:: At the moment only single objective function optimization is allowed.

FinalCartesianComponentObjectiveFunctionSettings
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This class allows to use one of the final cartesian components of a propagated body as objective function.

**Constructor:**
 .. code-block:: cpp

    FinalCartesianComponentObjectiveFunctionSettings( 
            orbital_elements::CartesianElements cartesianComponent,
            const std::string associatedBody, 
            const double objectiveValue, 
            const double tolerance = 1e-3,
            const unsigned int maxNumberOfEvolutions = 100 )


    FinalCartesianComponentObjectiveFunctionSettings( 
            orbital_elements::CartesianElements cartesianComponent,
            const std::string associatedBody, 
            const bool minimizeOrMaximize = true, 
            const unsigned int maxNumberOfEvolutions = 100 )


The first constructor allows to define an objective value towards which the objective function converges. The second constructor will
continue to optimize the configuration until the maximum number of evolutions is reached.

 * The ``cartesianComponent`` can be chosed as either of the following enumerations:
   
   * ``xCartesianPosition``
   * ``yCartesianPosition``
   * ``zCartesianPosition``
   * ``xCartesianVelocity``
   * ``yCartesianVelocity``
   * ``zCartesianVelocity``

 * The ``associatedBody`` is the name of the propagated body, as determined in the ``bodyMap`` of the dynamics simulator to which the 
   final cartesian component belongs.
 * The ``objectiveValue`` is the value towards which the optimization converges.
 * The ``tolerance`` is the error between the final value and the objective function for which the optimization is allowed to stop.
 * The ``maxNumberOfEvolutions`` is the maximum number of iterations allowed in the optimization process if the convergence is not reached.
 * The ``minimizeOrMaximize`` flag is set to:
   
   * ``true`` if you want to minimize the objective function,
   * ``false`` if you want to maximize the objective function.

FinalKeplerElementObjectiveFunctionSettings
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This class allows to use one of the final Kepler orbital elements of a propagated body as objective function.

**Constructor:**
 .. code-block:: cpp

    FinalKeplerElementObjectiveFunctionSettings( 
            orbital_elements::KeplerianElements keplerElement,
            const std::string associatedBody, 
            std::string centralBody, 
            const double objectiveValue,
            const double tolerance = 1e-3, 
            const unsigned int maxNumberOfEvolutions = 100 )

    FinalKeplerElementObjectiveFunctionSettings( 
            orbital_elements::KeplerianElements keplerElement,
            const std::string associatedBody, 
            std::string centralBody, 
            const bool minimizeOrMaximize = true,
            const unsigned int maxNumberOfEvolutions = 100 )

The first constructor allows to define an objective value towards which the objective function converges. The second constructor will
continue to optimize the configuration until the maximum number of evolutions is reached.

 * The ``keplerElement`` can be chosen as either of the following enumerations:

   * ``semiMajorAxis`` 
   * ``eccentricity``
   * ``inclination``
   * ``argumentOfPeriapsis``
   * ``longitudeOfAscendingNode``
   * ``trueAnomaly``

 * The ``associatedBody`` is the name of the propagated body, as determined in the ``bodyMap`` of the dynamics simulator to which the 
   final cartesian component belongs.
 * The ``centralBody`` is the name of the massive body containing the orbital parameter used to calculated the
   Kepler elements. Beware that at the moment there is no option to use a different reference frame than the one used in
   the simulation, therefore the ``centralBody`` should be the name of the associated central body of the
   ``associatedBody`` in the propagator settings.
 * The ``objectiveValue`` is the value towards which the optimization converges.
 * The ``tolerance`` is the error between the final value and the objective function for which the optimization is allowed to stop.
 * The ``maxNumberOfEvolutions`` is the maximum number of iterations allowed in the optimization process if the convergence is not reached.
 * The ``minimizeOrMaximize`` flag is set to:
   
   * ``true`` if you want to minimize the objective function,
   * ``false`` if you want to maximize the objective function.

FinalSphericalOrbitalElementObjectiveFunctionSettings
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This class allows to use one of the final spherical orbital elements of a propagated body as objective function.

**Constructor:**
 .. code-block:: cpp

    FinalSphericalOrbitalElementObjectiveFunctionSettings( 
            orbital_elements::SphericalOrbitalStateElements spericalOrbitalElement,
            const std::string associatedBody, 
            const double objectiveValue, 
            const double tolerance = 1e-3,
            const unsigned int maxNumberOfEvolutions = 100 )

    FinalSphericalOrbitalElementObjectiveFunctionSettings( 
            orbital_elements::SphericalOrbitalStateElements spericalOrbitalElement,
            const std::string associatedBody, 
            const bool minimizeOrMaximize = true,
            const unsigned int maxNumberOfEvolutions = 100 )

The first constructor allows to define an objective value towards which the objective function converges. The second constructor will
continue to optimize the configuration until the maximum number of evolutions is reached.

 * The ``sphericalOrbitalElement`` can be chosen as either of the following enumerations:
   
   * ``radius``
   * ``latitude``
   * ``longitude``
   * ``speed``
   * ``flightPath``
   * ``headingAngle``

 * The ``associatedBody`` is the name of the propagated body, as determined in the ``bodyMap`` of the dynamics simulator to which the 
   final cartesian component belongs.
 * The ``objectiveValue`` is the value towards which the optimization converges.
 * The ``tolerance`` is the error between the final value and the objective function for which the optimization is allowed to stop.
 * The ``maxNumberOfEvolutions`` is the maximum number of iterations allowed in the optimization process if the convergence is not reached.
 * The ``minimizeOrMaximize`` flag is set to:
   
   * ``true`` if you want to minimize the objective function,
   * ``false`` if you want to maximize the objective function.

ObjectiveFunctionFromFinalDependentVariableSettings
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This class allows to use the final value of one of the dependent variables saved by the propagator as objective function.

 .. warning:: You need to have set the DependentVariableSaveSettings inside the PropagatorSettings in order to use this class!

**Constructor:**
 .. code-block:: cpp

    ObjectiveFunctionFromFinalDependentVariableSettings(
            boost::shared_ptr< propagators::SingleDependentVariableSaveSettings > dependentVariableSaveSettings,
            const double objectiveValue, 
            const int indexInVectorialVariable = 0, 
            const double tolerance = 1e-3,
            const unsigned int maxNumberOfEvolutions = 100 )

    ObjectiveFunctionFromFinalDependentVariableSettings(
            boost::shared_ptr< propagators::SingleDependentVariableSaveSettings > dependentVariableSaveSettings,
            const bool minimizeOrMaximize, 
            const int indexInVectorialVariable = 0,
            const unsigned int maxNumberOfEvolutions = 100 )

The first constructor allows to define an objective value towards which the objective function converges. The second constructor will
continue to optimize the configuration until the maximum number of evolutions is reached.

 * The ``dependentVariableSettings`` is the ``SingleDependentVariableSaveSettings`` object belonging to the saved variable whose final value 
   you want to use as objective function.
 * The ``objectiveValue`` is the value towards which the optimization converges.
 * The ``indexInVectorialVariable`` is the index pointing to the position that you need in the dependent variable in case the latter
   is not scalar.
 * The ``tolerance`` is the error between the final value and the objective function for which the optimization is allowed to stop.
 * The ``maxNumberOfEvolutions`` is the maximum number of iterations allowed in the optimization process if the convergence is not reached.
 * The ``minimizeOrMaximize`` flag is set to:
   
   * ``true`` if you want to minimize the objective function,
   * ``false`` if you want to maximize the objective function.

ObjectiveFunctionFromMinOrMaxDependentVariableSettings
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This class allows to use the minimum or maximum value in the propagation history of one of the dependent variables saved by the propagator as objective function.

 .. warning:: You need to have set the DependentVariableSaveSettings inside the PropagatorSettings in order to use this class!

**Constructor:**
 .. code-block:: cpp

    ObjectiveFunctionFromMinOrMaxDependentVariableSettings(
            boost::shared_ptr< propagators::SingleDependentVariableSaveSettings > dependentVariableSaveSettings,
            const double objectiveValue, 
            const int indexInVectorialVariable = 0, 
            const bool minimumOrMaximum = true,
            const double tolerance = 1e-3, 
            const unsigned int maxNumberOfEvolutions = 100 )

    ObjectiveFunctionFromMinOrMaxDependentVariableSettings(
            boost::shared_ptr< propagators::SingleDependentVariableSaveSettings > dependentVariableSaveSettings,
            const int indexInVectorialVariable = 0, 
            const bool minimizeOrMaximize = true, 
            const bool minimumOrMaximum = true,
            const unsigned int maxNumberOfEvolutions = 100 )

The first constructor allows to define an objective value towards which the objective function converges. The second constructor will
continue to optimize the configuration until the maximum number of evolutions is reached.

 * The ``dependentVariableSettings`` is the ``SingleDependentVariableSaveSettings`` object belonging to the saved variable whose final value 
   you want to use as objective function.
 * The ``objectiveValue`` is the value towards which the optimization converges.
 * The ``indexInVectorialVariable`` is the index pointing to the position that you need in the dependent variable in case the latter
   is not scalar.
 * The ``minimumOrMaximum`` flag is set to:
   
   * ``true`` if you want to use the **minimum** in the propagation history of the variable,
   * ``false`` if you want to use the **maximum** in the propagation history of the variable.

 * The ``tolerance`` is the error between the final value and the objective function for which the optimization is allowed to stop.
 * The ``maxNumberOfEvolutions`` is the maximum number of iterations allowed in the optimization process if the convergence is not reached.
 * The ``minimizeOrMaximize`` flag is set to:
   
   * ``true`` if you want to minimize the objective function,
   * ``false`` if you want to maximize the objective function.


UserDefinedObjectiveFunctionSettings
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This class allows to setup an objective function which is not available among the above. It could also be used to define
a weighted multi-objective optimization problem.

See the paragraph on :ref:`howToSetupUserDefinedObjectiveFunction` for information on how to use this class.

**Constructor:**
 .. code-block:: cpp

    UserDefinedObjectiveFunctionSettings( 
            boost::function< double () > userDefinedFunction,
            const double objectiveValue = 0.0, 
            const double tolerance = 1e-3,
            const unsigned int maxNumberOfEvolutions = 100 )

    UserDefinedObjectiveFunctionSettings( 
            boost::function< double () > userDefinedFunction,
            const bool minimizeOrMaximize,
            const unsigned int maxNumberOfEvolutions = 100 )

The first constructor allows to define an objective value towards which the objective function converges. The second constructor will
continue to optimize the configuration until the maximum number of evolutions is reached.

 * The ``maxNumberOfEvolutions`` is the maximum number of iterations allowed in the optimization process if the convergence is not reached.
 * The ``minimizeOrMaximize`` flag is set to:
   
   * ``true`` if you want to minimize the objective function,
   * ``false`` if you want to maximize the objective function.

.. _howToSetupUserDefinedObjectiveFunction:

How to set up a user defined objective function
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You might have noticed the use of ``boost::function< double () >`` and the fact that the ``UserDefinedObjectiveFunction`` class does not allow
an objective function with parameters. How do you setup an objective function without parameters?

First of all you will need a support class that takes all the parameters that you need, possibly as pointers, and a ``double`` method that does
not accept parameters. You can store in this object pointers to e.g. your dynamics simulator or propagator settings or whatever.

 .. code-block:: cpp

    struct ClassToDefineFunction{
        
        //Constructor. Make sure to input all the needed values
        //as pointers or boost::shared_ptr
        ClassToDefineFunction( argsType args ) :
            members_( args )
        {
            ...
        }
        
        ~ClassToDefineFunction( ){ }
   
        //This is the function that provides the calculated
        //value of the objective function.
        double myFunction( void )
        {
            double value;

            ...

            return value;
	}
     
        memberType members_

     };


Afterwards, when creating your objective function from ``UserDefinedObjectiveFunction`` you will need to use the ``boost::bind()`` statement to bind the method to a
``boost::function< double () >`` object. All you need to know is the following sintax, which is the only sintax that you are allowed to use in such a case:

 .. code-block:: cpp
   
     //Create object from class to define function and input all
     //the needed arguments
     boost::shared_ptr< ClassToDefineFunction > objectToDefineFunction = 
         boost::make_shared< ClassToDefineFunction >( args );

     //Create boost::function object from the above object
     boost::function< double () > myObjectiveFunction = 
             boost::bind( &ClassToDefineFunction::myFunction, 
                     objectToDefineFunction )

     //Create object to retrieve objective function from simulation
     boost::shared_ptr< UserDefinedObjectiveFunction > myObjectiveFunctionObject =
             boost::make_shared< UserDefinedObjectiveFunction >( myObjectiveFunction )


Et voilà. The user defined objective function is set.



