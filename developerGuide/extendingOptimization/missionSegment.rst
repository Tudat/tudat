.. _extendMissionSegment:

Mission Segment Settings
========================

What is it
~~~~~~~~~~

A **mission segment** in this framework is defined by a ``SingleArcDynamicsSimulator``,
a ``DecisionVariableSettings`` object which establishes what variables in the dynamics simulator
have to be tuned in order to optimize the mission and an ``ObjectiveFunctionSettings`` (or derived classes)
object used to define the function to optimize.

The class is ``MissionSegmentSettings``. See the file ``missionSegmentSettings.h`` to visualize the code.

Linkage aiding methods
~~~~~~~~~~~~~~~~~~~~~~

Most of the methods used to link together two consecutive mission segments can be found in the
``MissionSegmentSettings`` class. Their codes can be found in the file ``missionSegmentSettings.cpp``.

They all use :ref:`optimization-linearInterpolation` to retrieve the final states.

Their implementation is pretty straightforward. I will leave the developer to explore the code by their own.

Decision variable settings
~~~~~~~~~~~~~~~~~~~~~~~~~~

Most of the methods to modify decision variables are also found in ``missionSegmentSettings.cpp``.
Their implementation is also pretty straightforward. The developer must have a previous knowledge of Tudat and
:ref:`optimization-dynamicPointerCast`.


.. _getObjectiveFunctionValue:

The method getObjectiveFunctionValue()
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The retrieval of the value is performed in the class ``MissionSegmentSettings`` using the function 

    .. code-block:: cpp
        
        double getObjectiveFunctionValue( void )

    .. warning:: The objective function based on final integration values must always be interpolated to the final conditions imposed by the user in the termination settings.

Let us take into consideration the case where the objective function is of the type ``FinalSphericalOrbitalElementObjectiveFunctionSettings``. 
In this case the code is going to retrieve, after the propagation, the final spherical orbital element defined by the user in the objective function settings
in the state map of the desired body.
Here's the code:

    .. code-block:: cpp
        :caption: :class:`Tudat/Optimization/missionSegmentSettings.cpp`
        :name: missionSegmentSettings-cpp
        
        double MissionSegmentSettings::getObjectiveFunctionValue( void )
        {

            double finalTime = getFinalSimulationEpoch();

            double objectiveFunctionValue;

            ...

            else if( objectiveFunctionSettings_->objectiveVariable_ == final_spherical_orbital_element_objective_function )
            {
                boost::shared_ptr< FinalSphericalOrbitalElementObjectiveFunctionSettings > objectiveFunctionSettings =
                        boost::dynamic_pointer_cast< FinalSphericalOrbitalElementObjectiveFunctionSettings >(
                             objectiveFunctionSettings_);
                std::map< std::string, Eigen::Vector6d > finalStates = getFinalOrbitalStates( finalTime );
                Eigen::Vector6d sphericalState = orbital_element_conversions::convertCartesianToSphericalOrbitalState(
                        finalStates[objectiveFunctionSettings->associatedBody_] );
                objectiveFunctionValue = sphericalState[objectiveFunctionSettings->sphericalOrbitalElement_];
            }
        
            ...
        
        }

First of all notice that the method ``getFinalSimulationEpoch()`` has been used to retrieve the variable finalTime. This is a method of the ``MissionSegmentSettings`` class,
which calculates the final epoch according to the user defined termination settings. Please visit the relative documentation for more information.
The member ``objectiveFunctionSettings_`` in ``MissionSegmentSettings`` is stored as a ``boost::shared_ptr< ObjectiveFunctionSettings >`` object.
In order to retrieve it as a ``FinalSphericalOrbitalElementObjectiveFunctionSettings``, which is a derived class, a dynamic pointer cast of the
member is performed, after constatating that

    .. code-block:: cpp

        objectiveFunctionSettings_->objectiveVariable_ == final_spherical_orbital_element_objective_function


The rest of the code is pretty straightforward. The final orbital state is recovered with the ``MissionSegmentSettings``method
``getFinalOrbitalStates( double )`` which accepts the epoch to which the states are interpolated.

    .. code-block:: cpp

        std::map< std::string, Eigen::Vector6d > finalStates = getFinalOrbitalStates( finalTime );

conveniently the method maps the orbital states in a map with the name of the body as key.

All the methods shown above to retrieve the final states are based on linear interpolation. In some cases there is not a method directly
defined, in any case the final result must always be interpolated to the user defined final conditions. For example, for the case of
``ObjectiveFunctionFromFinalDependentVariableSettings``, the code is the following:

    .. code-block:: cpp
        :caption: :class:`Tudat/Optimization/missionSegmentSettings.cpp`
        :name: missionSegmentSettings-cpp


        double MissionSegmentSettings::getObjectiveFunctionValue( void )
        {
        
            double finalTime = getFinalSimulationEpoch();

            double objectiveFunctionValue;
            ...

            else if( objectiveFunctionSettings_->objectiveVariable_ == final_dependent_variable_objective_function )
            {
                boost::shared_ptr< ObjectiveFunctionFromFinalDependentVariableSettings > objectiveFunctionSettings =
                         boost::dynamic_pointer_cast< ObjectiveFunctionFromFinalDependentVariableSettings >(
                            objectiveFunctionSettings_ );
                int positionInVector = objectiveFunctionSettings->getPositionInDependentVariablesSaveVector(
                            dynamicsSimulator_ );
                std::map< double, Eigen::VectorXd > dependentVariableHistory
                        = dynamicsSimulator_->getDependentVariableHistory();
                std::map< double, Eigen::VectorXd >::reverse_iterator it = dependentVariableHistory.rbegin();

                double time1 = it->first;
                double value1 = it->second[positionInVector];

                it++;

                double time2 = it->first;
                double value2 = it->second[positionInVector];

                // If the last two epochs do not contain the final time in their boundaries
                // scout all the epochs
                while( !( time2 <= finalTime && time1 >= finalTime) )
                {
                    time1 = time2;
                    value1 = value2;

                    it++;

                    time2 = it->first;
                    value2 = it->second[positionInVector];

                }

                objectiveFunctionValue = value2 + ( value1 - value2 )*(finalTime - time2)/(time1 - time2);
            }
            
            ...
        }

The method ``ObjectiveFunctionFromFinalDependentVariableSettings::getPositionInDependentVariablesSaveVector()`` is used in this case
to retrieve the position of the dependent variable settings defined in the objective function object from the dynamics simulator.

The method ends with the lines:

    .. code-block:: cpp
        :caption: :class:`Tudat/Optimization/missionSegmentSettings.cpp`
        :name: missionSegmentSettings-cpp

        ...

        if( objectiveFunctionSettings_->objectiveValueIsSet_ == true )
            return fabs( objectiveFunctionSettings_->objectiveValue_ - objectiveFunctionValue );
        else if( objectiveFunctionSettings_->minimizeOrMaximize_ )
            return objectiveFunctionValue;
        else
            return -objectiveFunctionValue;

        ...

The first if-case is for an objective function settings defined with the first constructor, in which case the return value is the
absolute difference of the retrieved value and the objective value. The other two cases are for an objective function settings defined 
with the second constructor: for the second if-case (minimize) simply the retrieved value is returned, otherwise (maximize) the
negative retrieved values is retrieved. That's because PaGMO is set to minimize the function.
