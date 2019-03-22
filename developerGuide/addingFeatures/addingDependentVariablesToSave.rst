.. _addingDependentVariablesToSave:

Adding Dependent Variables to Save
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
As the execution of the simulator happens "behind closed doors" in the :literal:`DynamicsSimulator` function, it can be non-trivial to retrieve the history of dependent variables. Luckily, as explained in :ref:`tudatFeaturesPropagatorSettingsDependentVariables`, it is possible to create a list of dependent variables and use it as an input to the propagation settings to be able to retrieve a history of these dependent variables after the simulation is done. A list of the available dependent variables is given in: :ref:`tudatFeaturesPropagatorSettingsDependentVariables`. It is possible that another dependent variable is needed by the user, which is not given in this list. This guide will show how a new dependent variable can be added to the code.

To add a new dependent variable, four files will need to be changed, these are all located in the following directory:

        .../tudatBundle/tudat/Tudat/SimulationSetup/PropagationSetup/

For this guide, the :literal:`mach_number_dependent_variable` will be used as an example of how a dependent variable is implemented. Be aware that there are some slight differences between the implementation of the :literal:`mach_number_dependent_variable` and other variables, but this will be explained in the guide.

The first addition needs to be made in the :literal:`propagationOutputSettings.h` file. Here, an :literal:`enum` called :literal:`PropagationDependentVariables` is given which shows the names of all the available dependent variables that can be saved. The list looks as follows;

.. code-block:: cpp
   
        enum PropagationDependentVariables
        {
            mach_number_dependent_variable = 0,
            altitude_dependent_variable = 1,
            airspeed_dependent_variable = 2,
            ...
            ...
            ...
            body_fixed_groundspeed_based_velocity_variable = 31,
            keplerian_state_dependent_variable = 32,
            modified_equinocial_state_dependent_variable = 33
        };


The name of the new variable can thus be added to the end of the list, with a number assigned to it.

The next addition needs to be made in the :literal:`propagationOutputSettings.cpp` file. In this file, a method is present called: :literal:`getDependentVariableName`. This method contains a switch statement which contains cases for all the different dependent variables. For a new variable, a new case must be made with the name that was added to the :literal:`PropagationDependentVariables` list. Inside the case, a string called :literal:`variableName` must be assigned with the name of the to be added variable. In the case of the mach number, it looks as follows:

.. code-block:: cpp
   
        case mach_number_dependent_variable:
            variableName = "Mach number ";
            break;

It is important to remember that the case should be the same as the name given in :literal:`PropagationDependentVariables`, but the :literal:`variableName` can be chosen to the preference of the user.

The third addition needs to be made in :literal:`propagationOutput.cpp`. Here a method called :literal:`getDependentVariableSize` is located. In this function, another switch statement is made with the same cases as before. Again a new case should be made for the new variable, but now, inside the case, another variable called :literal:`variableSize` should be changed to the size of the added variable. Thus this could be 1 for a :literal:`double`, or 3 for an :literal:`Eigen::Vector3d`. In the case of the mach number, it looks as follows:

.. code-block:: cpp
   
        case mach_number_dependent_variable:
            variableSize = 1;
            break;

The final change should then be made to the :literal:`propagationOutput.h` file. In a method called :literal:`getDoubleDependentVariableFunction`, the actual calculation of the variable is done. Again, a switch stamentent in the same manner as before is made, with in every case the calculation of the variable. This method gets the :literal:`bodyMap` as input, thus methods available inside the :literal:`FlightConditions` could, for example, be used to calculate the dependent variable. For the new variable, a new case needs to be made in which a :literal:`std::function< double( )>` is returned. This function can take several variables as input and should return the dependent variable. The implementation for the mach number is given here:

.. code-block:: cpp
   
        case mach_number_dependent_variable:
        {
            if( bodyMap.at( bodyWithProperty )->getFlightConditions( ) == NULL )
            {
                std::string errorMessage = "Error, no flight conditions available when requesting Mach number output of " +
                        bodyWithProperty + "w.r.t." + secondaryBody;
                throw std::runtime_error( errorMessage );
            }

            std::function< double( const double, const double ) > functionToEvaluate =
                    std::bind( &aerodynamics::computeMachNumber, std::placeholders::_1, std::placeholders::_2 );

            // Retrieve functions for airspeed and speed of sound.
            std::function< double( ) > firstInput =
                    std::bind( &aerodynamics::FlightConditions::getCurrentAirspeed,
                                 bodyMap.at( bodyWithProperty )->getFlightConditions( ) );
            std::function< double( ) > secondInput =
                    std::bind( &aerodynamics::FlightConditions::getCurrentSpeedOfSound,
                                 bodyMap.at( bodyWithProperty )->getFlightConditions( ) );


            variableFunction = std::bind( &evaluateBivariateFunction< double, double >,
                                            functionToEvaluate, firstInput, secondInput );
            break;
      }

If the variable is a vector (or matrix) and not a double, the case should not be added to :literal:`getDoubleDependentVariableFunction`, but to: :literal:`getVectorDependentVariableFunction`

If this is all done, the dependent variable name can be added to the dependent variables save list.
