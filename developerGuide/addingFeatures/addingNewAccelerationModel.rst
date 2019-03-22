.. _addingNewAccelerationModel:

Adding a New Acceleration Model
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
For some applications, the various acceleration models in Tudat might not be sufficient. This guide will explain how acceleration models are set-up within Tudat, and how these models can be modified, or how to add a completely new model.

Just as before, it is important to understand the framework of the acceleration models. This will be explained using an example, namely the aerodynamic acceleration model. The first time the acceleration models will be needed is in the user's own code. An example of how this can be done is shown below:

.. code-block:: cpp

    // Define propagator settings variables.
    SelectedAccelerationMap accelerationMap;
    std::vector< std::string > bodiesToPropagate;
    std::vector< std::string > centralBodies;

    // Define acceleration model settings.
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfWaverider;
    accelerationsOfWaverider[ "Mars" ].push_back( std::make_shared< AccelerationSettings >( aerodynamic ) );
    accelerationMap[  "Vehicle" ] = accelerationsOfWaverider;

    bodiesToPropagate.push_back( "Vehicle" );
    centralBodies.push_back( "Mars" );

    // Create acceleration models
    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
    bodyMap, accelerationMap, bodiesToPropagate, centralBodies );

The user will define a special map, with type :literal:`SelectedAccelerationMap` that works as follows: the first key defines the body exerting the acceleration, the second key the body undergoing the acceleration, and the value in the map is a vector of pointers to a type class called :literal:`AccelerationSettings`, which takes as input an :literal:`enum` called :literal:`AvailableAcceleration`, which lists all the available acceleration models. This map is then used as an input for the :literal:`createAccelerationModelsMap` function, which is a function that creates a special acceleration map that can be used as an input to the propagator settings.

There are two functions, located in :literal:`CreateAccelerationModels.cpp`, called :literal:`createAccelerationModelsMap`. The function that is called will first order the map, and then return the value of the other function, which takes as input the ordered map, made in the first :literal:`createAccelerationModelsMap`. The second function looks as follows:

.. code-block:: cpp

    //! Function to create a set of acceleration models from a map of bodies and acceleration model types.
    basic_astrodynamics::AccelerationMap createAccelerationModelsMap(
           const NamedBodyMap& bodyMap,
           const SelectedAccelerationMap& selectedAccelerationPerBody,
           const std::map< std::string, std::string >& centralBodies )
    {
        ...
        ...
        ...
            currentAcceleration = createAccelerationModel( bodyMap.at( bodyUndergoingAcceleration ),
                                                               bodyMap.at( bodyExertingAcceleration ),
                                                               accelerationsForBody.at( i ).second,
                                                               bodyUndergoingAcceleration,
                                                               bodyExertingAcceleration,
                                                               currentCentralBody,
                                                               currentCentralBodyName,
                                                               bodyMap );


            // Create acceleration model.
            mapOfAccelerationsForBody[ bodyExertingAcceleration ].push_back(
                         currentAcceleration );

        ...
        ...
        ...
            // Put acceleration models on current body in return map.
            accelerationModelMap[ bodyUndergoingAcceleration ] = mapOfAccelerationsForBody;
        }
        return accelerationModelMap;
    }

This function does several checks if the right models are present, and then calls the :literal:`createAccelerationModel( )` function. This function is where the acceleration models are called and looks as follows:

.. code-block:: cpp

    //! Function to create acceleration model object.
    std::shared_ptr< AccelerationModel< Eigen::Vector3d > > createAccelerationModel(
            const std::shared_ptr< Body > bodyUndergoingAcceleration,
            const std::shared_ptr< Body > bodyExertingAcceleration,
            const std::shared_ptr< AccelerationSettings > accelerationSettings,
            const std::string& nameOfBodyUndergoingAcceleration,
            const std::string& nameOfBodyExertingAcceleration,
            const std::shared_ptr< Body > centralBody,
            const std::string& nameOfCentralBody,
            const NamedBodyMap& bodyMap )
    {
        // Declare pointer to return object.
        std::shared_ptr< AccelerationModel< Eigen::Vector3d > > accelerationModelPointer;

        // Switch to call correct acceleration model type factory function.
        switch( accelerationSettings->accelerationType_ )
        {
        ...
        ...
        ...
        case aerodynamic:
            accelerationModelPointer = createAerodynamicAcceleratioModel(
                        bodyUndergoingAcceleration,
                        bodyExertingAcceleration,
                        nameOfBodyUndergoingAcceleration,
                        nameOfBodyExertingAcceleration );
            break;
        ...
        ...
        ...
        default:
            throw std::runtime_error(
                        std::string( "Error, acceleration model ") +
                        std::to_string( accelerationSettings->accelerationType_ ) +
                        " not recognized when making acceleration model of" +
                        nameOfBodyExertingAcceleration + " on " +
                       nameOfBodyUndergoingAcceleration );
            break;
        }
        return accelerationModelPointer;
    }

Ths function uses a switch statement on the :literal:`accelerationSettings->accelerationType_` (which was set by the user in the beginning), to determine which acceleration model needs to be called. The :literal:`createAerodynamicAcceleratioModel( )` then does several checks to determine if all the models that are used to calculate the aerodynamic force are present, after which it calls the constructor of the :literal:`AerodynamicAcceleration` class using the input values created in :literal:`createAerodynamicAcceleratioModel( )`. The structure of the :literal:`AerodynamicAcceleration` class, and how one can make a new acceleration model will be discussed in the next section.

When making a new acceleration model, there are two things that need to be added first. The first one is a new class containing functions to calculate the accelerations acting on the respective body. It is important that this class is derived from the base class: :literal:`basic_astrodynamics::AccelerationModel< dataType >`, where :literal:`dataType` is the type of the acceleration variable. Furthermore, there are two functions that need to be present in the new acceleration model class:


- :literal:`void updateMembers(const double currentTime = TUDAT_NAN)`. This function updates all the member variables to the current situation so they can be used by the other necessary function.
- :literal:`dataType getAcceleration( )`. This function uses the recently updated member variables to calculate the acceleration acting on the body.

For the case of the aerodynamic acceleration, these two functions look as follows:

.. code-block:: cpp


    Eigen::Vector3d getAcceleration( )
    {
        return computeAerodynamicAcceleration(
                    0.5 * currentDensity_ * currentAirspeed_ * currentAirspeed_,
                    currentReferenceArea_, currentForceCoefficients_, currentMass_ );
    }


    void updateMembers( const double currentTime = TUDAT_NAN )
    {
        if( !( this->currentTime_ == currentTime ) )
        {
            currentForceCoefficients_ = coefficientMultiplier_ * this->coefficientFunction_( );
            currentDensity_ = this->densityFunction_( );
            currentMass_ = this->massFunction_( );
            currentAirspeed_ = this->airSpeedFunction_( );
            currentReferenceArea_ = this->referenceAreaFunction_( );

            currentTime_ = currentTime;
        }
    }

where :literal:`computeAcceleration()` is a function that uses the well known lift, side, and drag acceleration equations to calculate the acceleration vector.

The second piece of code that needs to be added is the new create function inside :literal:`createAccelerationModel.cpp`. This function should return a pointer to the respective acceleration class, and takes as input the various models needed for the calculation of the acceleration. The main goal of this function is to do some checks if the necessary models are present, and either throw an error, or initialize the model there. When this is done, it uses all the available models as input to the acceleration model, and returns a pointer to this class.

Once the acceleration class is made, and the create function is in place, the acceleration model should be implemented into the acceleration framework. The first step for this is to add the a new model to the :literal:`AvailableAcceleration` enum. This name could be anything, as long as it relates to the actual acceleration model. Nothing has to be added to the :literal:`createAccelerationModelsMap`, however, the :literal:`createAccelerationModel` function has a switch statement that needs to be altered. This switch statement checks the :literal:`accelerationSettings` of the specific body to see which acceleration is acting on the body. When a new acceleration model is made, a new case should be made for the new model. This case should have the same name as the name added to the :literal:`AvailableAcceleration` enum. Inside this case, the variable :literal:`accelerationModelPointer` should be assigned to the create acceleration model function, which was made before. When this is done, the new acceleration model should be incorporated into the acceleration framework of tudat.


