.. _addingNewStateDerivativeModel:

Adding a New State Derivative Model
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Another type of model which allows the user to add extra features to it, is the state derivative model. These models are used by the dynamics simulator to determine how the state equations are solved. Currently there are several available, e.g. cowell state derivative model, encke state derivate model, and more, but if a special model is needed for a certain application, the user can add this to tudat by following this guide.

As before, this guide will start by looking at the framework of how the state derivative model is implemented. First, the state derivative model is chosen by the user from an enum called :literal:`TranslationalPropagatorType`, located in :literal:`nBodyStateDerivativeModel.h`. The model picked from this enum is then used as an input into the constructor of the :literal:`TranslationalPropagatorSettings`, where it is assigned to a specific member variable called: :literal:`propagator_`. The :literal:`TranslationalPropagatorSettings` is used as an input to the :literal:`SingleArcDynamicsSimulator`, where it will be used further. An example of this is shown below:

.. code-block:: cpp


     // Create propagation settings.
     std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
             std::make_shared< TranslationalStatePropagatorSettings< double > >
             ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState,
               terminationSettings, cowell, dependentVariablesToSave );

    // Create simulation object and propagate dynamics.
    SingleArcDynamicsSimulator< > dynamicsSimulator(
             bodyMap, integratorSettings, propagatorSettings );


Inside :literal:`SingleArcDynamicsSimulator`, the propagator settings are passed to various functions, but the state derivative model is only needed in :literal:`createStateDerivativeModels`. This function returns a vector of :literal:`SingleStateTypeDerivative`. If a hybrid state derivative model is used, this function will fill the vector with state derivative models. If not, this vector will only contain one state derivative model, created by the function (equally) called :literal:`createStateDerivativeModels`. This function checks what kind of state is propagated (translational, rotational, etc.) and creates the specific state derivative model. For example, for the translational state it will call the :literal:`createTranslationalStateDerivativeModel` function. This function looks as follows:

.. code-block:: cpp


     template< typename StateScalarType = double, typename TimeType = double >
     std::shared_ptr< SingleStateTypeDerivative< StateScalarType, TimeType > >
     createTranslationalStateDerivativeModel(
            const std::shared_ptr< TranslationalStatePropagatorSettings< StateScalarType > >
            translationPropagatorSettings,
            const simulation_setup::NamedBodyMap& bodyMap,
            const TimeType propagationStartTime )
    {

        // Create object for frame origin transformations.
        std::shared_ptr< CentralBodyData< StateScalarType, TimeType > > centralBodyData =
                createCentralBodyData< StateScalarType, TimeType >(
                    translationPropagatorSettings->centralBodies_,
                    translationPropagatorSettings->bodiesToIntegrate_,
                    bodyMap );

        std::shared_ptr< SingleStateTypeDerivative< StateScalarType, TimeType > > stateDerivativeModel;

        // Check propagator type and create corresponding state derivative object.
        switch( translationPropagatorSettings->propagator_ )
        {
        case cowell:
        {
            stateDerivativeModel = std::make_shared<
                    NBodyCowellStateDerivative< StateScalarType, TimeType > >
                    ( translationPropagatorSettings->getAccelerationsMap( ), centralBodyData,
                      translationPropagatorSettings->bodiesToIntegrate_ );
            break;
        }
        ...
        ...
        ...
        }
        default:
            throw std::runtime_error(
                        "Error, did not recognize translational state propagation type: " +
                        std::to_string( translationPropagatorSettings->propagator_ ) );
        }
        return stateDerivativeModel;
    }

This function is responsible for making the state derivate model based on the :literal:`TranslationalPropagatorType` variable decalred by the user, using a switch statement.

When a new state derivative model is implemented in tudat, two additions need to be made to the framework (besides making the model itself). The first part is to add the name of the model to the enum: :literal:`TranslationalPropagatorType`. This will allow the user to pick this model and use it as an input to the translational propagator settings. The second step is to add this model to the switch statement in :literal:`createTranslationalStateDerivativeModel`. A new case should be made for the name added to the enum before, and in it, the variable :literal:`stateDerivativeModel` should be assigned to the new model.

When building the new model, it is advised to use a state derivative model that is already available as an example for a place to start. The classes which contain these models are derived from the base class :literal:`NBodyStateDerivative`, and should contain a template for the :literal:`TimeType` that will be used. The :literal:`NBodyStateDerivative` class is again derived from another class called the :literal:`SingleStateTypeDerivative`. This class contains several pure virtual functions, which all should be added to the new model class in order for the new model to work. The specific names and input parameters of these functions can be found in the :literal:`SingleStateTypeDerivative` class. Once this is done, and the new model is implemented in the state derivative model framework, the new model should be available for the user.

.. note:: Don't forget to put the include statement in :literal:`createStateDerivativeModel.h` if the new class is made in a seperate file.

