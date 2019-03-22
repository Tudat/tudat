.. _addingNewEnvironmentModel:

Adding a New Environment Model
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Tudat contains various pre-built environment models that users are able to implement in their simulations. However, for some applications, it might be necesarry to alter an exisiting environment model, or add a completely new one. This section will describe how a new environment model can be added to Tudat.

First, it is important to understand how environment models are used in Tudat, this will give some insight in what needs to be added/changed if a new environment model is integrated into Tudat. For this tutorial, an existing environment model will be used as an example: the tabulated atmosphere model.

A user will define a list of bodies to be used in the simulation, an example is given here:

.. code-block:: cpp

    // Define simulation body settings.
    std::map< std::string, std::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( { "Mars" }, simulationStartEpoch - 10.0 * fixedStepSize,
                                    simulationEndEpoch + 10.0 * fixedStepSize );
    bodySettings[ "Mars" ]->ephemerisSettings = std::make_shared< simulation_setup::ConstantEphemerisSettings >(
                Eigen::Vector6d::Zero( ), "SSB", "J2000" );
    bodySettings[ "Mars" ]->rotationModelSettings->resetOriginalFrame( "J2000" );
    std::string atmosphereFile = ... ;
    bodySettings[ "Mars" ]->atmosphereSettings = std::make_shared< TabulatedAtmosphereSettings >( atmosphereFile );


The :literal:`BodySettings` class contains a member variable called :literal:`atmosphereSettings` of type :literal:`AtmosphereSettings`. As can be seen in the last line of the example code above, this member variable is set to a :literal:`TabulatedAtmosphereSettings` by the user. This is a valid statement (eventhough :literal:`atmosphereSettings` is not of type :literal:`TabulatedAtmosphereSettings`) because :literal:`TabulatedAtmosphereSettings` is derived from the base class :literal:`AtmosphereSettings`. These concepts are called polymorphism and inheritance, and should be understood by the reader before continuing.

When all the :literal:`bodySettings` are defined by the user, the following command will be executed:

.. code-block:: cpp

    simulation_setup::NamedBodyMap bodyMap = simulation_setup::createBodies( bodySettings );

The :literal:`simulation_setup::createBodies()` will set all the environment models up to be used later in the simulation, and stores them in the variable :literal:`bodyMap`. If the :literal:`simulation_setup::createBodies()` is looked at, the following can be seen:

.. code-block:: cpp

    // Create atmosphere model objects for each body (if required).
    for( unsigned int i = 0; i < orderedBodySettings.size( ); i++ )
    {
        if( orderedBodySettings.at( i ).second->atmosphereSettings != NULL )
        {
            bodyMap[ orderedBodySettings.at( i ).first ]->setAtmosphereModel(
                        createAtmosphereModel( orderedBodySettings.at( i ).second->atmosphereSettings,
                                               orderedBodySettings.at( i ).first ) );
        }
    }

Where :literal:`orderedBodySettings` is an ordered map of all the bodies. This code goes through all the bodies and checks if the bodies have an atmosphere model set or not. If not, the atmosphere model is set by the :literal:`setAtmosphereModel()` function, which takes as input a pointer to the atmosphereModel of the body. This pointer is returned by the function: :literal:`createAtmosphereModel()`, which creates the atmosphere model using the settings defined by the user in the variable :literal:`atmosphereSettings`. The :literal:`createAtmosphereModel()` looks as follows:

.. code-block:: cpp

        std::shared_ptr< aerodynamics::AtmosphereModel > createAtmosphereModel(
           const std::shared_ptr< AtmosphereSettings > atmosphereSettings,
           const std::string& body )
        {
            using namespace tudat::aerodynamics;

            // Declare return object.
            std::shared_ptr< AtmosphereModel > atmosphereModel;

            // Check which type of atmosphere model is to be created.
            switch( atmosphereSettings->getAtmosphereType( ) )
            {
            case exponential_atmosphere:
            {
                ...
                ...
                ...
            }
            case tabulated_atmosphere:
            {
                // Check whether settings for atmosphere are consistent with its type
                std::shared_ptr< TabulatedAtmosphereSettings > tabulatedAtmosphereSettings =
                        std::dynamic_pointer_cast< TabulatedAtmosphereSettings >( atmosphereSettings );
                if( tabulatedAtmosphereSettings == NULL )
                {
                    throw std::runtime_error(
                                "Error, expected tabulated atmosphere settings for body " + body );
                }
                else
                {
                    // Create and initialize tabulatedl atmosphere model.
                    atmosphereModel = std::make_shared< TabulatedAtmosphere >(
                                tabulatedAtmosphereSettings->getAtmosphereFile( ) );
                }
                break;
            }
            case nrlmsise00:
            {
                ...
                ...
                ...
            }
            default:
                throw std::runtime_error(
                            "Error, did not recognize atmosphere model settings type " +
                            std::to_string( atmosphereSettings->getAtmosphereType( ) ) );
            }

            return atmosphereModel;
        }

This function checks which atmosphere model is used by using a switch statement with the :literal:`atmosphereSettings->getAtmosphereType( )`. The three cases are defined using the enum :literal:`AtmosphereTypes` (located in :literal:`createAtmosphereModel.h`). In the case of the tabulated atmosphere, first a check is made if tabulated atmosphere setting are actually initialized. If they are, the :literal:`atmosphereModel` is set to a tabulated atmosphere using the settings from :literal:`atmosphereSettings`. It can be seen in the example code above that the variable :literal:`atmosphereModel` is of type :literal:`AtmosphereModel`, but it is set to a variable of type :literal:`TabulatedAtmosphere` in the tabulated atmosphere case. This is again valid due to the fact that the class :literal:`TabulatedAtmosphere` is derived from the base class :literal:`AtmosphereModel`.

After the atmosphere model is set in the body map, it can be used whenever a certain quantity, e.g. the density or temperature, is needed by another part of the simulation. How these quantities are calculated is defined in the :literal:`TabulatedAtmosphere` class.

Now that the creation of an environment model is understood, it can be discussed what should change in the Tudat code when something is added. There are two options when changing the code: the user can modify an existing environment model, or they can add a new model to an existing type of environment model. These type of modifications require different amount of changes made to Tudat and are thus explained seperately here:

- **Modify existing environment model:** this option is the easiest as it only requires changes in the specific environment type files. If a function is added to an environment model, it is important to also include this function in the base class file of that specific environment model, with the virtual statement. For example, take the :literal:`getDensity()` function. This function is put in the file :literal:`AtmosphereModel.h` as a pure virtual function (don't need a function definition), by putting a :literal:`=0` after the function definition and including the virtual statement. Now, in every derived class, this function should return something, depending on the implementation. If only something inside an already existing function needs to be changed, it (most of the time) shouldn't be changed in all the derived classes. If a variable is added to the constructor of the class, it is important that all the cases that this class is called should be changed in the entire Tudat code to prevent errors (use the :literal:`find usages` option in Qt to find them).

- **Add a new environment model:** this modification requires some extra changes to the framework of the environment model implementation. First, when the new model is made, make sure that it is derived from the base class, and that it contains some of the basic functions. Once the new model is made, a :literal:`Settings` class should be made in the same way as :literal:`TabulatedAtmosphereSettings`. This class should be added to the :literal:`create...` files, and should contain functions that store variables that are used to call the constructor of the corresponding class (again, make sure it is derived from the proper base class). Then, in the corresponding :literal:`.cpp` file, a new case for the new environment model should be added to the :literal:`create...` function, just as in the tabulated atmosphere example. Make sure that in the :literal:`getAtmosphereType( )` function, the name of the new model is also included. In the case statement, make sure to add checks, and throw runtime errors if they are violated. Another step that needs to be taken is to update the :literal:`createEnvironmentUpdater.cpp` file. This file includes several switch statements that need to have the new acceleration model in it. Use the existing code to determine how the new case should be made.

