.. _walkthroughsUnperturbedEarthOrbitingSatellite:

Unperturbed Earth-orbiting Satellite
====================================
The example described on this page is that of Asterix, a single satellite that is orbiting the Earth. The code for this tutorial is given on Github, and is also located in your tudat bundle at::

   tudatBundle/tudatExampleApplications/satellitePropagatorExamples/SatellitePropagatorExamples/singleSatellitePropagator.cpp

For this example, we have the following problem statement:

   *Given the position and velocity of the Asterix satellite at a certain point in time with respect to the Earth, what will its position and velocity be after a Julian day has passed?*

There are a number of assumptions we can make in answering this question, which is often the case when translating a problem statement into a simulation. In our first tutorial, we will approximate the problem by the following:

    - The (initial) position and velocity of Asterix are given in the form of a Keplerian state vector.
    - The Earth is modeled as a point mass, and the spacecraft is modeled as massless.
    - All perturbations, such as solar radiation and third-body gravity, etc. are neglected.

.. note:: In a followup tutorial (see :ref:`walkthroughsPerturbedEarthOrbitingSatellite`), we will show how to extend the dynamical model of the spacecraft, slightly relaxing the last assumption.

The main focus of this tutorial is to provide insight into the use of the Tudat library for your application. Using this test case, we will show how you can use some of the elements of the Tudat library, namely:

    - The transformation from Keplerian elements (provided) to Cartesian elements (required by propagator).
    - The settings for the environment.
    - The settings for the acceleration models.
    - The settings for the state derivative model (propagator).
    - The settings for the numerical integrator.
    - The use of 'create' functions to parse these settings.

In particular, it will show you how these elements can work together to provide you with a powerful satellite propagator tool. It is not the intention to give you a C++ tutorial. Only the parts of the source code specific to Tudat will be discussed on this page. Detailed documentation on the environment models, acceleration models and propagators/integrators is available, but we recommend working through the tutorials here before delving into the details of all options that are available to you.

In Tudat, the environment, acceleration models and propagator are not created directly by the user. Instead, the user defines a number of 'Settings' objects that describe how the actual models are to be put together. These Settings objects are then passed to 'create' functions to build up the objects/data used in the simulations, and link them all together. In linking these objects together, information is automatically updated, and where necessary transformed and processed in the manner specified by the user in the Settings objects.

Set Up the Environment
~~~~~~~~~~~~~~~~~~~~~~
First we discuss the setup of the environment, which is stored in a list of :class:`Body` objects, each of which represents either a natural celestial body, or a spacecraft (orbiter, entry vehicle, etc.). The entire environment (gravity fields, atmospheres, ephemerides, rotation models, etc.) are then defined as members of the corresponding :class:`Body` object. For this example, the only environment models we need are:

    - Earth gravity field (point-mass).
    - Earth ephemeris model.

The block of code used to generate the environment that we require is:

.. code-block:: cpp

   // Create body objects.
   std::vector< std::string > bodiesToCreate;
   bodiesToCreate.push_back( "Earth" );
   std::map< std::string, std::shared_ptr< BodySettings > > bodySettings =
           getDefaultBodySettings( bodiesToCreate );    
   bodySettings[ "Earth" ]->ephemerisSettings = std::make_shared< ConstantEphemerisSettings >(
               Eigen::Vector6d::Zero( ) );

   // Create Earth object
   NamedBodyMap bodyMap = createBodies( bodySettings );

Creating an environment starts by creating the settings for all required environment models, which are stored in a list of :class:`BodySettings` objects. Each :class:`BodySettings` object has a list of objects for environment model settings, which are empty upon creation of the :class:`BodySettings`, and may be set by the user to their desired specifications. To save you the trouble of having to define all settings manually, we provide default options for the environment, many of which will often be sufficient for a given simulation. A list of all possible properties of the :class:`BodySettings`, as well as the default settings, can be found :ref:`here <tudatFeaturesCreatingTheEnvironment>`. In our example, we only need environment models for the Earth, so we start by:

.. code-block:: cpp

     std::vector< std::string > bodiesToCreate;
     bodiesToCreate.push_back( "Earth" );
     std::map< std::string, std::shared_ptr< BodySettings > > bodySettings =
             getDefaultBodySettings( bodiesToCreate );    

With the settings that are now stored in the :literal:`bodySettings` map, we could generate our environment and move on to the next step of the simulation. However, both to showcase one of the options of setting up the environment, and for numerical accuracy, we make one modification to the default ephemeris settings, defining the Earth to be fixed at the center of the Solar system:

.. code-block:: cpp

   bodySettings[ "Earth" ]->ephemerisSettings = std::make_shared< ConstantEphemerisSettings >(
           Eigen::Vector6d::Zero( ) );

.. warning:: By using this code, we 'cheat' a little bit, since we put the Earth in the main inertial frame for Tudat: the Solar System Barycenter (something that is done by default by the :class:`ConstantEphemerisSettings`, since the second argument is not specified). This approach should **only** be used when considering no third-body perturbations. Note that we could have used any ephemeris setting for the Earth (since we propagate our satellite w.r.t. the Earth origin) without changing anything in our dynamical model, but at a slight loss of numerical precision.

Now, we have created the settings we need for the environment, and we can move on to creating the environment models themselves, which are stored in a set of :class:`Body` objects. These :class:`Body` objects are stored in a :class:`NamedBodyMap`, which is a :literal:`typedef` (shorthand name) for:

.. code-block:: cpp

   std::unordered_map< std::string, std::shared_ptr< simulation_setup::Body > >

This :literal:`unordered_map` may be accessed as a regular map. The :literal:`std::string` keys represent the names of the bodies in the list, and the value :literal:`std::shared_ptr` the corresponding :class:`Body` object containing all the environment models. With the settings of our ephemeris and gravity field, we now create the :literal:`bodyMap` by:

.. code-block:: cpp

   NamedBodyMap bodyMap = createBodies( bodySettings );

Create the Vehicle
~~~~~~~~~~~~~~~~~~
Our environment is now missing only one aspect: the spacecraft. Our spacecraft (called Asterix) requires no specific properties, it merely needs to exist (its initial state is defined in a subsequent part of the example). Therefore, we can simply add the Asterix satellite to our environment by creating a new empty :class:`Body` object:

.. code-block:: cpp

   bodyMap[ "Asterix" ] = std::make_shared< simulation_setup::Body >( );

Although not required in this simulation, it is good practice to call the following function following the complete setup of the :literal:`bodyMap`:

.. code-block:: cpp

   setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );

Calling this function will allow hierarchical ephemerides to be properly used in the simulation (i.e. orbiter ephemeris w.r.t. Moon, Moon w.r.t. Earth, Earth w.r.t. Sun, Sun w.r.t. barycenter).

Set Up the Acceleration  odels
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
To define the settings of the propagation of the orbit, we start by defining the required acceleration models. The block of code that performs the required operations is:

.. code-block:: cpp

   // Define propagator settings variables.
   SelectedAccelerationMap accelerationMap;
   std::vector< std::string > bodiesToPropagate;
   std::vector< std::string > centralBodies;

   bodiesToPropagate.push_back( "Asterix" );
   centralBodies.push_back( "Earth" );

   // Define propagation settings.
   std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfAsterix;
   accelerationsOfAsterix[ "Earth" ].push_back( std::make_shared< AccelerationSettings >(
                                                    basic_astrodynamics::central_gravity ) );
   accelerationMap[ "Asterix" ] = accelerationsOfAsterix;

   // Create acceleration models and propagation settings.
   basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
               bodyMap, accelerationMap, bodiesToPropagate, centralBodies );

As with the environment models, there is no need to manually create the models. The user must only define the properties of the acceleration models that are desired, which are:

    - List of bodies that are to be numerically propagated.
    - Origin of reference frame in which they are to be propagated (may be different for each body). Note that propagation in Tudat is always done in a non-rotating reference frame, only the origin of the frames can be varied.
    - A list of settings of the accelerations model(s) acting on the bodies in the simulation.

These properties are to be defined in the following variables:

.. code-block:: cpp

   // Define propagator settings variables.
   SelectedAccelerationMap accelerationMap;
   std::vector< std::string > bodiesToPropagate;
   std::vector< std::string > centralBodies;

The list of propagated bodies and the reference frame origins (central bodies) are simply lists of strings, which for our case are:

.. code-block:: cpp

   bodiesToPropagate.push_back( "Asterix" );
   centralBodies.push_back( "Earth" );

The settings for the accelerations require some more structure, though, and are stored in a :class:`SelectedAccelerationMap`. This is :literal:`typedef` for:

.. code-block:: cpp

   std::map< std::string, std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > >

This is a double map (with twice a string as a key). The two levels correspond to the names of bodies undergoing an acceleration (first key), and those for bodies exerting an acceleration (second key). This allows any number of bodies to be propagated, undergoing any number (and type) of accelerations. Mutual acceleration between bodies being propagated, as is the case for Solar system dynamics for instance, is automatically handled by the code and requires no specific consideration.

In our example, we have only a single point-mass acceleration due to Earth, acting on Asterix. We define the settings for the acceleration as follows:

.. code-block:: cpp

   std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfAsterix;
   accelerationsOfAsterix[ "Earth" ].push_back( std::make_shared< AccelerationSettings >(
                                                    basic_astrodynamics::central_gravity ) );
   accelerationMap[ "Asterix" ] = accelerationsOfAsterix;

A single acceleration, of type :literal:`central_gravity` to be exerted by body :literal:`"Earth"` on body :literal:`"Asterix"` is now defined. The list of the actual acceleration models is now created by:

.. code-block:: cpp

   basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
               bodyMap, accelerationMap, bodiesToPropagate, centralBodies );

which automatically links together all required objects and functions.

Set Up the Propagation Settings
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Now that we have both our environment models and our acceleration model, we can create the full settings for the propagation. These settings are stored in a :class:`PropagatorSettings` object. For this example, we will only consider the propagation of translational dynamics, which is stored in the derived class :class:`TranslationalStatePropagatorSettings`. The settings for the propagator are the following:

    - The acceleration models.
    - The list of bodies that are to be propagated.
    - The origins w.r.t. which these bodies are to be propagated.
    - The initial Cartesian state that is to be used.
    - Termination conditions for the propagation (here, a fixed final time).

The above settings are provided in the following block of code:

.. code-block:: cpp

   // Set Keplerian elements for Asterix.
   Vector6d asterixInitialStateInKeplerianElements;
   asterixInitialStateInKeplerianElements( semiMajorAxisIndex ) = 7500.0E3;
   asterixInitialStateInKeplerianElements( eccentricityIndex ) = 0.1;
   asterixInitialStateInKeplerianElements( inclinationIndex ) = convertDegreesToRadians( 85.3 );
   asterixInitialStateInKeplerianElements( argumentOfPeriapsisIndex ) =
           convertDegreesToRadians( 235.7 );
   asterixInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex ) =
           convertDegreesToRadians( 23.4 );
   asterixInitialStateInKeplerianElements( trueAnomalyIndex ) = convertDegreesToRadians( 139.87 );
    
   // Convert Asterix state from Keplerian elements to Cartesian elements.
   double earthGravitationalParameter = bodyMap.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( );
   Eigen::VectorXd systemInitialState = convertKeplerianToCartesianElements(
               asterixInitialStateInKeplerianElements,
               earthGravitationalParameter );

   // Create propagator settings.
   std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
           std::make_shared< TranslationalStatePropagatorSettings< double > >
           ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState, simulationEndEpoch );

If the body that is being propagated has a pre-existing ephemeris, the initial state may be retrieved automatically. In this example, however, we manually define our initial state from the Keplerian state:

.. code-block:: cpp

   // Set Keplerian elements for Asterix.
   Vector6d asterixInitialStateInKeplerianElements;
   asterixInitialStateInKeplerianElements( semiMajorAxisIndex ) = 7500.0E3;
   asterixInitialStateInKeplerianElements( eccentricityIndex ) = 0.1;
   asterixInitialStateInKeplerianElements( inclinationIndex ) = convertDegreesToRadians( 85.3 );
   asterixInitialStateInKeplerianElements( argumentOfPeriapsisIndex )
           = convertDegreesToRadians( 235.7 );
   asterixInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex )
           = convertDegreesToRadians( 23.4 );
   asterixInitialStateInKeplerianElements( trueAnomalyIndex ) = convertDegreesToRadians( 139.87 );

   // Convert Asterix state from Keplerian elements to Cartesian elements.
   double earthGravitationalParameter = bodyMap.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( );
   Eigen::VectorXd systemInitialState = convertKeplerianToCartesianElements(
               asterixInitialStateInKeplerianElements,
               earthGravitationalParameter );

Note that we use the Earth gravity field inside the :literal:`bodyMap` for the conversion from Keplerian to Cartesian coordinates. Now, we can create our propagator settings with:

.. code-block:: cpp

   // Create propagator settings.
   std::shared_ptr< TranslationalStatePropagatorSettings< > > propagatorSettings =
       std::make_shared< TranslationalStatePropagatorSettings< > >
           ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState, simulationEndEpoch);

where we have passed exactly the five aspects listed above as input to the :class:`TranslationalStatePropagatorSettings`. If you have a look at the code for the :class:`TranslationalStatePropagatorSettings`, you will notice that there are multiple constructors for the class, each with a number of additional input arguments (for which we use the default values). These more advanced options are discussed in following tutorials.

A final piece of information needed to propagate the orbit is the settings object for the numerical integration. We use a Runge-Kutta 4 integrator, with a 10 second time step, starting the numerical integration at :math:`t = 0`:

.. code-block:: cpp

   // Create numerical integrator settings.
   double simulationStartEpoch = 0.0;
   const double fixedStepSize = 10.0;
   std::shared_ptr< IntegratorSettings< > > integratorSettings =
       std::make_shared< IntegratorSettings< > >
           ( rungeKutta4, simulationStartEpoch, fixedStepSize );

Perform the Orbit Propagation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Now, we have defined all the information needed to propagate the orbit of our satellite, which are stored in the :literal:`bodyMap` (environment), :literal:`propagatorSettings` (settings for the full state derivative model) and :literal:`integratorSettings` (settings on how to obtain the numerical solution). The propagation is done by an object of type (derived from) :class:`DynamicsSimulator`:

.. code-block:: cpp

   SingleArcDynamicsSimulator< > dynamicsSimulator(
           bodyMap, integratorSettings, propagatorSettings );

Upon creating this class, the numerical propagation is performed, and the output is stored in the class. Various options exist for parsing the output of the numerical propagation, which will be discussed in the next tutorials. The numerical solution of the orbit can be retrieved as follows:

.. code-block:: cpp

   std::map< double, Eigen::VectorXd > integrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );

where the retrieved result is a :literal:`std::map` where the key :literal:`double` is the time at each step in the integration, and the value, :literal:`Eigen::VectorXd`, is the corresponding Cartesian state of Asterix w.r.t. the Earth, in the ECLIPJ2000 reference frame. To analyze/plot your numerical results further using e.g. MATLAB, you can print the output to a text file as follows:

.. code-block:: cpp

   // Write satellite propagation history to file.
   input_output::writeDataMapToTextFile( integrationResult,
                                         "singleSatellitePropagationHistory.dat",
                                         tudat_applications::getOutputPath( ),
                                         "",
                                         std::numeric_limits< double >::digits10,
                                         std::numeric_limits< double >::digits10,
                                         "," );

For more details about the input and output basics go to: :ref:`tudatFeaturesInputOutput`. 

Results
~~~~~~~
The output of the program should look similar to the output below::

   Starting ../tudatBundle/tudatExampleApplications/satellitePropagatorExamples/bin/applications/application_SingleSatellitePropagator... 
   Single Earth-Orbiting Satellite Example.
   The initial position vector of Asterix is [km]:
   7037.48
   3238.06
   2150.72
   The initial velocity vector of Asterix is [km/s]:
     -1.46566
   -0.0409584
       6.6228
   After 86400 seconds, the position vector of Asterix is [km]:
   -4560.45
   -1438.32
    5973.99
   And the velocity vector of Asterix is [km/s]:
   -4.55021
   -2.41254
   -4.95063
   ../tudatBundle/tudatExampleApplications/satellitePropagatorExamples/bin/applications/application_SingleSatellitePropagator exited with code 0





