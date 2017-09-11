.. _tudatFeaturesFrameworkAccelerations:

Implementation Framework of the Acceleration Settings
=====================================================

Acceleration settings
~~~~~~~~~~~~~~~~~~~~~
The settings for accelerations are defined and stored by the user in a :class:`SelectedAccelerationMap`, which is :literal:`typedef` for:

.. code-block:: cpp

        std::map< std::string, std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > >

This is a double map (with twice a string as a key). The two levels correspond to the names of bodies undergoing an acceleration (first key) , and those for bodies exerting an acceleration (second key). This allows any number of bodies to be propagated, undergoing any number (and type) of accelerations from any number of bodies.
For a given environment, most acceleration models are completely defined by:

    - Type of acceleration model (a list is provided in the AvailableAcceleration enum in Tudat/Astrodynamics/BasicAstrodynamics/accelerationModelTypes.h.
    - Name of body undergoing acceleration
    - Name of body exerting acceleration

For instance, when using the following:

.. code-block:: cpp

    SelectedAccelerationMap accelerationSettings;
    accelerationSettings[ "Apollo" ][ "Earth" ].push_back( boost::make_shared< AccelerationSettings >( central_gravity ) );
    accelerationSettings[ "Apollo" ][ "Earth" ].push_back( boost::make_shared< AccelerationSettings >( aerodynamic ) );

We have defined a point-mass Earth gravity model and an aerodynamic acceleration due to Earth's atmosphere to be used. In this example, we have only defined the type of the acceleration, without the need for any additional information. All required variables used on the computations of the accelerations are uniquely defined in the Apollo and Earth entries of the body map (provided that the required environment models have been set).

However, as was the case for the settings of the environment models, certain types of accelerations require additional information. An important example is the spherical harmonic acceleration. We cannot replace :literal:`central_gravity` with :literal:`spherical_harmonic_gravity` in the above, as there is still an ambiguity in how the acceleration model is defined. In particular, we now also need the maximum degree and order of the gravity field that is to be used in addition to the three properties listed above. Consequently, we have created a dedicated :class:`AccelerationSettings` derived class for defining spherical harmonic acceleration settings. Updating the above example to use J_2, J_3 and J_4 (maximum degree = 4; maximum order = 0), we now have:

.. code-block:: cpp

    SelectedAccelerationMap accelerationSettings;
    accelerationSettings[ "Apollo" ][ "Earth" ].push_back( boost::make_shared< SphericalHarmonicAccelerationSettings >( 4, 0 ) );
    accelerationSettings[ "Apollo" ][ "Earth" ].push_back( boost::make_shared< AccelerationSettings >( aerodynamic ) );

A full list of the available acceleration models, as well as their required input and environment models, is given at the end of this page. One KEY point to understand is that, when creating an object of the :class:`AccelerationSettings` type (or its derived class), you must not provide any of the third body acceleration types (:literal:`third_body_central_gravity`, :literal:`third_body_spherical_harmonic_gravity`, :literal:`third_body_mutual_spherical_harmonic_gravity`) as input. If you wish to use a third-body gravity acceleration (typically from a point mass), simply provide :literal:`central_gravity` as input. Depending on the settings for your central bodies, the code will automatically create the corresponding acceleration object (central or third-body).

Having defined all the required settings for the accelerations in your :class:`SelectedAccelerationMap`, you can create the actual acceleration models by using the :literal:`createAccelerationModelsMap` function. This function requires four input parameters:

    - Full environment, as defined by a NamedBodyMap
    - Settings for the acceleration models, given by SelectedAccelerationMap
    - A list of bodies to numerically propagate
    - A list of central bodies (one for each numerically propagated body).

The list of central bodies defines the reference frame origins in which the bodies are propagated. The use of a hierarchical system is perfectly acceptable. For instance, one can propagate the Earth and Mars w.r.t. the Sun, the Sun w.r.t. the barycenter, the Moon w.r.t the Earth, etc. For this case, the central bodies and propagated bodies are defined as:

.. code-block:: cpp

    std::map< std::string, std::string > centralBodyMap;
    centralBodyMap[ "Moon" ] = "Earth";
    centralBodyMap[ "Earth" ] = "Sun";
    centralBodyMap[ "Mars" ] = "Sun";
    centralBodyMap[ "Sun" ] = "SSB";

There is no hardcoded limit to the number of permitted levels in the frame hierarchy, but it is not allowed to include circular dependencies, i.e. body A w.r.t. body B, body B w.r.t. body C and body C w.r.t. body A. More information of the acceleration models is discussed in the Propagators section. The following gives an example on how to create the acceleration model objects:

.. code-block:: cpp

    NamedBodyMap bodyMap;
    ....
    // Create environment here
    ....
    std::map< std::string, std::string > centralBodyMap;
    ....
    // Set central and propagated bodies here
    ....
    AccelerationMap accelerationModelMap = createAccelerationModelsMap( bodyMap, accelerationMap, centralBodyMap )

Mutual acceleration between bodies being propagated (i.e body A exerting acceleation on body B and vice versa), as is the case for solar system dynamics, is automatically handled by the :literal:`createAccelerationModelsMap` code and requires no specific consideration. Moreover, when creating a gravitational acceleration, the code checks whether it is a direct or a third-body gravitational acceleration and creates the acceleration models accordingly. Similarly, the code automatically checks which value of the gravitational parameter mu to use in such computations. For instance, when computing the gravitational acceleration due to the Sun acting on the Earth, :literal:`mu_Sun` is used when propagating w.r.t. the barycenter, whereas :literal:`mu_Sun + mu_Earth` is used when propagating w.r.t. the Sun.

For every acceleration, a model for the current state of the body exerting the acceleration must be available (the state of the body undergoing the acceleration is taken from the numerically propagated state). This means that, in the above example of the Apollo capsule entering Earth's atmosphere, we must include one of the following:

    - An ephemeris member for Earth.
    - Numerically integrate the Earth concurrently with our Apollo vehicle.

For this example, the second option is of course a bit 'non-standard'. However, for cases where entire planetary systems are propagated, such an approach is typically taken (for certain applications, the numerically propagated body must also have a particular ephemeris member object, see Propagators).

Available acceleration models
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
As stated above, the :literal:`createAccelerationModelsMap` function uses your environment and settings for the accelerations to automatically retrieve and put together all functions used to calculate the accelerations during each function evaluation of the numerical scheme. For reference, we provide a list of available acceleration models, below, including example of how to add settings for the model to the :class:`SelectedAccelerationMap`. In addition, we define the list of environment models required for their creation.

    **Point mass gravity:**
        No derived class of :class:`AccelerationSettings`, accessed by feeding :literal:`central_gravity` to the constructor. Added to :class:`SelectedAccelerationMap` as follows, for example of acceleration exerted on "Apollo" by "Earth":

        .. code-block:: cpp

            SelectedAccelerationMap accelerationSettings;
            accelerationSettings[ "Apollo" ][ "Earth" ].push_back( boost::make_shared< AccelerationSettings >( central_gravity ) );

        Requires the following environment models to be defined:

        - Gravity field for body exerting acceleration.
        - Current state of bodies undergoing and exerting acceleration, either from an Ephemeris model or from the numerical propagation.

    **Spherical harmonic gravity:**
        Accessed by means of the derived class :class:`SphericalHarmonicAccelerationSettings`. Added to :class:`SelectedAccelerationMap` as follows, for example of acceleration exerted on "Apollo" by "Earth":

        .. code-block:: cpp

            SelectedAccelerationMap accelerationSettings;
            int maximumDegree = 12;
            int maximumOrder = 12;
                accelerationSettings[ "Apollo" ][ "Earth" ].push_back( boost::make_shared< SphericalHarmonicAccelerationSettings >( maximumDegree, maximumOrder ) );

        where the gravity field will be expanded up to degree and order 12 in the acceleration model. Requires the following environment models to be defined:

        - Spherical harmonic gravity field for body exerting acceleration.
        - Rotation model from the inertial frame to the body-fixed frame.
        - Current state of bodies undergoing and exerting acceleration, either from an ephemeris model or from the numerical propagation.

    **Mutual spherical harmonic gravity:**
        Accessed by means of the derived class :class:`MutualSphericalHarmonicAccelerationSettings`. This model is typically only used for detailed propagation of planetary systems, and discussed in more detail here. It is added to :class:`SelectedAccelerationMap` as follows, for example of acceleration exerted on "Io" by "Jupiter":

        .. code-block:: cpp

            SelectedAccelerationMap accelerationSettings;
            int maximumDegreeOfIo = 12;
            int maximumOrderOfIo = 12;
            int maximumDegreeOfJupiter = 4;
            int maximumOrderOfJupiter = 4;
            accelerationSettings[ "Io" ][ "Jupiter" ].push_back( boost::make_shared< MutualSphericalHarmonicAccelerationSettings >( 
                maximumDegreeOfJupiter, maximumOrderOfJupiter, maximumDegreeOfIo, maximumOrderOfIo ) );

        where the gravity fields of Io and Jupiter will be expanded up to degree and order 12 and 4, respectively, in the acceleration model. Requires the following environment models to be defined:

        - Spherical harmonic gravity field for body exerting acceleration and body undergoing acceleration.
        - Rotation model from the inertial frame to the body-fixed frame and body undergoing acceleration.
        - Current state of bodies undergoing and exerting acceleration, either from an Ephemeris model or from the numerical propagation.

        For the case where a third-body mutual spherical harmonic acceleration (e.g. Ganymede on Io when propagating w.r.t. Jupiter), additional parameters have to be provided that denote the expansion degree/order of the central body, so:

        .. code-block:: cpp

            SelectedAccelerationMap accelerationSettings;
            int maximumDegreeOfIo = 12;
            int maximumOrderOfIo = 12;
            int maximumDegreeOfGanymede = 4;
            int maximumOrderOfGanymede = 4;
            int maximumDegreeOfJupiter = 4;
            int maximumOrderOfJupiter = 4;
            accelerationSettings[ "Io" ][ "Jupiter" ].push_back( boost::make_shared< MutualSphericalHarmonicAccelerationSettings >( 
                maximumDegreeOfJupiter, maximumOrderOfJupiter, maximumDegreeOfGanymede, maximumOrderOfGanymede, maximumDegreeOfIo, maximumOrderOfIo ) );

        where Jupiter now takes the role of central body, instead of body exerting the acceleration.

    **Aerodynamic acceleration:**
        No derived class of :class:`AccelerationSettings`, accessed by feeding :literal:`aerodynamic` to the constructor. Added to :class:`SelectedAccelerationMap` as follows, for example of acceleration exerted on "Apollo" by "Earth" (e.g. atmosphere model belonging to Earth):

        .. code-block:: cpp

            SelectedAccelerationMap accelerationSettings;
            accelerationSettings[ "Apollo" ][ "Earth" ].push_back( boost::make_shared< AccelerationSettings >( aerodynamic ) );

        Requires the following environment models to be defined:

        - Atmosphere model for body exerting acceleration.
        - Shape model for body exerting acceleration.
        - Aerodynamic coefficient interface for body undergoing acceleration. NOTE: In the case that the aerodynamic coefficients are defined as a function of the vehicle orientation (e.g. angle of attack and sideslip angle), these angles can be manually or automatically defined. We have dedicated a specific page to this here.
        - Mass model for body undergoing acceleration.
        - Current state of body undergoing and body with atmosphere.

        .. warning:: Defining settings for a vehicle's orientation, which may influence your aerodynamic force, is done after creating the acceleration models, as discused here.

    **Cannonball radiation pressure:**
        No derived class of :class:`AccelerationSettings`, accessed by feeding :literal:`cannon_ball_radiation_pressure` to the constructor. Added to :class:`SelectedAccelerationMap` as follows, for example of acceleration exerted on "Apollo" by "Sun":

        .. code-block:: cpp

            SelectedAccelerationMap accelerationSettings;
            accelerationSettings[ "Apollo" ][ "Sun" ].push_back( boost::make_shared< AccelerationSettings >( cannon_ball_radiation_pressure ) );

        Requires the following environment models to be defined:

        - Radiation pressure model for body undergoing acceleration (from source equal to body exerting acceleration)
        - Current state of body undergoing and body emitting radiation

    **Thrust acceleration:**
        Accessed by means of the derived class :class:`ThrustAccelerationSettings`, requiring:

    - Mass of body undergoing acceleration.
    - Settings for both the direction and magnitude of the thrust force. These models may in turn have additional environmental dependencies. The creation of thrust accelerations is discussed in more detail here.
    
    **Relativistic acceleration correction:**
        A first-order (in :math:`1/c^{2}`) correction to the acceleration due to the influence of relativity. It implements the model of Chapter 10, Section 3 of the IERS 2010 Conventions. These settings are defined by means of the derived class :class:`RelativisticAccelerationCorrectionSettings`, requiring:

    - Boolean whether to include the Schwarzschild correction term
    - Boolean whether to include the Lense-Thirring correction term
    - Boolean whether to include the de Sitter correction term
    - The name of the so-called 'primary body', for a planetary orbiter this should be set as the Sun (only relevant for de Sitter correction)
    - The angular momentum vector of the orbited body (only relevant for Lense-Thirring correction)
    
    **Empirical Acceleration**
    
       A constant/once-per-orbit acceleration, expressed in the RSW frame, for which the mangitude is determined empirically (typically during an orbit determination process). The acceleration components are defined according to Montenbruck and Gill (2000), with a total of 9 components: a constant, sine and cosine term (with true anomaly as argument) for each of the three independent directions of the RSW frame. The settings are defined by means of the derived class :class:`EmpiricalAccelerationSettings`, requiring:
       
    - Vector containing the constant terms of the accelerations in the R, S and W directions.
    - Vector containing the sine terms of the accelerations in the R, S and W directions.
    - Vector containing the cosine terms of the accelerations in the R, S and W directions.
    
Mass rate model setup
~~~~~~~~~~~~~~~~~~~~~
Although propagating a body's translational dynamics is the backbone of Tudat's simulations, it is also possible to propagate a vehicle's mass (either concurrently or separately). The manner in which the models that govern the 'mass dynamics', i.e. mass-rate models, are handled in the code is very similar to the acceleration models: a list of settings for the models is created by the user, which are then used to create the required objects. The list to be created by the user is:

.. code-block:: cpp

    std::map< std::string, std::vector< boost::shared_ptr< MassRateModelSettings > > > massRateModelSettings;

where the map key denotes the body of which the mass-rate is to be computed. At present, two mass-rate models are available, each with its own derived class of :class:`MassRateModelSettings`. These are:

    **Custom mass-rate:**
        Accessed by means of the derived class :class:`CustomMassRateModelSettings`. Using this class, the user must provide a :literal:`boost::function< double( const double ) > function`, i.e. a function returning a double, representing the mass-rate, and taking another double, representing time, as an input. The internal workings of this function are completely up to the user. If any help is required in setting up such a model please contact the Tudat support team.

    **From-thrust mass-rate:**
        Accessed by means of the derived class :class:`FromThrustMassModelSettings`. Using this mass-rate model, the change in vehicle mass due to the expulsion of propellant is taken into account when propagating a vehicle's dynamics. It retrieves the required data from a :class:`ThrustAcceleration` object, ensuring full consistency between the two. Two option are available when creating this type of mass-rate model:

        - Use all thrust forces acting on a single body, combined into a single mass-rate model. This will in most cases be the model of choice, as there is often no need to distinguish between thurst sources when computing the mass rate: only the total amount of propellant usage is relevant. This option is toggled by setting the :literal:`useAllThrustModels` input argument of the :class:`FromThrustMassModelSettings` constructor to true.
        - Use a single thrust model, defined by a string-identifier. When creating a thrust model, a :literal:`thrustOriginId` input is provided to the :class:`ThrustEngineSettings` settings constructor. Only in the :literal:`FromBodyThrustEngineSettings` derived type (see here for additional explanation) is this thrust origin id set to anything else than an empty string: it represents the engine name.

