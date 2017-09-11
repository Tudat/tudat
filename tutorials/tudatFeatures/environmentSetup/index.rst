.. _tudatFeaturesEnvironmentIndex:

Environment Set-up
==================
Models for the physical environment are one of the cornerstones of a numerical astrodynamics toolbox. Here, we define the environment in the broadest sense, including all physical properties of the solar system, such as atmospheres and gravity fields, but also any models for the orbits and rotations of these bodies.

On this page, we will give an overview of how the environment is represented in Tudat, which models have been implemented, and how to create the environment that is tailored to your needs.

Setting up the Environment
~~~~~~~~~~~~~~~~~~~~~~~~~~
In Tudat, the physical environment is defined by a list of :class:`Body` objects, each of which represents either a celestial body, or a manmade vehicle. Consequently, all properties that are required for computing e.g. accelerations are stored in :class:`Body` objects. Typically the entire environment is stored in a named list of :class:`Body` object, the standard typedef for which is the :class:`NamedBodyMap`. (It as an unordered map of shared pointers to body objects, see this wiki page for a discussion of shared pointers; don't worry if you're not sure what a shared pointer or unordered map is at this point.).

Manually creating the environment
*********************************
The following shows how to manually declare a :class:`NamedBodyMap`, and then create entries in this body map for a number of bodies:

.. code-block:: cpp

    NamedBodyMap bodyMap;
    bodyMap[ "Earth" ] = boost::make_shared< Body >( );
    bodyMap[ "Moon" ] = boost::make_shared< Body >( );
    bodyMap[ "Sun" ] = boost::make_shared< Body >( );
    bodyMap[ "Apollo" ] = boost::make_shared< Body >( );

This creates four body objects (representing three celestial bodies and one vehicle; Tudat does not distinguish between the two). However, these bodies do not yet have any physical properties, the :literal:`bodyMap` created above now only indicates the existence of these four bodies.

To actually define the physical properties of the environment, a :class:`Body` object may be endowed with any of a number of properties. In particular, the following properties may be set. A more extensive list of possible model types is given at the end of this tutorial page:

    - **Ephemeris:** defines the state of the body as a function of time (Dynamical Barycentric Time seconds since J2000 is default).
    - **Gravity field:** defines the gravity field of the body, in terms of its gravitational potential and associated quantities.
    - **Time-variations of the gravity field:** defines models for the time-dependency of this gravity field.
    - **Atmosphere model:** defines the atmospheric properties (density, temperature, etc.) as a function of relative position and time
    - **Shape model:** defines the shape of a body, from which for instance the altitude of another body can be computed
    - **Rotation model:** defines the instantaneous rotation matrix (and its time derivative) of the body-fixed frame, w.r.t. some inertial frame.
    - **Aerodynamic coefficient interface:** defines the aerodynamic properties of the body, such as its aerodynamic coeficients as a function of some set of independent variables.
    - **Radiation pressure interface:** defines the radiation pressure properties of the body.
    - **Mass model:** defines the mass of a body (possibly as a function of time). This separate function is typically used for vehicles only. For celestial bodies, the mass is typically derived from the gravity field member (if applicable).
    - **Vehicle system models:** This is a container object that stores properties of systems and physical properties of a vehicle. The options in this container are presently limited to propulsion systems and some physical characteristics related to entry heating.

These properties can be set manually or default settings can be used. For instance, to manually create and set an ephemeris (from Spice w.r.t. the barycenter) and gravity field (point-mass only) object in the 'Earth' entry of the body map, the following can be used:

.. code-block:: cpp

    bodyMap[ "Earth" ]->setEphemeris( boost::make_shared< SpiceEphemeris >( "Earth", "SSB", false, false, true, "J2000" ) ); 
    bodyMap[ "Earth" ]->setGravityFieldModel( boost::make_shared< GravityFieldModel >( 3.986004418E14 ) );  

This calls the constructors of the :class:`SpiceEphemeris` and :class:`GravityFieldModel` classes (see here and here for details), and assigns the objects that are constructed to the "Earth" entry of the ``bodyMap``.

Creating the environment from :class:`BodySettings`
***************************************************
Manually creating all objects defining the full environment is possible, but not recommended. In particular, various environment models are interdependent and these dependencies must be fully and consistently defined for the code to function properly. To this end, we provide a :class:`BodySettings` object, in which the general properties of each environment model can be set (see above for the list of the available types of environment models). We note that for :class:`Body` objects that represent vehicles, the manual creation is typically used, as the vehicle conditions may depend on the celestial bodies, but not vice versa.

In many cases, default properties of (celestial) bodies may be used by calling the :literal:`getDefaultBodySettings` function, so that the user does not need to define all required properties line-by-line. At present, the following default settings are used (none if not in list):

    - **Ephemeris:** Tabulated ephemeris created from Spice (valid in the interval that is specified by the input time-arguments to getDefaultBodySettings).
    - **Gravity field models:** Point mass gravity field models, with gravitational parameter from Spice (if available). Exceptions are the Earth and Moon, for which the EGM96 and GGLP spherical harmonic gravity fields are loaded, respectively.
    - **Rotation model:** For a given body (if available) the Spice rotation model, with ECLIPJ2000 as base frame, and for a body AAA frame IAU_AAA as target frame (the standard body-fixed frame for each body in Spice).
    - **Atmosphere model:** 1976 US Standard Atmosphere for Earth (using pregenerated tables). For other bodies, no default shape model is given.
    - **Shape model:** Spherical model with mean radius obtained from Spice (if avaiable).

The default settings for a body are loaded as follows:

.. code-block:: cpp

    std::vector< std::string > bodyNames;
    bodyNames.push_back( "Earth" );
    bodyNames.push_back( "Sun" );
    bodyNames.push_back( "Moon" );
    bodyNames.push_back( "Mars" );
    double initialEphemerisTime = 1.0E7;
    double finalEphemerisTime = 2.0E7;
    double buffer = 5000.0;
    std::map< std::string, boost::shared_ptr< BodySettings > > bodySettings =
                getDefaultBodySettings( bodyNames, initialEphemerisTime - buffer, finalEphemerisTime + buffer );

The reasons for passing the initial/final time (as well as the 'buffer') are discussed in more detail here:"http://tudat.tudelft.nl/projects/tudat/wiki/timeBuffer". As can be seen from the above, the settings for the environment are stored in a map of pointers to :class:`BodySettings` objects (with the key the name of the associated bodies). If you have a look at the definition of the :class:`BodySettings` class (in ``SimulationSetup/createBodies.h``), you will see that this type is simply a contained for a list of specific environment settings, which we discuss in more detail below. As a result, specifying settings for a given type of environment model requires the creation of an object of the correct type of class (derived class of :class:`EphemerisSettings` for defining an ephemeris; derived class of :class:`BodyShapeSettings` for defining a body shape etc.)

Often, one will wish to load the default settings, but make slight modifications or additions to it before creating the :class:`NamedBodyMap`. This can be achieved as follows for the example of a shape model: we want an oblate spheroid shape model instead of a spherical shape model for Earth.

.. code-block:: cpp

    bodySettings[ "Earth" ]->shapeModelSettings = boost::make_shared< OblateSphericalBodyShapeSettings >( 6378.0E3, 0.01 );

which changes the shape model settings of the Earth from the default spherical to the oblate spheroid. A list of available environment models, as well as the manner in which to provide settings for them, is provided at the end of this tutorial. The above appproach is identical for adding or modifying environment model settings (that is, it does not matter whether Earth already had ``shapeModelSettings`` or not). Once the settings for the environment model have been defined, the following creates the actual :class:`Body` objects and all associated environment models

.. code-block:: cpp

    NamedBodyMap bodyMap = createBodies( bodySettings );

It should be noted that default settings presently exist only for celestial bodies. The addition of objects to represent vehicles may be done either at the settings level (appending the ``bodySettings`` map) or at the body object level (appending the ``bodyMap``). Here, we give the example of directly appending the ``bodyMap``. For instance, creating an Apollo entry vehicle object, and adding a mass and aerodynamic properties is achieved as follows:

.. code-block:: cpp

    bodyMap[ "Apollo" ] = boost::make_shared< Body >( );
    bodyMap[ "Apollo" ]->setAerodynamicCoefficientInterface( getApolloCoefficientInterface( ) );
    bodyMap[ "Apollo" ]->setConstantBodyMass( 5.0E3 );

where the ``getApolloCoefficientInterface`` is a predefined function that generates an aerodynamic database from the Apollo capsule's shape. A final, but crucial step in creating the bodyMap is the following:

.. code-block:: cpp

    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );
 
This line of code allows the ephemerides and rotation models of the various bodies to be defined w.r.t. different origins (and even w.r.t. each other).

Available Settings for the Environment Models
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Here, we will provide a full list of the available properties of the :class:`BodySettings` object. Each type of environment model has one base class to define settings for the creation of the model). Often, a specific derived class is implemented for a specific environment model of a given class, in which any additional information that may be needed can be provided. For instance, when defining a gravity field model, one can simply use:

.. code-block:: cpp

    bodySettings[ "Earth" ]->gravityFieldSettings = boost::make_shared< GravityFieldSettings >( central_spice ); 

if you want to use a central gravity field with the gravitational parameter taken from Spice: no information is needed except the type of gravity field model that is created. On the other hand, if you want to use a spherical harmonic gravity field, you need to specify additional parameters yourself, which is done by using the specific derived class:

.. code-block:: cpp

    bodySettings[ "Earth" ]->gravityFieldSettings = boost::make_shared< SphericalHarmonicsGravityFieldSettings >( gravitationalParameter, referenceRadius, cosineCoefficients, sineCoefficients, associatedReferenceFrame ); 

To find out which input arguments must be provided to create a specific settings class, have a look at the documentation in the code (written above the code for the constructor of the settings class you are interested in). Below, we give examples of each type of environment model setting.

The full list of available environment model settings is:

    **Atmosphere model:** defined in the atmosphereSettings member (of type ``shared_ptr< AtmosphereSettings >``). Models currently available through the :class:`BodySettings` architecture are (with examples when defining settings for Earth):
    
        - **Exponential atmosphere:** Simple atmosphere model independent of time, latitude and longitude based on an exponentially decaying density profile with a constant temperature. Settings stored in :class:`ExponentialAtmosphereSettings` object.

        .. code-block:: cpp

            bodySettings[ "Earth" ]->atmosphereSettings = boost::make_shared< ExponentialAtmosphereSettings >( 7.2E3, 290.0, 1.225, 287.06 ); 

        for an exponential atmosphere with a scale height of 7200 m, a constant temperature of 290 K, a density at 0 m altitude of 1.225 kg/m^3 and a specific gas constant of 287.06 J/(kg K).

        - **Tabulated atmosphere:** Atmosphere model with properties (pressure, density, temperature) read in from a file. Current implementation is independent of time, latitude and longitude. Settings stored in TabulatedAtmosphereSettings object. 

        .. code-block:: cpp

            std::string atmosphereFile = ...
            bodySettings[ "Earth" ]->atmosphereSettings = boost::make_shared< TabulatedAtmosphereSettings >( atmosphereFile ); 

        which will read the atmospheric properties from the file atmosphereFile (with four columns altitude and associated presure, density and temperature).

        - **NRLMSISE00:** Atmosphere model using the NRLMSISE-00 atmosphere model. To use this model, the :literal:`USE_NRLMSISE` flag in your top-level :literal:`CMakeLists` must be set to true. No derived class of :class:`AtmosphereSettings` base class required, created by passing :literal:`nrlmsise00` as argument to base class constructor.

        .. code-block:: cpp

            bodySettings[ "Earth" ]->atmosphereSettings = boost::make_shared< AtmosphereSettings >( nrlmsise00 ); 

    **Ephemeris model:** defined in the ephemerisSettings (of type ``shared_ptr< EphemerisSettings >``). Models currently available through the :class:`BodySettings` architecture are:

        - **Approximate planet positions:** Highly simplified model of ephemerides of major Solar system bodies (model described here). Both a three-dimensional, and circular coplanar approximation may be used. Settings stored in ApproximatePlanetPositionSettings object. 

        .. code-block:: cpp

            bodySettings[ "Jupiter" ]->ephemerisSettings = boost::make_shared< ApproximatePlanetPositionSettings >( ephemerides::ApproximatePlanetPositionsBase::jupiter, false ); 

        where the first constructor argument is taken from the enum in approximatePlanetPositionsBase.h, and the second argument (false) denotes that the circular coplanar approximation is not made.

         - **Direct Spice ephemeris:** Ephemeris retrieved directly using the Spice toolbox (see this wiki page). Settings stored in DirectSpiceEphemerisSettings object. 

        .. code-block:: cpp

            std::string frameOrigin = "SSB";
            std::string frameOrientation = "J2000";
            bodySettings[ "Jupiter" ]->ephemerisSettings = boost::make_shared< DirectSpiceEphemerisSettings >( frameOrigin, frameOrientation ); 

        creating a barycentric (SSB) ephemeris with axes along J2000, with data directly from spice.

        - **Interpolated Spice ephemeris:** Using this option the state of the body is retrieved at regular intervals, and used to create an interpolator, before setting up environment. This has the advantage of only requiring calls to Spice outside of the propagation inner loop, reducing computation time. However, it has the downside of begin applicable only during a limited time interval. Settings stored in InterpolatedSpiceEphemerisSettings object. 

        .. code-block:: cpp

            double initialTime = 0.0;
            double finalTime = 1.0E8;
            double timeStep = 3600.0;
            std::string frameOrigin = "SSB";
            std::string frameOrientation = "J2000";
            bodySettings[ "Jupiter" ]->ephemerisSettings = boost::make_shared< InterpolatedSpiceEphemerisSettings >(
                initialTime, finalTime, timeStep, frameOrigin, frameOrientation ); 

        creating a barycentric (SSB) ephemeris with axes along J2000, with data retrieved from Spice at 3600 s intervals between t=0 and t=1.0E8, using a 6th order Lagrange interpolator. Settings for the interpolator (discussed here, can be added as a sixth argument if you wish to use a different interpolation method)

        - **Tabulated ephemeris:** Ephemeris created directly by interpolating user-specified states as a function of time. Settings stored in TabulatedEphemerisSettings object.

        .. code-block:: cpp

            std::map< double, Eigen::Vector6d > bodyStateHistory ...
            std::string frameOrigin = "SSB";
            std::string frameOrientation = "J2000";
            bodySettings[ "Jupiter" ]->ephemerisSettings = boost::make_shared< TabulatedEphemerisSettings >(
                bodyStateHistory, frameOrigin, frameOrientation ); 

        creating an ephemeris interpolated (with 6th order Lagrange interpolation) from the data in bodyStateHistory, which contains the Cartesian state (w.r.t. SSB; axes along J2000) for a given number of times (map keys).

        - **Kepler ephemeris:** Ephemeris modelled as being a perfect Kepler orbit. Settings stored in KeplerEphemerisSettings object.

        .. code-block:: cpp

            Eigen::Vector6d initialStateInKeplerianElements = ...
            double epochOfInitialState = ...
            double centralBodyGravitationalParameter = ...
            std::string frameOrigin = "SSB";
            std::string frameOrientation = "J2000";
            bodySettings[ "Jupiter" ]->ephemerisSettings = boost::make_shared< KeplerEphemerisSettings >(
                initialStateInKeplerianElements, epochOfInitialState, centralBodyGravitationalParameter, frameOrigin, frameOrientation ); 

        creating a Kepler orbit as ephemeris using the given kepler elements and associated initial time and gravitational parameter. See this page for more details on orbital elements in Tudat.

        - **Constant ephemeris:** Ephemeris modelled as being independent of time. Settings stored in ConstantEphemerisSettings object.

        .. code-block:: cpp

            Eigen::Vector6d constantCartesianState = ...
            std::string frameOrigin = "SSB";
            std::string frameOrientation = "J2000";
            bodySettings[ "Jupiter" ]->ephemerisSettings = boost::make_shared< ConstantEphemerisSettings >(
                constantCartesianState, frameOrigin, frameOrientation ); 

    **Gravity field model:** defined in the gravityFieldSettings (of type shared_ptr< GravityFieldSettings >). Models currently available through the BodySettings architecture are:

        - **Central (user-defined):** Point-mass gravity field model, with user-defined gravitational parameter. Settings stored in CentralGravityFieldSettings object.

        .. code-block:: cpp

            double gravitationalParameter = ...
            bodySettings[ "Earth" ]->gravityFieldSettings = boost::make_shared< CentralGravityFieldSettings >( gravitationalParameter );

        - **Central (from Spice):** Point-mass gravity field model, with gravitational parameter from Spice. No derived class of GravityFieldSettings base class required, created by passing central_spice as argument to base class constructor.

        .. code-block:: cpp

            bodySettings[ "Earth" ]->gravityFieldSettings = boost::make_shared< GravityFieldSettings >( central_spice ); 

        - **Spherical harmonic:** Gravity field model as a spherical harmonic expansion. Settings stored in SphericalHarmonicsGravityFieldSettings object. 

        .. code-block:: cpp

            double gravitationalParameter = ...
            double referenceRadius = ...
            Eigen::MatrixXd cosineCoefficients =  // NOTE: entry (i,j) denotes coefficient at degree i and order j
            Eigen::MatrixXd sineCoefficients =  // NOTE: entry (i,j) denotes coefficient at degree i and order j
            std::string associatedReferenceFrame = ...
            bodySettings[ "Earth" ]->gravityFieldSettings = boost::make_shared< SphericalHarmonicsGravityFieldSettings >( gravitationalParameter, referenceRadius, cosineCoefficients, sineCoefficients, associatedReferenceFrame ); 

        The associatedReferenceFrame reference frame must presently be the same frame as the target frame of the body's rotation model (see below). It represents the frame to which the spherical harmonic field is fixed.

    **Rotation model:** defined in the rotationModelSettings (of type shared_ptr< RotationModelSettings >). Models currently available through the BodySettings architecture are:

        - **Simple rotation model:** Rotation model with constant orientation of the rotation axis, and constant rotation rate about local z-axis. Settings stored in SimpleRotationModelSettings object.

        .. code-block:: cpp

            Eigen::Quaterniond initialOrientation = ...
            double initialTime = ...
            double rotationRate = ...
            std::string originalFrame = "J2000";
            std::string targetFrame = "IAU_Earth";
            bodySettings[ "Earth" ]->rotationModelSettings = boost::make_shared< SimpleRotationModelSettings >( 
                originalFrame, targetFrame , initialOrientation, initialTime, rotationRate ); 

        where the rotation from originalFrame to targetFrame at initialTime is given by the quaternion initialOrientation. This is mapped to other times using the rotation rate rotationRate.

        - **Spice rotation model:** Rotation model directly obtained from Spice. No derived class of RotationModelSettings base class required, created by passing spice_rotation_model as argument to base class constructor.

        .. code-block:: cpp

            std::string originalFrame = "J2000";
            std::string targetFrame = "IAU_Earth";
            bodySettings[ "Earth" ]->rotationModelSettings = boost::make_shared< RotationModelSettings >( spice_rotation_model, originalFrame, targetFrame ); 

    **Shape model:** defined in the shapeModelSettings (of type shared_ptr< BodyShapeSettings >). Models currently available through the BodySettings architecture are:

        - **Spherical shape model:** Model defining a body shape as a perfect sphere, with the sphere radius provided by the user. Settings stored in SphericalBodyShapeSettings object.

        .. code-block:: cpp

            double bodyRadius = 6378.0E3;
            bodySettings[ "Earth" ]->shapeModelSettings = boost::make_shared< SphericalBodyShapeSettings >( bodyRadius ); 

        - **Spherical shape model from spice:** Model defining a body shape as a perfect sphere, with the sphere radius retrieved from Spice. No derived class of BodyShapeSettings base class required, created by passing spherical_spice as argument to base class constructor.

        .. code-block:: cpp

            double bodyRadius = 6378.0E3;
            bodySettings[ "Earth" ]->shapeModelSettings = boost::make_shared< BodyShapeSettings >( spherical_spice ); 

        - **Oblate spheroid shape model:** Model defining a body shape as a flattened sphere, with the equatorial radius and flattening provided by the user. Settings stored in OblateSphericalBodyShapeSettings object.

        .. code-block:: cpp

            double bodyRadius = 6378.0E3;
            double bodyFlattening = 1.0 / 300.0;
            bodySettings[ "Earth" ]->shapeModelSettings = boost::make_shared< OblateSphericalBodyShapeSettings >( bodyRadius, bodyFlattening ); 

    **Radiation pressure model:** defined in the radiationPressureSettings (of type map< string, shared_ptr< RadiationPressureInterfaceSettings > >). A separate model can be used for different bodies emitting radiation (key valyes of radiationPressureSettings) Models currently available through the BodySettings architecture are:

        - **Cannon-ball model:** Properties for a cannonball radiation pressure model, i.e. effective force colinear with vector from source to target. Settings stored in CannonBallRadiationPressureInterfaceSettings object.

        .. code-block:: cpp

            std::string sourceBody = "Sun";
            double area = 20.0;
            const double radiationPressureCoefficient = 1.2;
            const std::vector< std::string > occultingBodies;
            occultingBodies.push_back( "Earth" );
            bodySettings[ "TestVehicle" ]->radiationPressureSettings[ sourceBody ] = boost::make_shared< CannonBallRadiationPressureInterfaceSettings >(
                sourceBody, area, radiationPressureCoefficient, occultingBodies ); 

        Creating cannonball radiation pressure settings for radiation due to the Sun, acting on the "TestVehicle" body, where the occultations due to the Earth are taken into account.

        .. note:: Occultations by multiple bodies are not yet supported. Please contact the Tudat suppport team if you wish to use multiple occultations.

    **Aerodynamic coefficients:** defined in the aerodynamicCoefficientSettings (of type shared_ptr< AerodynamicCoefficientSettings >). Models currently available through the BodySettings architecture are:

        - **Constant coefficients:** Settings for constant (not a function of any independent variables) aerodynamic coefficients. Settings stored in ConstantAerodynamicCoefficientSettings object.

        .. code-block:: cpp

            double referenceArea = 20.0;
            Eigen::Vector3d constantCoefficients;
            constantCoefficients( 0 ) = 1.5;
            constantCoefficients( 2 ) = 0.3;
            bodySettings[ "TestVehicle" ]->aerodynamicCoefficientSettings = boost::make_shared< ConstantAerodynamicCoefficientSettings >( 
                referenceArea, constantCoefficients, true, true ); 

        For constant drag coefficient of 1.5 and lift coefficient of 0.3.

        - **Tabulated:** Settings for tabulated aerodynamic coefficients as a function of given independent variables. Settings stored in TabulatedAerodynamicCoefficientSettings object. These tables can be defined either manually or loaded from a file, as discussed in more detail here Coefficients can be defined as a function of angle of sideslip, angle of attack, Mach number or altitude. If you simulation requires any other dependencies for the coefficients, please open an issue on Github requesting feature.

        - **Local inclination methods:** Settings for aerodynamic coefficients computed internally using a shape model of the vehicle, valid for hypersonic Mach numbers. Currently, this type of aerodynamic coefficients can only be set manually in the :class:`Body` object directly.

    **Gravity field variations:** defined in the gravityFieldVariationSettings (of type vector< shared_ptr< GravityFieldVariationSettings > >). Any number of gravity field variations may be used (hence the use of a vector). NOTE: You can only use gravity field variations for bodies where you have defined a spherical harmonic gravity field (through the use of SphericalHarmonicsGravityFieldSettings). Models currently available through the BodySettings architecture are:

        - **First-order solid body tide:** Tidal variation of the gravity field using first-order tidal theory. Settings stored in BasicSolidBodyGravityFieldVariationSettings object.

        - **Tabulated gravity field variations:** Variations in spherical harmonic coefficients tabulated as a function of time. Settings stored in TabulatedGravityFieldVariationSettings object.

The Environment During Propagation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Each :class:`Body` object and its constituent members is updated to the current state and time automatically during the numerical propagation. We stress that only those models that are relevant for a given propagation are updated every time step (this is handled automatically, without user intervention). Most time-dependent properties of the body are set in the environment models themselves. However, a number are updated and stored directly in the :class:`Body` object. These are:

    - The current translational state of the body
    - The current orientation of the body (and its time derivative)
    - The current mass of the body

.. note:: As a user, you will typically not access these variables directly.

The Environment Valid Time-Range
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Most of the environment models are valid for any time, but there is a key exception. In particular, the default settings do not directly use the Spice ephemerides, but retrieve the state for each body from Spice, and then create a :class:`TabulatedEphemeris` (which is only valid in the given time range), as opposed to a :class:`SpiceEphemeris`, which is valid for the entire time interval that the Spice kernels contain data. This approach is taken for computational reasons: retrieving a state from Spice is very time-consuming, much more so than retrieving it from a 6th- or 8th-order Lagrange interpolator that is used here for the tabulated ephemeris. An additional consequence of this is that the start and end time of the environment must be slightly (3 times the integration time step) larger than that which is used for the actual propagation, as a Lagrange interpolator can be unreliable at the edges of its domain. It is also possible to use the :class:`SpiceEphemeris` directly, at the expense of longer runtimes, by creating the ``bodySettings`` amd ``bodyMap`` as:

.. code-block:: cpp

    std::map< std::string, boost::shared_ptr< BodySettings > > bodySettings = getDefaultBodySettings( bodiesToCreate )
    NamedBodyMap bodyMap = createBodies( bodySettings );

