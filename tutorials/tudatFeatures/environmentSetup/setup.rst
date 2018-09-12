.. _tudatFeaturesEnvironmentSetUp:

Setting up the Environment
~~~~~~~~~~~~~~~~~~~~~~~~~~
In Tudat, the physical environment is defined by a list of :class:`Body` objects, each of which represents either a celestial body, or a manmade vehicle. Consequently, all properties that are required for computing e.g. accelerations are stored in :class:`Body` objects.

.. class:: Body
  
   Container for all properties of a body required for computing e.g. accelerations.

Typically the entire environment is stored in a named list of :class:`Body` object, the standard typedef for which is the :class:`NamedBodyMap`.

.. class:: NamedBodyMap

   An unordered map of shared pointers to :class:`Body` objects, see :ref:`this <externalBoost>` wiki page for a discussion of shared pointers; don't worry if you're not sure what a shared pointer or unordered map is at this point.

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

These properties can be set manually or default settings can be used. For instance, to manually create and set an ephemeris (from Spice w.r.t. the barycenter) and gravity field (point-mass only) object in the ``"Earth"`` entry of the body map, the following can be used:

.. code-block:: cpp

    bodyMap[ "Earth" ]->setEphemeris( boost::make_shared< SpiceEphemeris >( "Earth", "SSB", false, false, true, "J2000" ) ); 
    bodyMap[ "Earth" ]->setGravityFieldModel( boost::make_shared< GravityFieldModel >( 3.986004418E14 ) );  

This calls the constructors of the :class:`SpiceEphemeris` and :class:`GravityFieldModel` classes, and assigns the objects that are constructed to the "Earth" entry of the ``bodyMap``.

.. _tudatFeaturesCreatingTheEnvironment:

Creating the environment from :class:`BodySettings`
***************************************************
Manually creating all objects defining the full environment is possible, but not recommended. In particular, various environment models are interdependent and these dependencies must be fully and consistently defined for the code to function properly. To this end, we provide a :class:`BodySettings` object.

.. class:: BodySettings

   Class in which the general properties of each environment model can be set (see above for the list of the available types of environment models). We note that for :class:`Body` objects that represent vehicles, the manual creation is typically used, as the vehicle conditions may depend on the celestial bodies, but not vice versa.

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

The reasons for passing the initial/final time (as well as the 'buffer') are discussed in more detail at the end of this page. As can be seen from the above, the settings for the environment are stored in a map of pointers to :class:`BodySettings` objects (with the key the name of the associated bodies). If you have a look at the definition of the :class:`BodySettings` class (in ``SimulationSetup/createBodies.h``), you will see that this type is simply a container for a list of specific environment settings, which we discuss in more detail below. As a result, specifying settings for a given type of environment model requires the creation of an object of the correct type of class (derived class of :class:`EphemerisSettings` for defining an ephemeris; derived class of :class:`BodyShapeSettings` for defining a body shape etc.)

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