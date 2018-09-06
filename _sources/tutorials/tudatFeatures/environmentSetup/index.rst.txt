.. _tudatFeaturesEnvironmentIndex:

Environment Set-up
==================
Models for the physical environment are one of the cornerstones of a numerical astrodynamics toolbox. Here, we define the environment in the broadest sense, including all physical properties of the solar system, such as atmospheres and gravity fields, but also any models for the orbits and rotations of these bodies.

On this page, we will give an overview of how the environment is represented in Tudat, which models have been implemented, and how to create the environment that is tailored to your needs. A graphical representation of the basic structure for implementing the environment is shown below. For the options within these settings object, please click on the object in the figure.

.. graphviz::

   digraph
   {
      # General diagram settings
      rankdir = "LR";
      splines = ortho;    
      compound = true;  


      # general node settings 
      node [shape = box, style = filled, width = 1.25, fixedsize = true, color = lightgrey, fontname = FontAwesome, fontsize = 9];


      # specific node color settings
      NamedBodyMap [color = lightgreen];
      BodySettings [color = lightblue];


      # Hyperlinks (Sphinx auto referencing not working here, need to link to exact web adres)
      "GravityField\nVariationSettings" [href = "http://tudat.tudelft.nl/tutorials/tudatFeatures/environmentSetup/index.html#GravityFieldVariationSettings", target = "_top"];
      "RotationModelSettings" [href = "http://tudat.tudelft.nl/tutorials/tudatFeatures/environmentSetup/index.html#RotationModelSettings", target = "_top"];
      "Aerodynamic\nCoefficientSettings" [href = "http://tudat.tudelft.nl/tutorials/tudatFeatures/environmentSetup/index.html#AerodynamicCoefficientSettings", target = "_top"];
      "BodyShapeSettings" [href = "http://tudat.tudelft.nl/tutorials/tudatFeatures/environmentSetup/index.html#BodyShapeSettings", target = "_top"];
      "RadiationPressure\nInterfaceSettings" [href = "http://tudat.tudelft.nl/tutorials/tudatFeatures/environmentSetup/index.html#RadiationPressureInterfaceSettings", target = "_top"];
      "AtmosphereSettings" [href = "http://tudat.tudelft.nl/tutorials/tudatFeatures/environmentSetup/index.html#AtmosphereSettings", target = "_top"];
      "GravityFieldSettings" [href = "http://tudat.tudelft.nl/tutorials/tudatFeatures/environmentSetup/index.html#GravityFieldSettings", target = "_top"];
      "EphemerisSettings" [href = "http://tudat.tudelft.nl/tutorials/tudatFeatures/environmentSetup/index.html#EphemerisSettings", target = "_top"];
      BodySettings [href = "http://tudat.tudelft.nl/tutorials/tudatFeatures/environmentSetup/index.html#BodySettings", target = "_top"];
      NamedBodyMap [href = "http://tudat.tudelft.nl/tutorials/tudatFeatures/environmentSetup/index.html#NamedBodyMap", target = "_top"];


      # NamedBodyMap input
      BodySettings -> NamedBodyMap;
      "User-defined \nbodies" [style = dotted, fillcolor = lightgrey, color = black];
      "User-defined \nbodies" -> NamedBodyMap;
      {rank = same; BodySettings, NamedBodyMap};


      # BodySettings input
      EphemerisSettings -> BodySettings [ltail = clusterEnvironmentSettings];
      getDefaultBodySettings -> BodySettings;
      

      # getDefaultBodySettings input
      bodyNames -> getDefaultBodySettings;
      "(initial/final) time" -> getDefaultBodySettings;
      "(initial/final) time" [style = dotted, fillcolor = lightgrey, color = black];

		#point0 [shape = point, style = vis, width = 0.1];
      #point0 -> "setGlobalFrame\nBodyEphemerides";
      #point0 -> "NamedBodyMap"; 
      #BodySettings -> point0;

      # Cluster all environment settings derived classes
      subgraph clusterEnvironmentSettings
      {
         # cluster settings
         label = "User-defined environment settings";
         fontsize = 9;
         style = dashed;
         rank = min;

         # Make three collumns
         {rank = same; AtmosphereSettings, EphemerisSettings, GravityFieldSettings};
         {rank = same; "GravityField\nVariationSettings", "RadiationPressure\nInterfaceSettings"};
         {rank = same; "Aerodynamic\nCoefficientSettings", BodyShapeSettings, RotationModelSettings};

         BodyShapeSettings -> EphemerisSettings [style = invis];
         "RadiationPressure\nInterfaceSettings" -> BodyShapeSettings [style = invis];
         "GravityField\nVariationSettings" -> RotationModelSettings [style = invis];
      }
      "setGlobalFrame\nBodyEphemerides" -> NamedBodyMap [dir = both];
      "Frame origin\n& orientation" -> "setGlobalFrame\nBodyEphemerides";
		
		
      {"RadiationPressure\nInterfaceSettings", "GravityField\nVariationSettings" [fillcolor = lightcoral]};
   }

.. graphviz::

   digraph
   {
      # General diagram settings
      rankdir = "LR";
      splines = ortho;    
      compound = true;  

      subgraph clusterLegend
      {
      rank = min;
      style = dashed;


     	# general node settings 
     	node [shape = box, style = filled, width = 1.25, fixedsize = true, color = lightgrey, fontname = FontAwesome, fontsize = 9];


   	"List of settings" [ fillcolor = lightcoral];
     	"Main block" [fillcolor = lightgreen];
     	"Optional input" [style = dotted, fillcolor = lightgrey, color = black];
     	"Input for \nmain block" [fillcolor = lightblue];
     	"Optional input"-> "List of settings" -> "Input for \nmain block" -> "Main block" [style = invis];
      }
   }


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

Available Settings for the Environment Models
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Here, we will provide a full list of the available properties of the :class:`BodySettings` object. Each type of environment model has one base class to define settings for the creation of the model). Often, a specific derived class is implemented for a specific environment model of a given class, in which any additional information that may be needed can be provided. For instance, when defining a gravity field model, one can simply use:

.. code-block:: cpp

    bodySettings[ "Earth" ]->gravityFieldSettings = boost::make_shared< GravityFieldSettings >( central_spice ); 

if you want to use a central gravity field with the gravitational parameter taken from Spice: no information is needed except the type of gravity field model that is created. On the other hand, if you want to use a spherical harmonic gravity field, you need to specify additional parameters yourself, which is done by using the specific derived class:

.. code-block:: cpp

    bodySettings[ "Earth" ]->gravityFieldSettings = boost::make_shared< SphericalHarmonicsGravityFieldSettings >( gravitationalParameter, referenceRadius, cosineCoefficients, sineCoefficients, associatedReferenceFrame ); 

To find out which input arguments must be provided to create a specific settings class, have a look at the documentation in the code (written above the code for the constructor of the settings class you are interested in). Below, we give examples of each type of environment model setting.

The full list of available environment model settings is described below.

Atmosphere model
****************

.. class:: AtmosphereModel

   Base class for all atmosphere models. This model is constructed using the settings classes described below.

.. class:: AtmosphereSettings

   The base class for atmosphere settings. Models currently available through the :class:`BodySettings` architecture are (with examples when defining settings for Earth):
    
.. class:: ExponentialAtmosphereSettings

   Simple atmosphere model independent of time, latitude and longitude based on an exponentially decaying density profile with a constant temperature.

   .. code-block:: cpp

      bodySettings[ "Earth" ]->atmosphereSettings = boost::make_shared< ExponentialAtmosphereSettings >( 7.2E3, 290.0, 1.225, 287.06 ); 

   for an exponential atmosphere with a scale height of 7200 m, a constant temperature of 290 K, a density at 0 m altitude of 1.225 kg/m^3 and a specific gas constant of 287.06 J/(kg K).

.. class:: TabulatedAtmosphereSettings

   Atmosphere model with properties (pressure, density, temperature) read in from a file. Current implementation is independent of time, latitude and longitude. 

   .. code-block:: cpp

      std::string atmosphereFile = ...
      bodySettings[ "Earth" ]->atmosphereSettings = boost::make_shared< TabulatedAtmosphereSettings >( atmosphereFile ); 

   which will read the atmospheric properties from the file ``atmosphereFile`` (with four columns altitude and associated presure, density and temperature).

.. method:: NRLMSISE-00 
    
   This can be used to select the NRLMSISE-00 atmosphere model. To use this model, the :literal:`USE_NRLMSISE` flag in your top-level :literal:`CMakeLists` must be set to true. No derived class of :class:`AtmosphereSettings` base class required, the model can be created by passing :literal:`nrlmsise00` as argument to base class constructor. 

   .. code-block:: cpp

      bodySettings[ "Earth" ]->atmosphereSettings = boost::make_shared< AtmosphereSettings >( nrlmsise00 );  

.. class:: CustomWindModelSettings

   Custom wind model which can be used to retrieve a wind vector. This wind vector is in the body-fixed, body-centered reference frame. 

   .. code-block:: cpp
   
      bodySettings[ "Earth" ]->atmosphereSettings = boost::make_shared< CustomWindModelSettings >(  windFunction )
   
   where ``windFunction`` is a ``boost::function`` with inputs; altitude, longitude, latitude and time (for more details about boost: :ref:`externalBoost`).

.. _ephemerisModel:

Ephemeris model
****************  

.. class:: Ephemeris
  
   Base class for the ephemeris. It is constructed using one of the settings classes below.

.. class:: EphemerisSettings

   Base class for the ephemeris settings. Models currently available through the :class:`BodySettings` architecture and set by their respective derived classes are:

.. class:: ApproximatePlanetPositionSettings

   Highly simplified model of ephemerides of major Solar system bodies (model described here). Both a three-dimensional, and circular coplanar approximation may be used. 

   .. code-block:: cpp

       bodySettings[ "Jupiter" ]->ephemerisSettings = boost::make_shared< ApproximatePlanetPositionSettings >( ephemerides::ApproximatePlanetPositionsBase::jupiter, false ); 

   where the first constructor argument is taken from the enum in approximatePlanetPositionsBase.h, and the second argument (false) denotes that the circular coplanar approximation is not made.

.. class:: DirectSpiceEphermerisSettings

   Ephemeris retrieved directly using :ref:`tudatFeaturesSpice`.

   .. code-block:: cpp

       std::string frameOrigin = "SSB";
       std::string frameOrientation = "J2000";
       bodySettings[ "Jupiter" ]->ephemerisSettings = boost::make_shared< DirectSpiceEphemerisSettings >( frameOrigin, frameOrientation ); 

   creating a barycentric (SSB) ephemeris with axes along J2000, with data directly from spice.

.. class:: InterpolatedSpiceEphemerisSettings 
      
   Using this option the state of the body is retrieved at regular intervals, and used to create an interpolator, before setting up environment. This has the advantage of only requiring calls to Spice outside of the propagation inner loop, reducing computation time. However, it has the downside of begin applicable only during a limited time interval.

   .. code-block:: cpp

       double initialTime = 0.0;
       double finalTime = 1.0E8;
       double timeStep = 3600.0;
       std::string frameOrigin = "SSB";
       std::string frameOrientation = "J2000";
       bodySettings[ "Jupiter" ]->ephemerisSettings = boost::make_shared< InterpolatedSpiceEphemerisSettings >(
           initialTime, finalTime, timeStep, frameOrigin, frameOrientation ); 

   creating a barycentric (SSB) ephemeris with axes along J2000, with data retrieved from Spice at 3600 s intervals between t=0 and t=1.0E8, using a 6th order Lagrange interpolator. Settings for the interpolator (discussed here, can be added as a sixth argument if you wish to use a different interpolation method)

.. class:: TabulatedEphemerisSettings

   Ephemeris created directly by interpolating user-specified states as a function of time.

   .. code-block:: cpp

       std::map< double, Eigen::Vector6d > bodyStateHistory ...
       std::string frameOrigin = "SSB";
       std::string frameOrientation = "J2000";
       bodySettings[ "Jupiter" ]->ephemerisSettings = boost::make_shared< TabulatedEphemerisSettings >(
           bodyStateHistory, frameOrigin, frameOrientation ); 

   creating an ephemeris interpolated (with 6th order Lagrange interpolation) from the data in bodyStateHistory, which contains the Cartesian state (w.r.t. SSB; axes along J2000) for a given number of times (map keys, valid time range between first and last time in this map). 

.. class::  KeplerEphemerisSettings

   Ephemeris modelled as being a perfect Kepler orbit. 

   .. code-block:: cpp

       Eigen::Vector6d initialStateInKeplerianElements = ...
       double epochOfInitialState = ...
       double centralBodyGravitationalParameter = ...
       std::string frameOrigin = "SSB";
       std::string frameOrientation = "J2000";
       bodySettings[ "Jupiter" ]->ephemerisSettings = boost::make_shared< KeplerEphemerisSettings >(
           initialStateInKeplerianElements, epochOfInitialState, centralBodyGravitationalParameter, frameOrigin, frameOrientation ); 

   creating a Kepler orbit as ephemeris using the given kepler elements and associated initial time and gravitational parameter. See :ref:`tudatFeaturesFrameStateTransformations` for more details on orbital elements in Tudat.

.. class:: ConstantEphemerisSettings

   Ephemeris modelled as being independent of time.

   .. code-block:: cpp

       Eigen::Vector6d constantCartesianState = ...
       std::string frameOrigin = "SSB";
       std::string frameOrientation = "J2000";
       bodySettings[ "Jupiter" ]->ephemerisSettings = boost::make_shared< ConstantEphemerisSettings >(
           constantCartesianState, frameOrigin, frameOrientation ); 

.. method:: Multi-arc ephemeris

   An ephemeris model (for translational state) that allows the body’s state to be defined by distinct ephemeris models over different arcs. Class is implemented to support multi-arc propagation/estimation. No derived class of :class:`EphemerisSettings` base class required, the created ephemeris can be made multi-arc by using the ``resetMakeMultiArcEphemeris`` function of the :class:`EphemerisSettings` class. The resulting :class:`Ephemeris` object will then be :class:`MultiArcEphemeris` (with the same ephemeris model for each arc when created, according to the settings in the :class:`EphemerisSettings` object)

   .. code-block:: cpp

      bodySettings[ "Earth" ]->ephemerisSettings-> resetMakeMultiArcEphemeris( true );   

.. class:: CustomEphemerisSettings

   Allows user to provide arbitrary boost function as ephemeris model. 

   .. code-block:: cpp

      boost::shared_ptr< EphemerisSettings > customEphemerisSettings =
                   boost::make_shared< CustomEphemerisSettings >(
                      customBoostFunction, frameOrigin, frameOrientation );


Gravity field model
*******************

.. class:: GravityFieldModel

   Base class for the gravity field model, set using the settings classes described below.

.. class:: GravityFieldSettings

   Base class for the gravity field settings. Models currently available through the :class:`BodySettings` architecture can be called by the following:

.. class:: CentralGravityFieldSettings

   Point-mass gravity field model, with user-defined gravitational parameter. 

   .. code-block:: cpp

       double gravitationalParameter = ...
       bodySettings[ "Earth" ]->gravityFieldSettings = boost::make_shared< CentralGravityFieldSettings >( gravitationalParameter );

.. method:: Point-mass gravity field model from Spice

   Point-mass gravity field model, with gravitational parameter from Spice. No derived class of :class:`GravityFieldSettings` base class required, created by passing ``central_spice`` as argument to base class constructor.

   .. code-block:: cpp

       bodySettings[ "Earth" ]->gravityFieldSettings = boost::make_shared< GravityFieldSettings >( central_spice ); 

.. class:: SphericalHarmonicsGravityFieldSettings

   Gravity field model as a spherical harmonic expansion. 

   .. code-block:: cpp

       double gravitationalParameter = ...
       double referenceRadius = ...
       Eigen::MatrixXd cosineCoefficients =  // NOTE: entry (i,j) denotes coefficient at degree i and order j
       Eigen::MatrixXd sineCoefficients =  // NOTE: entry (i,j) denotes coefficient at degree i and order j
       std::string associatedReferenceFrame = ...
       bodySettings[ "Earth" ]->gravityFieldSettings = boost::make_shared< SphericalHarmonicsGravityFieldSettings >( gravitationalParameter, referenceRadius, cosineCoefficients, sineCoefficients, associatedReferenceFrame ); 

   The associatedReferenceFrame reference frame must presently be the same frame as the target frame of the body's rotation model (see below). It represents the frame to which the spherical harmonic field is fixed.

Rotational model
****************

.. class:: RotationalEphemeris

   Base class for the rotational ephemeris model, set using the settings classes described below.

.. class:: RotationModelSettings

   Base class for the rotational model settings. Models currently available through the :class:`BodySettings` architecture are:


.. class:: SimpleRotationModelSettings

   Rotation model with constant orientation of the rotation axis, and constant rotation rate about local z-axis. 

   .. code-block:: cpp

       Eigen::Quaterniond initialOrientation = ...
       double initialTime = ...
       double rotationRate = ...
       std::string originalFrame = "J2000";
       std::string targetFrame = "IAU_Earth";
       bodySettings[ "Earth" ]->rotationModelSettings = boost::make_shared< SimpleRotationModelSettings >( 
           originalFrame, targetFrame , initialOrientation, initialTime, rotationRate ); 

   where the rotation from originalFrame to targetFrame at initialTime is given by the quaternion initialOrientation. This is mapped to other times using the rotation rate rotationRate.

.. method:: Spice Rotation model

   Rotation model directly obtained from Spice. No derived class of :class:`RotationModelSettings` base class required, created by passing ``spice_rotation_model`` as argument to base class constructor.

   .. code-block:: cpp

       std::string originalFrame = "J2000";
       std::string targetFrame = "IAU_Earth";
       bodySettings[ "Earth" ]->rotationModelSettings = boost::make_shared< RotationModelSettings >( spice_rotation_model, originalFrame, targetFrame ); 

.. method:: Tabulated RotationalEphemeris model

   Rotation model obtained from an interpolator, with dependent variable a ``Eigen::VectorXd`` of size 7. Currently the settings interface is not yet implemented but the functionality is implemented in :class:`TabulatedRotationalEphemeris`. The tabulated rotational ephemeris can be implemented as follows:

   .. code-block:: cpp

      // Create tabulated rotational model
      boost::shared_ptr< TabulatedRotationalEphemeris< double, double > > tabulatedEphemeris =
              boost::make_shared< TabulatedRotationalEphemeris<  double, double > >( rotationInterpolator );

.. method:: Constant Rotation Model

   Rotation model with a constant value for the rotation. Currently the settings interface is not yet implemented. 

Body shape model
****************

.. class:: BodyShapeModel

   Base class for body shape models. It is constructed using the settings described below.

.. class:: BodyShapeSettings

   Base class for the body shape settings. Models currently available through the :class:`BodySettings` architecture are:

.. class:: SphericalBodyShapeSettings

   Model defining a body shape as a perfect sphere, with the sphere radius provided by the user. 

   .. code-block:: cpp

       double bodyRadius = 6378.0E3;
       bodySettings[ "Earth" ]->shapeModelSettings = boost::make_shared< SphericalBodyShapeSettings >( bodyRadius ); 

.. method:: Perfect sphere

   Model defining a body shape as a perfect sphere, with the sphere radius retrieved from Spice. No derived class of :class:`BodyShapeSettings` base class required, created by passing ``spherical_spice`` as argument to base class constructor.

   .. code-block:: cpp

       double bodyRadius = 6378.0E3;
       bodySettings[ "Earth" ]->shapeModelSettings = boost::make_shared< BodyShapeSettings >( spherical_spice ); 

.. class:: OblateSphericalBodyShapeSettings  

   Model defining a body shape as a flattened sphere, with the equatorial radius and flattening provided by the user. 

   .. code-block:: cpp

       double bodyRadius = 6378.0E3;
       double bodyFlattening = 1.0 / 300.0;
       bodySettings[ "Earth" ]->shapeModelSettings = boost::make_shared< OblateSphericalBodyShapeSettings >( bodyRadius, bodyFlattening ); 

Radiation pressure interface
****************************

.. class:: RadiationPressureInterface

   Class containing the properties of a solar radiation pressure acceleration model. It is constructed using the settings classes below. 

.. class:: RadiationPressureInterfaceSettings

   Base class for the radiation pressure interface settings. A separate model can be used for different bodies emitting radiation (key values of radiationPressureSettings) Models currently available through the :class:`BodySettings` architecture are:

.. class:: CannonBallRadiationPressureInterfaceSettings

   Properties for a cannonball radiation pressure model, i.e. effective force colinear with vector from source to target.

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

   .. _aerodynamicCoefficientOptions:


Aerodynamic coefficient interface
*********************************

.. class:: AerodynamicCoefficientInterface

   Base class containing the aerodynamic coefficient interface set by the settings classes below.

.. class:: AerodynamicCoefficientSettings

   Base class for the aerodynamic coefficient settings. Models currently available through the :class:`BodySettings` architecture are:
         
.. class:: ConstantAerodynamicCoefficientSettings

   Settings for constant (not a function of any independent variables) aerodynamic coefficients. 

   .. code-block:: cpp

       double referenceArea = 20.0;
       Eigen::Vector3d constantCoefficients;
       constantCoefficients( 0 ) = 1.5;
       constantCoefficients( 2 ) = 0.3;
       bodySettings[ "TestVehicle" ]->aerodynamicCoefficientSettings = boost::make_shared< ConstantAerodynamicCoefficientSettings >( 
           referenceArea, constantCoefficients, true, true ); 

   For constant drag coefficient of 1.5 and lift coefficient of 0.3.

.. class:: TabulatedAerodynamicCoefficientSettings

   Settings for tabulated aerodynamic coefficients as a function of given independent variables. These tables can be defined either manually or loaded from a file, as discussed in more detail :ref:`here <tudatFeaturesAerodynamicGuidanceReadingAerodynamicCoefficients>`. Coefficients can be defined as a function of angle of sideslip, angle of attack, Mach number or altitude. If you simulation requires any other dependencies for the coefficients, please open an issue on Github requesting feature.

.. method:: Local Inclination methods

   Settings for aerodynamic coefficients computed internally using a shape model of the vehicle, valid for hypersonic Mach numbers. Currently, this type of aerodynamic coefficients can only be set manually in the :class:`Body` object directly.

Time-variations of the gravity field
************************************

.. class:: GravityFieldVariations

   Virtual base class for spherical harmonic gravity field variations. Constructed using the settings classes below.

.. class:: GravityFieldVariationSettings

   Base class for the gravity field variation settings. Any number of gravity field variations may be used (hence the use of a vector). NOTE: You can only use gravity field variations for bodies where you have defined a spherical harmonic gravity field (through the use of :class:`SphericalHarmonicsGravityFieldSettings`). Models currently available through the :class:`BodySettings` architecture are:

.. class:: BasicSolidBodyGravityFieldVariationSettings

   Tidal variation of the gravity field using first-order tidal theory. 

.. class:: TabulatedGravityFieldVariationSettings

   Variations in spherical harmonic coefficients tabulated as a function of time. 

The Environment During Propagation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Each :class:`Body` object and its constituent members is updated to the current state and time automatically during the numerical propagation. We stress that only those models that are relevant for a given propagation are updated every time step (this is handled automatically, without user intervention). Most time-dependent properties of the body are set in the environment models themselves. However, a number are updated and stored directly in the :class:`Body` object. These are:

    - The current translational state of the body
    - The current orientation of the body (and its time derivative)
    - The current mass of the body

.. note:: As a user, you will typically not access these variables directly.

The Environment Valid Time-Range
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Most of the environment models are valid for any time, but there is a key exception. In particular, the default settings do not directly use the Spice ephemerides, but retrieve the state for each body from Spice, and then create a :class:`TabulatedEphemeris` (which is only valid in the given time range, of which settings are explained in :class:`TabulatedEphemerisSettings`), as opposed to a :class:`SpiceEphemeris` (as discussed in :class:`DirectSpiceEphermerisSettings`), which is valid for the entire time interval that the Spice kernels contain data. This approach is taken for computational reasons: retrieving a state from Spice is very time-consuming, much more so than retrieving it from a 6th- or 8th-order Lagrange interpolator that is used here for the tabulated ephemeris. An additional consequence of this is that the start and end time of the environment must be slightly (3 times the integration time step) larger than that which is used for the actual propagation, as a Lagrange interpolator can be unreliable at the edges of its domain. It is also possible to use the :class:`SpiceEphemeris` directly, at the expense of longer runtimes, by creating the ``bodySettings`` and ``bodyMap`` as:

.. code-block:: cpp

    std::map< std::string, boost::shared_ptr< BodySettings > > bodySettings = getDefaultBodySettings( bodiesToCreate )
    NamedBodyMap bodyMap = createBodies( bodySettings );

.. _tudatFeaturesAerodynamicCoefficients:

Aerodynamic Coefficients
~~~~~~~~~~~~~~~~~~~~~~~~
For a body that experiences aerodynamic forces, the aerodynamic coefficients of that body need to be defined. Aerodynamic coefficients are created in Tudat by creating an object of type :class:`AerodynamicCoefficientSettings`:

.. code-block:: cpp
    
    boost::shared_ptr< AerodynamicCoefficientSettings > aerodynamicCoefficientSettings = .....
    bodyMap[ "VehicleName" ]->setAerodynamicCoefficientInterface(
                        createAerodynamicCoefficientInterface( aerodynamicCoefficientSettings, "VehicleName" ) );

The aerodynamic coefficients can be implemented in several different ways that are discussed in this section.

Local Incliniation Methods
**************************
One of the options to determine the aerodynamic coefficients during hypersonic flight is to use a hypersonic local inclination method. These methods calculate the aerodynamic coefficients using the orientation of the vehicle with respect to the incoming flow. A selection of these methods are implemented into Tudat and are able to calculate the aerodynamic coefficients for several pre-defined geometries and user defined grids.

The class which can be used is defined below.

.. class:: HypersonicLocalInclinationAnalysis

The settings for this class are constructed as follows:

.. code-block:: cpp
    
    HypersonicLocalInclinationAnalysis(
            const std::vector< std::vector< double > >& dataPointsOfIndependentVariables,
            const boost::shared_ptr< SurfaceGeometry > inputVehicleSurface,
            const std::vector< int >& numberOfLines,
            const std::vector< int >& numberOfPoints,
            const std::vector< bool >& invertOrders,
            const std::vector< std::vector< int > >& selectedMethods,
            const double referenceArea,
            const double referenceLength,
            const Eigen::Vector3d& momentReferencePoint );

where:

- :literal:`dataPointsOfIndependentVariables`:
      A 2-dimensional vector of doubles that contains the data points of each independent variable. Index 0 corresponds to the Mach number, index 1 to the angle of attack, and index 2 to the
      sideslip angle.

- :literal:`inputVehicleSurface`:
      A :literal:`SurfaceGeometry` type variable that defines the geometry of the vehicle used for the determination of the aerodynamic coefficients. 

- :literal:`numberOfLines`:
      The number of discretization points in the first independent surface for each subpart of the vehicle defined in :literal:`inputVehicleSurface`. 

- :literal:`numberOfPoints`:
      The number of discretization points in the second independent surface for each subpart of the vehicle defined in :literal:`inputVehicleSurface`. 

- :literal:`invertOrders`:
      A boolean that determines the orientation of the normal vectors of the panels (the discretized units) of each subpart of the vehicle defined in :literal:`inputVehicleSurface`. This can either be
      inward- or outward-facing.

- :literal:`selectedMethods`:
      A two dimensional vector that selects which specific method is used for the local inclination analysis. The first index determines if a compression (0) or a expansion (1) method is used, and
      the second index determines the specific method (a list of available methods is given below).

- :literal:`referenceArea`:
      A :literal:`double` that gives the reference area of the vehicle that is used for the analysis.

- :literal:`referencelength`:
      A :literal:`double` that gives the reference length of the vehicle that is used for the analysis.

- :literal:`momentReferencePoint`:
      A :literal:`Eigen::Vector3d` that gives the location of the moment reference point, which is used as a reference point to calculate the moments acting on the vehicle.


The :literal:`HypersonicLocalInclinationAnalysis` has several methods incorporated in it, that can be selected by the user through the :literal:`selectedMethods` setting. The first index of this vector determines which kind of method is used: a compression or expansion method. The second index then determines the specific method.

For the compression methods, the following are available:

- 0: Newtonian Method.
- 1: Modified Newtonian.
- 2 and 3: not available at this moment.
- 4: tangent-wedge method.
- 5: tangent-cone method.
- 6: modified Dahlem-Buck method.
- 7: VanDyke unified pressure method.
- 8: Smyth Delta Wing method.
- 9: Hankey flat surface method.

The expansion method has the following options:

- 0: Vacuum Pressure coefficient method.
- 1: Zero Pressure function.
- 4: High Mach base pressure method.
- 3 or 5: Prandtl-Meyer method.
- 6: ACM empirical pressure coefficient. 

To get a better idea of how this can be implemented, an example contained in the aerodynamics unit test of Tudat is used. This example describes how the aerodynamic coefficients are generated for the Apollo capsule. A function called :literal:`getApolloCoefficientInterface()` is made which returns a variable of the type: :literal:`HypersonicLocalInclinationAnalysis`.

First the capsule geometry is made using a pre-defined capsule shape:

.. code-block:: cpp
    
    boost::shared_ptr< geometric_shapes::Capsule > capsule
            = boost::make_shared< geometric_shapes::Capsule >(
               4.694, 1.956, 2.662, -1.0 * 33.0 * PI / 180.0, 0.196 );

There is also an option to define a grid by the user, but this is not explained here. To define the panels on the vehicle, which are used for the calculation of the local inclination, the number of lines and points are easily definable:

.. code-block:: cpp
    
    std::vector< int > numberOfLines;
    std::vector< int > numberOfPoints;
    std::vector< bool > invertOrders;
    numberOfLines.resize( 4 );
    numberOfPoints.resize( 4 );
    invertOrders.resize( 4 );

    // Set number of analysis points.
    numberOfLines[ 0 ] = 31;
    numberOfPoints[ 0 ] = 31;
    numberOfLines[ 1 ] = 31;
    numberOfPoints[ 1 ] = 31;
    numberOfLines[ 2 ] = 31;
    numberOfPoints[ 2 ] = 10;
    numberOfLines[ 3 ] = 11;
    numberOfPoints[ 3 ] = 11;
    invertOrders[ 0 ] = 0;
    invertOrders[ 1 ] = 0;
    invertOrders[ 2 ] = 0;
    invertOrders[ 3 ] = 0;

A moment reference point is then defined. Afterwards, the independent variable points are defined, which are used to calculate the coefficients in the simulation during specific flight conditions. These points can be defined by the user, but there is also an option to get default values:

.. code-block:: cpp
    
    std::vector< std::vector< double > > independentVariableDataPoints;
    independentVariableDataPoints.resize( 3 );
    independentVariableDataPoints[ 0 ] = getDefaultHypersonicLocalInclinationMachPoints( "Full" );
    std::vector< double > angleOfAttackPoints;
    angleOfAttackPoints.resize( 15 );

    for ( int i = 0; i < 15; i++ )
    {
        angleOfAttackPoints[ i ] = static_cast< double >( i - 6 ) * 5.0 * PI / 180.0;
    }

    independentVariableDataPoints[ 1 ] = angleOfAttackPoints;
    independentVariableDataPoints[ 2 ] =
    getDefaultHypersonicLocalInclinationAngleOfSideslipPoints( );

For the Mach number and angle of sideslip, default values are used, whereas for the angle of attack, the values are defined by the user. The final part of the code is the selection of methods. In this case there are several methods used to determine the coefficients. 

.. code-block:: cpp
    
    std::vector< std::vector< int > > selectedMethods;
    selectedMethods.resize( 2 );
    selectedMethods[ 0 ].resize( 4 );
    selectedMethods[ 1 ].resize( 4 );

    selectedMethods[ 0 ][ 0 ] = 1;
    selectedMethods[ 0 ][ 1 ] = 5;
    selectedMethods[ 0 ][ 2 ] = 5;
    selectedMethods[ 0 ][ 3 ] = 1;
    selectedMethods[ 1 ][ 0 ] = 6;
    selectedMethods[ 1 ][ 1 ] = 3;
    selectedMethods[ 1 ][ 2 ] = 3;
    selectedMethods[ 1 ][ 3 ] = 3;


Finally, in the return statement, the local inclination analysis is made, which can be used to define an :literal:`AerodynamicCoefficientInterface`:

.. code-block:: cpp
    
    // Create analysis object and capsule database.
    return boost::make_shared< HypersonicLocalInclinationAnalysis >(
                independentVariableDataPoints, capsule, numberOfLines, numberOfPoints,
                invertOrders, selectedMethods, PI * pow( capsule->getMiddleRadius( ), 2.0 ),
                3.9116, momentReference ); 

.. _tudatFeaturesAerodynamicGuidanceReadingAerodynamicCoefficients:
 
Reading aerodynamic coefficients from Files
*******************************************
For many simulations/analyses involving atmospheric flight, the aerodynamic coefficients will be provided in tabulated form. If put into the correct file format, these files can be read into Tudat used during the orbit propagation. By specifying the physical meaning of the independent variables of the aerodynamic coefficients, no action on the side of the user is required to update the aerodynamic coefficients to their correct values during propagation. Here, we give an overview and some examples on how to load aerodynamic coefficients from a file.

As a reminder, aerodynamic coefficients are created in Tudat by creating an object of type :class:`AerodynamicCoefficientSettings`:

.. code-block:: cpp
    
    boost::shared_ptr< AerodynamicCoefficientSettings > aerodynamicCoefficientSettings = .....
    bodyMap[ "VehicleName" ]->setAerodynamicCoefficientInterface(
                        createAerodynamicCoefficientInterface( aerodynamicCoefficientSettings, "VehicleName" ) );

To create an :class:`AerodynamicCoefficientSettings` object from data in files, we provide two functions named :literal:`readTabulatedAerodynamicCoefficientsFromFiles` (in :literal:`createFlightConditions.h`). For one of the functions, only force coefficients are loaded (with moment coefficients set to zero at all times). The other function allows both force and moment coefficients to be loaded.

For the situation where only force coefficients are considered, several pieces information are needed:

   - A list of files for any of the three aerodynamic coefficients (e.g. C\ :sub:`D`, C\ :sub:`S`, C\ :sub:`L` or C\ :sub:`X`, C\ :sub:`Y`, C\ :sub:`Z`). Note that the behaviour of each coefficient must be provided in a separate file. Note that not every coefficient needs to be defined. If a file is not provided for one of the coefficients, as will often be the case for C\ :sub:`S`, zeros are assumed at all points in the propagation.
   - The physical meaning of each of the independent variables of the coefficients.
   - The reference area for the aerodynamics. This is not read from the file and must be provided as an input to the function.
   - Two booleans denoting the orientation and direction of the aerodynamic coefficients. For instance C\ :sub:`D`, C\ :sub:`S`, C\ :sub:`L` denote the strength of the aerodynamic force in the aerodynamic reference frame, in a direction opposite to the axes of that frame. The C\ :sub:`X`, C\ :sub:`Y`, C\ :sub:`Z` coefficients, on the other hand, are defined in the body-fixed frame.

As an example, the following can be used to create :class:`AerodynamicCoefficientSettings` for force coefficients only from a file:

    .. code-block:: cpp
    
        double referenceArea = 50.0; // Define reference area

        // Define physical meaning of independent variables, in this case Mach number and angle of attack
        std::vector< aerodynamics::AerodynamicCoefficientsIndependentVariables > independentVariableNames; 
        independentVariableNames.push_back( aerodynamics::mach_number_dependent );
        independentVariableNames.push_back( aerodynamics::angle_of_attack_dependent ); 

        // Define list of files for force coefficients. Entry 0 denotes the x-direction (C ~D~/C ~X~), 1 the y-direction (C ~S~/C ~Y~) and 2 the z-direction (C ~L~/C ~Z~)
        std::map< int, std::string > forceCoefficientFiles; 
        forceCoefficientFiles[ 0 ] = tudat::input_output::getTudatRootPath( ) + "Astrodynamics/Aerodynamics/UnitTests/aurora_CD.txt"; // Set drag coefficient file
        forceCoefficientFiles[ 2 ] = tudat::input_output::getTudatRootPath( ) + "Astrodynamics/Aerodynamics/UnitTests/aurora_CL.txt"; // Set lift coefficient file

        // Define reference frame in which the loaded coefficients are defined.
        bool areCoefficientsInAerodynamicFrame = true;
        bool areCoefficientsInNegativeAxisDirection = true;

        // Load and parse files; create coefficient settings.
        boost::shared_ptr< AerodynamicCoefficientSettings > aerodynamicCoefficientSettings = 
            readTabulatedAerodynamicCoefficientsFromFiles( forceCoefficientFiles, referenceArea, independentVariableNames, areCoefficientsInAerodynamicFrame,         areCoefficientsInNegativeAxisDirection );

        // Create and set aerodynamic coefficients
        bodyMap[ "VehicleName" ]->setAerodynamicCoefficientInterface(
                            createAerodynamicCoefficientInterface( aerodynamicCoefficientSettings, "VehicleName" ) );

Note that in the above, no side force coefficient file (entry 1 for forceCoefficientFiles) is given, so that C\ :sub:`S`=0 always. Two independent variables have been defines (Mach number and angle of attack). If either the lift or drag coefficient files encounter a different number of independent variables, the program will terminate with an appropriate error message. Also, if the independent variables used for the lift and drag coefficients are not identical, the program is terminated.

.. tip:: Moment coefficients are added in a completely analogous manner (with separate files for the x-, y- and z-components).

In addition to defining aerodynamic coefficients for the vehicle itself, the influence of control surface deflections on the values of the coefficients will be needed for certain applications. In Tudat, any number of control surfaces may be defined for a vehicle, the deflection of which may be set by your particular :class:`AerodynamicGuidance` derived class. Loading the aerodynamic coefficient increments of the control surface is done in a manner similar to those of the total vehicle, but:

    - Reference area and reference frame in which the coefficients are defined are not provided. These are implcitily assumed to be equal to those of the aerodynamic coefficients of the 'main body'. If your application requires these quantities to be different for the body and control surface deflections, please open an issue on Github requesting the functionality.
    - Exactly one of the independent variables of the coefficient increments must be a control surface deflection.

Below, an example is given on how to load the aerodynamic coefficient increments:

.. code-block:: cpp
    
    // Create coefficient settings for body.
    boost::shared_ptr< AerodynamicCoefficientSettings > aerodynamicCoefficientSettings = ...

    // Define physical meaning of independent variables for control surface increments, in this case Mach number, angle of attack and control surface deflection
    std::vector< aerodynamics::AerodynamicCoefficientsIndependentVariables > controlSurfaceIndependentVariableNames; 
    controlSurfaceIndependentVariableNames.push_back( aerodynamics::mach_number_dependent );
    controlSurfaceIndependentVariableNames.push_back( aerodynamics::angle_of_attack_dependent ); 
    controlSurfaceIndependentVariableNames.push_back( aerodynamics::control_surface_deflection_dependent ); 

    // Define name of control surface
    std::string controlSurfaceName = "Elevon";

    // Define list of files for force coefficients. 
    std::map< int, std::string > controlSurfaceForceCoefficientFiles; 
    controlSurfaceForceCoefficientFiles[ 0 ] = tudat::input_output::getTudatRootPath( ) + "Astrodynamics/Aerodynamics/UnitTests/dCDwTest.txt"; // Set drag coefficient file

    // Add settings for control surface increments to main aerodynamic coefficients
    aerodynamicCoefficientSettings->setControlSurfaceSettings( 
        readTabulatedControlIncrementAerodynamicCoefficientsFromFiles( controlSurfaceForceCoefficientFiles, controlSurfaceIndependentVariableNames, controlSurfaceName ) );

    // Create and set aerodynamic coefficients
    bodyMap[ "VehicleName" ]->setAerodynamicCoefficientInterface(
        createAerodynamicCoefficientInterface( aerodynamicCoefficientSettings, "VehicleName" ) );

For this example, only the drag coefficient is affected by the control surface deflections.

.. _tudatFeaturesCustomAerodynamicCoefficients:
 
Custom Aerodynamic Coefficient Settings
****************************************
If a user specific aerodynamic coefficient interface is needed, the :literal:`CustomAerodynamicCoefficientInterface` can be used. This class allows a generic aerodynamic and moment coefficient function input, and allows the user to use all the aerodynamic coefficient interface methods. 

.. class:: CustomAerodynamicCoefficientInterface

The constructor for this class looks as follows:

.. code-block:: cpp
    
    CustomAerodynamicCoefficientInterface(
            const boost::function< Eigen::Vector3d( const std::vector< double >& ) >
            forceCoefficientFunction,
            const boost::function< Eigen::Vector3d( const std::vector< double >& ) >
            momentCoefficientFunction,
            const double referenceLength,
            const double referenceArea,
            const double lateralReferenceLength,
            const Eigen::Vector3d& momentReferencePoint,
            const std::vector< AerodynamicCoefficientsIndependentVariables >
            independentVariableNames,
            const bool areCoefficientsInAerodynamicFrame = true,
            const bool areCoefficientsInNegativeAxisDirection = true )

where:

- :literal:`forceCoefficientFunction`:
      A function outputting the aerodynamic force coefficients as a function of the independent variables, listed in :literal:`independentVariableNames`. This requires a :literal:`boost::function
      and :literal:`boost::bind`, which is explained in: :ref:`externalBoost`.

- :literal:`momentCoefficientFunction`:
      A function outputting the aerodynamic moment coefficients as a function of the independent variables, listed in :literal:`independentVariableNames`. This requires a :literal:`boost::function
      and :literal:`boost::bind`, which is explained in: :ref:`externalBoost`. 

- :literal:`referenceLength`:
      The length with which the coefficients are non-dimensionalized (about the z- and x-axis). 

- :literal:`referenceArea`:
      The area with which the coefficients are non-dimensionalized. 

- :literal:`lateralReferenceLength`:
      The length with which the coefficients are non-dimensionalized (about the y-axis). 

- :literal:`momentReferencePoint`:
      The point with respect to the aerodynamic moments are calculated.

- :literal:`independentVariableNames`:
      A vector containing identifiers for the independent variables that are used in the :literal:`forceCoefficientFunction`.

- :literal:`areCoefficientsInAerodynamicFrame`:
      A :literal:`bool` that determines if the aerodynamic coefficients are in the aerodynamic reference frame or not.

- :literal:`areCoefficientsInNegativeAxisDirection`:
      A :literal:`bool` that determines if the aerodynamic coefficients are in the negative axis direction.


