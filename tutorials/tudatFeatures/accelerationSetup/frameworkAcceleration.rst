.. _tudatFeaturesFrameworkAccelerations:

Acceleration Model Set-up
=========================

Below is a schemaric overview of the various acceleration options in Tudat, and the top-level architecture for how to set up acceleration models.

Acceleration settings
~~~~~~~~~~~~~~~~~~~~~
The settings for accelerations are defined and stored by the user in a :class:`SelectedAccelerationMap`. 

.. class:: SelectedAccelerationMap

   This is a  :literal:`typedef` for:

   .. code-block:: cpp

         std::map< std::string, std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > >

   This is a double map (with twice a string as a key). The two levels correspond to the names of bodies undergoing an acceleration (first key) , and those for bodies exerting an acceleration (second key). This allows any number of bodies to be propagated, undergoing any number (and type) of accelerations from any number of bodies. In this manner, settings for each required acceleration model are stored in an object of type :class:`AccelerationSettings`.
 
For a given environment, most acceleration models are completely defined by:

    - Type of acceleration model (a list is provided in the :class:`AvailableAcceleration` enum in ``Tudat/Astrodynamics/BasicAstrodynamics/accelerationModelTypes.h``).
    - Name of body undergoing acceleration
    - Name of body exerting acceleration

For instance, when using the following from the :ref:`walkthroughsUnguidedCapsuleEntry`:

.. code-block:: cpp

    SelectedAccelerationMap accelerationSettings;
    accelerationSettings[ "Apollo" ][ "Earth" ].push_back( std::make_shared< AccelerationSettings >( central_gravity ) );
    accelerationSettings[ "Apollo" ][ "Earth" ].push_back( std::make_shared< AccelerationSettings >( aerodynamic ) );

We have defined a point-mass Earth gravity model and an aerodynamic acceleration due to Earth's atmosphere to be used. In this example, we have only defined the type of the acceleration, without the need for any additional information. All required variables used on the computations of the accelerations are uniquely defined in the Apollo and Earth entries of the body map (provided that the required environment models have been set as discussed in :ref:`tudatFeaturesEnvironmentIndex`).

Below a graphical representation of the acceleration setup, with the different model types and top-level architecture.

.. graphviz::

   digraph
   {
      # General diagram settings
      rankdir = "TB";
      splines = ortho;    
      compound = true;   


      # general node settings 
      node [shape = box, style = filled, width = 1.5, fixedsize = true, color = lightgrey, fontname = FontAwesome, fontsize = 9];


      # specific node color settings
      NamedBodyMap, bodiesToPropagate, centralBodies, AccelerationSettings [color = lightblue];
      AccelerationModel [color= lightgreen];
      Empirical, "Relativistic correction", "Mutual spherical \nharmonic gravity", "Tidal Dissipation", "Thrust" [color =  darkturquoise]; 
      "Additional \ninformation" [style = dotted, fillcolor = lightgrey, color = black];

      subgraph clusterAccelerationType
      {
         label = AccelerationType;
         fontsize = 9;
         style = dashed;
      
         {rank = same; "Cannon ball \nradiation pressure", "Central gravity", Aerodynamic};
         {rank = same; "Thrust", "Tidal Dissipation", "Spherical harmonic \ngravity" };
         {rank = same; "Mutual spherical \nharmonic gravity", "Relativistic correction", Empirical};
         
         "Mutual spherical \nharmonic gravity" -> "Spherical harmonic \ngravity" -> "Central gravity" [style = invis];
         Empirical -> Thrust -> Aerodynamic [style = invis];
         "Cannon ball \nradiation pressure" -> "Central gravity" [style = invis];

      }
      
      Aerodynamic -> "Additional \ninformation"-> NamedBodyMap [style = invis];
      "Mutual spherical \nharmonic gravity" -> bodiesToPropagate [style = invis];
      AccelerationSettings -> bodiesToPropagate [style = invis];

      # AccelerationSettings input
      "Additional \ninformation" -> AccelerationSettings;
      "Central gravity" -> AccelerationSettings [ltail = clusterAccelerationType];    

      
      # AccelerationModel input
      AccelerationSettings -> AccelerationModel;
      bodiesToPropagate -> AccelerationModel;
      centralBodies -> AccelerationModel;
      NamedBodyMap -> AccelerationModel;
     

      # Structure the layout
      {rank = same; NamedBodyMap, AccelerationModel, centralBodies, bodiesToPropagate};
      {rank = same; AccelerationSettings, "Additional \ninformation"};

      # Hyperlinks (Sphinx auto referencing not working here, need to link to exact web adres)
      "AccelerationSettings" [href = "http://tudat.tudelft.nl/tutorials/tudatFeatures/accelerationSetup/frameworkAcceleration.html#AccelerationSettings", target = "_top"];
      "Central gravity" [href = "http://tudat.tudelft.nl/tutorials/tudatFeatures/accelerationSetup/frameworkAcceleration.html#AccelerationSettings", target = "_top"];
      "AccelerationModel" [href = "http://tudat.tudelft.nl/tutorials/tudatFeatures/accelerationSetup/frameworkAcceleration.html##SelectedAccelerationMap", target = "_top"];
      NamedBodyMap [href = "http://tudat.tudelft.nl/tutorials/tudatFeatures/environmentSetup/index.html#NamedBodyMap", target = "_top"];
      "Spherical harmonic \ngravity" [href = "http://tudat.tudelft.nl/tutorials/tudatFeatures/accelerationSetup/frameworkAcceleration.html#SphericalHarmonicAccelerationSettings", target = "_top"];
      "Mutual spherical \nharmonic gravity" [href = "http://tudat.tudelft.nl/tutorials/tudatFeatures/accelerationSetup/frameworkAcceleration.html#MutualSphericalHarmonicAccelerationSettings", target = "_top"];
      "Aerodynamic" [href = "http://tudat.tudelft.nl/tutorials/tudatFeatures/accelerationSetup/frameworkAcceleration.html#MutualSphericalHarmonicAccelerationSettings", target = "_top"];
      "Cannon ball \nradiation pressure" [href = "http://tudat.tudelft.nl/tutorials/tudatFeatures/accelerationSetup/frameworkAcceleration.html#MutualSphericalHarmonicAccelerationSettings", target = "_top"];
      "Thrust" [href = "http://tudat.tudelft.nl/tutorials/tudatFeatures/accelerationSetup/frameworkAcceleration.html#ThrustAccelerationSettings", target = "_top"];
      "Tidal Dissipation" [href = "http://tudat.tudelft.nl/tutorials/tudatFeatures/accelerationSetup/frameworkAcceleration.html#DissipationAccelerationSettings", target = "_top"];
      "Relativistic correction" [href = "http://tudat.tudelft.nl/tutorials/tudatFeatures/accelerationSetup/frameworkAcceleration.html#RelativisticAccelerationCorrectionSettings", target = "_top"];
      "Empirical" [href = "http://tudat.tudelft.nl/tutorials/tudatFeatures/accelerationSetup/frameworkAcceleration.html#EmpiricalAccelerationSettings", target = "_top"];

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


   	"Object requiring \nadditional information" [ fillcolor = darkturquoise];
     	"Main block" [fillcolor = lightgreen];
     	"Optional input" [style = dotted, fillcolor = lightgrey, color = black];
     	"Input for \nmain block" [fillcolor = lightblue];
     	"Optional input"-> "Object requiring \nadditional information" -> "Input for \nmain block" -> "Main block" [style = invis];
      }
   }

As was the case for the settings of the environment models, certain types of accelerations require additional information. An important example is the spherical harmonic acceleration. We cannot replace :literal:`central_gravity` with :literal:`spherical_harmonic_gravity` in the above, as there is still an ambiguity in how the acceleration model is defined. In particular, we now also need the maximum degree and order of the gravity field that is to be used in addition to the three properties listed above. Consequently, we have created a dedicated :class:`AccelerationSettings` derived class for defining spherical harmonic acceleration settings: :class:`SphericalHarmonicAccelerationSettings`. Updating the above example to use J\ :sub:`2`, J\ :sub:`3` and J\ :sub:`4` (maximum degree = 4; maximum order = 0), we now have:

.. code-block:: cpp

    SelectedAccelerationMap accelerationSettings;
    accelerationSettings[ "Apollo" ][ "Earth" ].push_back( std::make_shared< SphericalHarmonicAccelerationSettings >( 4, 0 ) );
    accelerationSettings[ "Apollo" ][ "Earth" ].push_back( std::make_shared< AccelerationSettings >( aerodynamic ) );

A full list of the available acceleration models, as well as their required input and environment models, is given at the end of this page. 

.. note:: When creating an object of the :class:`AccelerationSettings` type (or its derived class), you must not provide any of the third body acceleration types (:literal:`third_body_central_gravity`, :literal:`third_body_spherical_harmonic_gravity`, :literal:`third_body_mutual_spherical_harmonic_gravity`) as input. If you wish to use a third-body gravity acceleration (typically from a point mass), simply provide :literal:`central_gravity` as input. Depending on the settings for your central bodies, the code will automatically create the corresponding acceleration object (central or third-body).


Having defined all the required settings for the accelerations in your :class:`SelectedAccelerationMap`, you can create the actual acceleration models by using the :literal:`createAccelerationModelsMap` function. This function requires four input parameters:

    - Full environment, as defined by a :class:`NamedBodyMap`.
    - Settings for the acceleration models, given by :class:`SelectedAccelerationMap`.
    - A list of bodies to numerically propagate.
    - A list of central bodies (one for each numerically propagated body).

The list of central bodies defines the reference frame origins in which the bodies are propagated. The use of a hierarchical system is perfectly acceptable. For instance, one can propagate the Earth and Mars w.r.t. the Sun, the Sun w.r.t. the barycenter, the Moon w.r.t the Earth, etc. For this case, the central bodies and propagated bodies are defined as:

.. code-block:: cpp

    std::vector< std::string > propagatedBodies;
    propagatedBodies.push_back( "Moon" );
    propagatedBodies.push_back( "Earth" );
    propagatedBodies.push_back( "Mars" );
    propagatedBodies.push_back( "Sun" );

    std::vector< std::string > centralBodies;
    centralBodies.push_back( "Earth" );
    centralBodies.push_back( "Sun" );
    centralBodies.push_back( "Sun" );
    centralBodies.push_back( "SSB" );

There is no hardcoded limit to the number of permitted levels in the frame hierarchy, but it is not allowed to include circular dependencies, i.e. body A w.r.t. body B, body B w.r.t. body C and body C w.r.t. body A. More information of the acceleration models is discussed in :ref:`tudatFeaturesPropagatorSettings`. The following gives an example on how to create the acceleration model objects:

.. code-block:: cpp

    NamedBodyMap bodyMap;
    ....
    // Create environment here
    ....
    std::vector< std::string > propagatedBodies;
    std::vector< std::string > centralBodies;
    ....
    // Set central and propagated bodies here
    ....
    AccelerationMap accelerationModelMap = createAccelerationModelsMap( bodyMap, accelerationMap, propagatedBodies, centralBodies );

Mutual acceleration between bodies being propagated (i.e body A exerting acceleation on body B and vice versa), as is the case for solar system dynamics, is automatically handled by the :literal:`createAccelerationModelsMap` code and requires no specific consideration. Moreover, when creating a gravitational acceleration, the code checks whether it is a direct or a third-body gravitational acceleration and creates the acceleration models accordingly. Similarly, the code automatically checks which value of the gravitational parameter "mu" to use in such computations. For instance, when computing the gravitational acceleration due to the Sun acting on the Earth, :literal:`mu_Sun` is used when propagating w.r.t. the barycenter, whereas :literal:`mu_Sun + mu_Earth` is used when propagating w.r.t. the Sun.

For every acceleration, a model for the current state of the body exerting the acceleration must be available (the state of the body undergoing the acceleration is taken from the numerically propagated state). This means that, in the above example of the Apollo capsule entering Earth's atmosphere (:ref:`walkthroughsUnguidedCapsuleEntry`), we must include one of the following:

    - An ephemeris member for Earth.
    - Numerically integrate the Earth concurrently with our Apollo vehicle.

For this example, the second option is of course a bit 'non-standard'. However, for cases where entire planetary systems are propagated, such an approach is typically taken (for certain applications, the numerically propagated body must also have a particular ephemeris member object, as discussed in :ref:`tudatFeaturesPropagatorSettings`).

Available acceleration models
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
As stated above, the :literal:`createAccelerationModelsMap` function uses your environment and settings for the accelerations to automatically retrieve and put together all functions used to calculate the accelerations during each function evaluation of the numerical scheme. For reference, we first provide a bried list of available acceleration models:

- Point-mass gravity (central of third-body)
- Spherical harmonic gravity (central of third-body)
- Mutual spherical harmonic gravity (central of third-body)
- Aerodynamic acceleration
- Cannonball radiation pressure
- Panelled radiation pressure
- Solar sailing acceleration     
- Thrust acceleration
- Quasi impulsive shot acceleration
- Relativistic acceleration correction (IERS 2010 Conventions)
- Empiricical accelerations (constant, sine and cosine of true anomaly components in RSW frame)
- Tidal effect on natural satellites (Lainey et al., 2007, 2012)

Subsequently, we provide details on how to add settings for the model to the :class:`SelectedAccelerationMap`. In addition, we define the list of environment models required for their creation.

.. class:: AccelerationSettings

   Base class for setting the accelerations on a body. Settings currently available are the following:

.. method:: Point mass gravity

   Settings for a point mass gravity acceleration. No derived class of :class:`AccelerationSettings` is required, this acceleration setting are constructed by feeding :literal:`central_gravity` to the constructor. Added to :class:`SelectedAccelerationMap` as follows, for example of acceleration exerted on "Apollo" by "Earth":

   .. code-block:: cpp

       SelectedAccelerationMap accelerationSettings;
       accelerationSettings[ "Apollo" ][ "Earth" ].push_back( std::make_shared< AccelerationSettings >( central_gravity ) );

   Requires the following environment models to be defined:

   - Gravity field for body exerting acceleration (set by :class:`GravityFieldSettings`).
   - Current state of bodies undergoing and exerting acceleration, either from an Ephemeris model (set by :class:`EphemerisSettings`) or from the numerical propagation.

.. class:: SphericalHarmonicAccelerationSettings

   Settings for the accelerations as set by :class:`SphericalHarmonicsGravityFieldSettings`. Added to :class:`SelectedAccelerationMap` as follows, for example of acceleration exerted on "Apollo" by "Earth":

   .. code-block:: cpp

      SelectedAccelerationMap accelerationSettings;
      int maximumDegree = 12;
      int maximumOrder = 12;
          accelerationSettings[ "Apollo" ][ "Earth" ].push_back( std::make_shared< SphericalHarmonicAccelerationSettings >( maximumDegree, maximumOrder ) );

   where the gravity field will be expanded up to degree and order 12 in the acceleration model. Requires the following environment models to be defined:

   - Spherical harmonic gravity field for body exerting acceleration (set by :class:`SphericalHarmonicsGravityFieldSettings`).
   - Rotation model from the inertial frame to the body-fixed frame (set by :class:`RotationModelSettings`).
   - Current state of bodies undergoing and exerting acceleration, either from an ephemeris model (set by :class:`EphemerisSettings`) or from the numerical propagation.

.. note:: The spherical harmonic acceleration up to degree N and order M includes the point-mass gravity acceleration (which is the degree and order 0 term).

.. class:: MutualSphericalHarmonicAccelerationSettings

   This model is typically only used for detailed propagation of planetary systems. It is added to :class:`SelectedAccelerationMap` as follows, for example of acceleration exerted on "Io" by "Jupiter":

   .. code-block:: cpp

      SelectedAccelerationMap accelerationSettings;
      int maximumDegreeOfIo = 12;
      int maximumOrderOfIo = 12;
      int maximumDegreeOfJupiter = 4;
      int maximumOrderOfJupiter = 4;
      accelerationSettings[ "Io" ][ "Jupiter" ].push_back( std::make_shared< MutualSphericalHarmonicAccelerationSettings >( 
          maximumDegreeOfJupiter, maximumOrderOfJupiter, maximumDegreeOfIo, maximumOrderOfIo ) );

   where the gravity fields of Io and Jupiter will be expanded up to degree and order 12 and 4, respectively, in the acceleration model. Requires the following environment models to be defined:

   - Spherical harmonic gravity field for body exerting acceleration and body undergoing acceleration set by :class:`SphericalHarmonicsGravityFieldSettings`).
   - Rotation model from the inertial frame to the body-fixed frame and body undergoing acceleration (set by :class:`RotationModelSettings`).
   - Current state of bodies undergoing and exerting acceleration, either from an Ephemeris model (set by :class:`EphemerisSettings`) or from the numerical propagation.

   For the case where a third-body mutual spherical harmonic acceleration (e.g. Ganymede on Io when propagating w.r.t. Jupiter), additional parameters have to be provided that denote the expansion degree/order of the central body, so:

   .. code-block:: cpp

      SelectedAccelerationMap accelerationSettings;
      int maximumDegreeOfIo = 12;
      int maximumOrderOfIo = 12;
      int maximumDegreeOfGanymede = 4;
      int maximumOrderOfGanymede = 4;
      int maximumDegreeOfJupiter = 4;
      int maximumOrderOfJupiter = 4;
      accelerationSettings[ "Io" ][ "Jupiter" ].push_back( std::make_shared< MutualSphericalHarmonicAccelerationSettings >( 
          maximumDegreeOfJupiter, maximumOrderOfJupiter, maximumDegreeOfGanymede, maximumOrderOfGanymede, maximumDegreeOfIo, maximumOrderOfIo ) );

   where Jupiter now takes the role of central body, instead of body exerting the acceleration.

.. method:: Aerodynamic acceleration

   No derived class of :class:`AccelerationSettings` required, accessed by feeding :literal:`aerodynamic` to the constructor. Added to :class:`SelectedAccelerationMap` as follows, for example of acceleration exerted on "Apollo" by "Earth" (e.g. atmosphere model belonging to Earth):

   .. code-block:: cpp

      SelectedAccelerationMap accelerationSettings;
      accelerationSettings[ "Apollo" ][ "Earth" ].push_back( std::make_shared< AccelerationSettings >( aerodynamic ) );

   Requires the following environment models to be defined:

   - Atmosphere model for body exerting acceleration (set by :class:`AtmosphereSettings`).
   - Shape model for body exerting acceleration (set by :class:`BodyShapeSettings`).
   - Aerodynamic coefficient interface for body undergoing acceleration (set by :class:`AerodynamicCoefficientSettings`). NOTE: In the case that the aerodynamic coefficients are defined as a function of the vehicle orientation (e.g. angle of attack and sideslip angle), these angles can be manually or automatically defined. 
   - Mass model for body undergoing acceleration.
   - Current state of body undergoing acceleration and body with atmosphere.

   .. warning:: Defining settings for a vehicle's orientation, which may influence your aerodynamic force, is done after creating the acceleration models, as discused here.

.. method:: Cannonball radiation pressure

   No derived class of :class:`AccelerationSettings` required, accessed by feeding :literal:`cannon_ball_radiation_pressure` to the constructor. Added to :class:`SelectedAccelerationMap` as follows, for example of acceleration exerted on "Apollo" by "Sun":

   .. code-block:: cpp

      SelectedAccelerationMap accelerationSettings;
      accelerationSettings[ "Apollo" ][ "Sun" ].push_back( std::make_shared< AccelerationSettings >( cannon_ball_radiation_pressure ) );

   Requires the following environment models to be defined:

   - Radiation pressure model for body undergoing acceleration (from source equal to body exerting acceleration) (set by :class:`RadiationPressureInterfaceSettings`).
   - Current state of body undergoing and body emitting radiation

.. method:: Panelled radiation pressure

   No derived class of :class:`AccelerationSettings` required, accessed by feeding :literal:`panelled_radiation_pressure_acceleration` to the constructor. Added to :class:`SelectedAccelerationMap` as follows, for example of acceleration exerted on "Apollo" by "Sun":

   .. code-block:: cpp

      SelectedAccelerationMap accelerationSettings;
      accelerationSettings[ "Apollo" ][ "Sun" ].push_back( std::make_shared< AccelerationSettings >( panelled_radiation_pressure_acceleration ) );

   Requires the following environment models to be defined:

   - Panelled radiation pressure model for body undergoing acceleration (from source equal to body exerting acceleration) (set by :class:`PanelledRadiationPressureInterfaceSettings`).
   - Current state of body undergoing and body emitting radiation

.. method:: Solar sailing acceleration

   No derived class of :class:`AccelerationSettings` required, accessed by feeding :literal:`solar_sail_acceleration` to the constructor. Added to :class:`SelectedAccelerationMap` as follows, for example of acceleration exerted on "Apollo" by "Sun":

   .. code-block:: cpp

      SelectedAccelerationMap accelerationSettings;
      accelerationSettings[ "Apollo" ][ "Sun" ].push_back( std::make_shared< AccelerationSettings >( solar_sail_acceleration ) );

   Requires the following environment models to be defined:

   - Solar sailing radiation pressure model for body undergoing acceleration (from source equal to body exerting acceleration) (set by :class:`SolarSailRadiationInterfaceSettings`).
   - Current state of body undergoing and body emitting radiation


.. class:: ThrustAccelerationSettings

   Used to define the resulting accerelations of a thrust force, requiring:

   - Mass of body undergoing acceleration.
   - Settings for both the direction and magnitude of the thrust force (set by :class:`ThrustMagnitudeSettings`). These models may in turn have additional environmental dependencies. 
   
Setting up a thrust acceleration is discussed in more detail on the page :ref:`tudatFeaturesThrustModels`.

.. class:: QuasiImpulsiveShotsAccelerationSettings

   Used to define the resulting acceleration of a quasi-impulsive shot, requiring:

   - Mass of the body undergoing acceleration.
   - Settings for the characteristics of the quasi-impulsive shots (total duration, rise time, associated deltaVs), as well as the times at which they are applied.

   .. code-block:: cpp

     SelectedAccelerationMap accelerationSettings;
     std::vector< double > thrustMidTimes = { 1.0 * 3600.0, 2.0 * 3600.0, 3.0 * 3600.0 };
     std::vector< Eigen::Vector3d > deltaVValues = { 1.0E-3 * ( Eigen::Vector3d( ) << 0.3, -2.5, 3.4 ).finished( ),
       1.0E-3 * ( Eigen::Vector3d( ) << 2.0, 5.9, -0.5 ).finished( ),
       1.0E-3 * ( Eigen::Vector3d( ) << -1.6, 4.4, -5.8 ).finished( ) };
     double totalManeuverTime = 90.0;
     double maneuverRiseTime = 15.0;

     accelerationSettings[ "Vehicle" ][ "Vehicle" ] = std::make_shared< QuasiImpulsiveShotsAccelerationSettings >( 
       thrustMidTimes, deltaVValues, totalManeuverTime, maneuverRiseTime );
    

where the input variables represent:

    - Midtimes of the quasi-impulsive shots (assumed to be the time at which an ideal impulsive shot would have been applied).
    - DeltaVs (threee-dimensional vectors) associated with the quasi-impulsive shots.
    - Total duration of the quasi-impulsive shots (same value for each of them).
    - Rise time, i.e. time required to reach the peak acceleration (same value for each impulsive shot).

.. class:: RelativisticAccelerationCorrectionSettings

   A first-order (in :math:`1/c^{2}`) correction to the acceleration due to the influence of relativity. It implements the model of Chapter 10, Section 3 of the IERS 2010 Conventions. It requires a specific derived class of :class:`AccelerationSettings`. Added to :class:`SelectedAccelerationMap` as follows, for example that includes all three contributions (Schwarzschild, Lense-Thirring and de Sitter)
   
   .. code-block:: cpp

    SelectedAccelerationMap accelerationSettings;
    bool calculateSchwarzschildCorrection = true;
    bool calculateLenseThirringCorrection = true;
    bool calculateDeSitterCorrection = true;
    std::string primaryBody = "Sun";
    const Eigen::Vector3d marsAngularMomentum = ...
    accelerationSettings[ "Orbiter" ][ "Mars" ] = std::make_shared< RelativisticAccelerationCorrectionSettings >( 
       calculateSchwarzschildCorrection, calculateLenseThirringCorrection,  calculateDeSitterCorrection, primaryBody,
       centralBodyAngularMomentum )

Here, the 'primary body' for a planetary orbiter should always be set as the Sun (only relevant for de Sitter correction). The angular momentum vector of the orbited body is only relevant for Lense-Thirring correction.
    
.. class:: EmpiricalAccelerationSettings
    
   A constant/once-per-orbit acceleration, expressed in the RSW frame, for which the mangitude is determined empirically (typically during an orbit determination process). The acceleration components are defined according to Montenbruck and Gill (2000), with a total of 9 components: a constant, sine and cosine term (with true anomaly as argument) for each of the three independent directions of the RSW frame. The settings object (for a vehicle called "Orbiter" around Mars) is created as:

   .. code-block:: cpp
   
      SelectedAccelerationMap accelerationSettings;
      Eigen::Vector3d constantAcceleration = ( Eigen::Vector3d( ) << 0.4, -0.1, 0.05 ).finished( );
      Eigen::Vector3d sineAcceleration = ( Eigen::Vector3d( ) << 0.0, 0.02, 0.0 ).finished( );
      Eigen::Vector3d cosineAcceleration = ( Eigen::Vector3d( ) << -0.01, 0.0, 0.0 ).finished( );
      accelerationSettings[ "Orbiter" ][ "Mars" ] = std::make_shared< EmpiricalAccelerationSettings >( 
         constantAcceleration, sineAcceleration, cosineAcceleration );

Where the three input variables represent:
       
    - Vector containing the constant terms of the accelerations in the R, S and W directions.
    - Vector containing the sine terms of the accelerations in the R, S and W directions.
    - Vector containing the cosine terms of the accelerations in the R, S and W directions.

    
.. _tudatFeaturesFrameworkAccelerationsMassRateModelSetup:

.. class:: DirectTidalDissipationAccelerationSettings
    
  The direct of tidal effects in a satellite system, applied directly as an acceleration (as opposed to a modification of spherical harmonic coefficients). The model is based on Lainey et al. (2007,2012). It can compute either the acceleration due to tides, and in particular tidal dissipation, on a planetary satellites. The accelertion can compute either the effect of tide raised on the satellite by the planet, or on the planet by the satellite. The satellite is assumed to be tidally locked to the planet.
   
   .. code-block:: cpp

      double loveNumber = 0.1;
      double timeLag = 100.0;
    
      SelectedAccelerationMap accelerationSettings;
      accelerationSettings[ "Io" ][ "Jupiter" ] = std::make_shared< DirectTidalDissipationAccelerationSettings >(
         loveNumber, timeLag, false, false );      

Where the three input variables represent:       
                   
    - Value of the k2 Love number (real value) that is used.
    - Value of the tidal time lag (in seconds) that is used.
    - Boolean denoting whether the term independent of the time lag is to be computed (default true)
    - Boolean denoting whether the tide raised on the planet is to be modelled (if true), or the tide raised on the satellite (if false). Default is true.


