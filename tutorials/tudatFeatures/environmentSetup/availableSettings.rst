Available Settings for the Environment Models
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Here, we will provide a full list of the available properties of the :class:`BodySettings` object. Each type of environment model has one base class to define settings for the creation of the model). Often, a specific derived class is implemented for a specific environment model of a given class, in which any additional information that may be needed can be provided. For instance, when defining a gravity field model, one can simply use:

.. code-block:: cpp

    bodySettings[ "Earth" ]->gravityFieldSettings = std::make_shared< GravityFieldSettings >( central_spice ); 

if you want to use a central gravity field with the gravitational parameter taken from Spice: no information is needed except the type of gravity field model that is created. On the other hand, if you want to use a spherical harmonic gravity field, you need to specify additional parameters yourself, which is done by using the specific derived class:

.. code-block:: cpp

    bodySettings[ "Earth" ]->gravityFieldSettings = std::make_shared< SphericalHarmonicsGravityFieldSettings >( gravitationalParameter, referenceRadius, cosineCoefficients, sineCoefficients, associatedReferenceFrame ); 

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

      bodySettings[ "Earth" ]->atmosphereSettings = std::make_shared< ExponentialAtmosphereSettings >( 7.2E3, 290.0, 1.225, 287.06 ); 

   for an exponential atmosphere with a scale height of 7200 m, a constant temperature of 290 K, a density at 0 m altitude of 1.225 kg/m^3 and a specific gas constant of 287.06 J/(kg K).

   If you want to model the exponential atmosphere for Earth or Mars, you can also simply input :literal:`aerodynamics::earth` or :literal:`aerodynamics::mars` to load the default settings, which are defined in the table below. 

      ============================  ==========  ==========  =================
      Property                      Earth       Mars        Units
      ============================  ==========  ==========  =================
      Scale Height                  7.2         1.11        km
      Density at Zero Altitude      1.225       0.02        kg/m :math:`{}^3`
      Constant Temperature          246.0       215.0       K
      Specific Gas Constant         287.0       197.0       J/kg/K
      Ratio of Specific Heats       1.4         1.3         --
      ============================  ==========  ==========  =================

   References for the values above are:

      - **Earth**: Lecture notes, Rocket Motion by Prof. Ir. B.A.C. Ambrosius, November 2009
      - **Mars**: Spohn, T., Breuer, D., and Johnson, T., Eds., Encyclopedia of the Solar System, 3rd ed. Elsevier, 2014

.. class:: TabulatedAtmosphereSettings

   Due to the extensive customization available for the tabulated atmosphere, you can find the settings for this class in a separate page: :ref:`tudatTabulatedAtmosphere`.

.. class:: CustomConstantTemperatureAtmosphereSettings

   With this class, you can define your own constant temperature atmosphere, which computes the atmospheric properties based on an input function. For instance, one can link a function to the settings as such:

   .. code-block:: cpp

      // Outside main
      double customDensityFunction( const double altitude, const double longitude, const double latitude, const double time )
      {
         // Return a linear combination of the input values
         return 0.5 * altitude + 0.25 * longitude + 0.15 * latitude + 0.1 * time;
      }

      int main( )
      {
        // ...

        // Define atmosphere settings
        double constantTemperature = 250.0;
        double specificGasConstant = 300.0;
        double ratioOfSpecificHeats = 1.4;
        bodySettings[ "Earth" ]->atmosphereSettings = std::make_shared< CustomConstantTemperatureAtmosphereSettings >( &customDensityFunction, constantTemperature, specificGasConstant, ratioOfSpecificHeats ); 

        // ...
      }

   As shown in the example above, the user-defined function, in this case :literal:`customDensityFunction`, is required to have those inputs, and in that specific order. The value of pressure is computed by assuming hydrostatic equilibrium, whereas temperature, gas constant and the ratio of specific heats are assumed to be constant. 

   .. tip ::
      Note that, by using :literal:`std::bind`, you can have more inputs than the ones in :literal:`customDensityFunction`. However, keep in mind that :literal:`std::bind` only allows up to 9 inputs. 

.. method:: NRLMSISE-00 
    
   This can be used to select the NRLMSISE-00 atmosphere model. To use this model, the :literal:`USE_NRLMSISE` flag in your top-level :literal:`CMakeLists` must be set to true. No derived class of :class:`AtmosphereSettings` base class required, the model can be created by passing :literal:`nrlmsise00` as argument to base class constructor. 

   .. code-block:: cpp

      bodySettings[ "Earth" ]->atmosphereSettings = std::make_shared< AtmosphereSettings >( nrlmsise00 );  

.. class:: CustomWindModelSettings

   Custom wind model which can be used to retrieve a wind vector. This wind vector is in the body-fixed, body-centered reference frame. 

   .. code-block:: cpp
   
      bodySettings[ "Earth" ]->atmosphereSettings = std::make_shared< CustomWindModelSettings >(  windFunction )
   
   where ``windFunction`` is a ``std::function`` with inputs; altitude, longitude, latitude and time.

Ephemeris model
****************  

.. class:: Ephemeris
  
   Base class for the ephemeris. It is constructed using one of the settings classes below.

.. class:: EphemerisSettings

   Base class for the ephemeris settings. Models currently available through the :class:`BodySettings` architecture and set by their respective derived classes are:

.. class:: ApproximatePlanetPositionSettings

   Highly simplified model of ephemerides of major Solar system bodies (model described here). Both a three-dimensional, and circular coplanar approximation may be used. 

   .. code-block:: cpp

       bodySettings[ "Jupiter" ]->ephemerisSettings = std::make_shared< ApproximatePlanetPositionSettings >( ephemerides::ApproximatePlanetPositionsBase::jupiter, false ); 

   where the first constructor argument is taken from the enum in approximatePlanetPositionsBase.h, and the second argument (false) denotes that the circular coplanar approximation is not made.

.. class:: DirectSpiceEphermerisSettings

   Ephemeris retrieved directly using :ref:`tudatFeaturesSpice`.

   .. code-block:: cpp

       std::string frameOrigin = "SSB";
       std::string frameOrientation = "J2000";
       bodySettings[ "Jupiter" ]->ephemerisSettings = std::make_shared< DirectSpiceEphemerisSettings >( frameOrigin, frameOrientation ); 

   creating a barycentric (SSB) ephemeris with axes along J2000, with data directly from spice.

.. class:: InterpolatedSpiceEphemerisSettings 
      
   Using this option the state of the body is retrieved at regular intervals, and used to create an interpolator, before setting up environment. This has the advantage of only requiring calls to Spice outside of the propagation inner loop, reducing computation time. However, it has the downside of begin applicable only during a limited time interval.

   .. code-block:: cpp

       double initialTime = 0.0;
       double finalTime = 1.0E8;
       double timeStep = 3600.0;
       std::string frameOrigin = "SSB";
       std::string frameOrientation = "J2000";
       bodySettings[ "Jupiter" ]->ephemerisSettings = std::make_shared< InterpolatedSpiceEphemerisSettings >(
           initialTime, finalTime, timeStep, frameOrigin, frameOrientation ); 

   creating a barycentric (SSB) ephemeris with axes along J2000, with data retrieved from Spice at 3600 s intervals between t=0 and t=1.0E8, using a 6th order Lagrange interpolator. Settings for the interpolator (discussed here, can be added as a sixth argument if you wish to use a different interpolation method)

.. class:: TabulatedEphemerisSettings

   Ephemeris created directly by interpolating user-specified states as a function of time.

   .. code-block:: cpp

       std::map< double, Eigen::Vector6d > bodyStateHistory ...
       std::string frameOrigin = "SSB";
       std::string frameOrientation = "J2000";
       bodySettings[ "Jupiter" ]->ephemerisSettings = std::make_shared< TabulatedEphemerisSettings >(
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
       bodySettings[ "Jupiter" ]->ephemerisSettings = std::make_shared< KeplerEphemerisSettings >(
           initialStateInKeplerianElements, epochOfInitialState, centralBodyGravitationalParameter, frameOrigin, frameOrientation ); 

   creating a Kepler orbit as ephemeris using the given kepler elements and associated initial time and gravitational parameter. See :ref:`tudatFeaturesFrameStateTransformations` for more details on orbital elements in Tudat.

.. class:: ConstantEphemerisSettings

   Ephemeris modelled as being independent of time.

   .. code-block:: cpp

       Eigen::Vector6d constantCartesianState = ...
       std::string frameOrigin = "SSB";
       std::string frameOrientation = "J2000";
       bodySettings[ "Jupiter" ]->ephemerisSettings = std::make_shared< ConstantEphemerisSettings >(
           constantCartesianState, frameOrigin, frameOrientation ); 

.. method:: Multi-arc ephemeris

   An ephemeris model (for translational state) that allows the body’s state to be defined by distinct ephemeris models over different arcs. Class is implemented to support multi-arc propagation/estimation. No derived class of :class:`EphemerisSettings` base class required, the created ephemeris can be made multi-arc by using the ``resetMakeMultiArcEphemeris`` function of the :class:`EphemerisSettings` class. The resulting :class:`Ephemeris` object will then be :class:`MultiArcEphemeris` (with the same ephemeris model for each arc when created, according to the settings in the :class:`EphemerisSettings` object)

   .. code-block:: cpp

      bodySettings[ "Earth" ]->ephemerisSettings-> resetMakeMultiArcEphemeris( true );   

.. class:: CustomEphemerisSettings

   Allows user to provide arbitrary function as ephemeris model. 

   .. code-block:: cpp

      std::shared_ptr< EphemerisSettings > customEphemerisSettings =
                   std::make_shared< CustomEphemerisSettings >(
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
       bodySettings[ "Earth" ]->gravityFieldSettings = std::make_shared< CentralGravityFieldSettings >( gravitationalParameter );

.. method:: Point-mass gravity field model from Spice

   Point-mass gravity field model, with gravitational parameter from Spice. No derived class of :class:`GravityFieldSettings` base class required, created by passing ``central_spice`` as argument to base class constructor.

   .. code-block:: cpp

       bodySettings[ "Earth" ]->gravityFieldSettings = std::make_shared< GravityFieldSettings >( central_spice ); 

.. class:: SphericalHarmonicsGravityFieldSettings

   Gravity field model as a spherical harmonic expansion. 

   .. code-block:: cpp

       double gravitationalParameter = ...
       double referenceRadius = ...
       Eigen::MatrixXd cosineCoefficients =  // NOTE: entry (i,j) denotes coefficient at degree i and order j
       Eigen::MatrixXd sineCoefficients =  // NOTE: entry (i,j) denotes coefficient at degree i and order j
       std::string associatedReferenceFrame = ...
       bodySettings[ "Earth" ]->gravityFieldSettings = std::make_shared< SphericalHarmonicsGravityFieldSettings >( gravitationalParameter, referenceRadius, cosineCoefficients, sineCoefficients, associatedReferenceFrame ); 

   The :literal:`associatedReferenceFrame` reference frame must presently be the same frame as the target frame of the body's rotation model (see below). It represents the frame to which the spherical harmonic field is fixed.

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
       bodySettings[ "Earth" ]->rotationModelSettings = std::make_shared< SimpleRotationModelSettings >( 
           originalFrame, targetFrame , initialOrientation, initialTime, rotationRate ); 

   where the rotation from originalFrame to targetFrame at initialTime is given by the quaternion initialOrientation. This is mapped to other times using the rotation rate rotationRate.

.. method:: Spice Rotation model

   Rotation model directly obtained from Spice. No derived class of :class:`RotationModelSettings` base class required, created by passing ``spice_rotation_model`` as argument to base class constructor.

   .. code-block:: cpp

       std::string originalFrame = "J2000";
       std::string targetFrame = "IAU_Earth";
       bodySettings[ "Earth" ]->rotationModelSettings = std::make_shared< RotationModelSettings >( spice_rotation_model, originalFrame, targetFrame ); 

.. method:: Tabulated RotationalEphemeris model

   Rotation model obtained from an interpolator, with dependent variable a ``Eigen::VectorXd`` of size 7. Currently the settings interface is not yet implemented but the functionality is implemented in :class:`TabulatedRotationalEphemeris`. The tabulated rotational ephemeris can be implemented as follows:

   .. code-block:: cpp

      // Create tabulated rotational model
      std::shared_ptr< TabulatedRotationalEphemeris< double, double > > tabulatedEphemeris =
              std::make_shared< TabulatedRotationalEphemeris<  double, double > >( rotationInterpolator );

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
       bodySettings[ "Earth" ]->shapeModelSettings = std::make_shared< SphericalBodyShapeSettings >( bodyRadius ); 

.. method:: Perfect sphere

   Model defining a body shape as a perfect sphere, with the sphere radius retrieved from Spice. No derived class of :class:`BodyShapeSettings` base class required, created by passing ``spherical_spice`` as argument to base class constructor.

   .. code-block:: cpp

       double bodyRadius = 6378.0E3;
       bodySettings[ "Earth" ]->shapeModelSettings = std::make_shared< BodyShapeSettings >( spherical_spice ); 

.. class:: OblateSphericalBodyShapeSettings  

   Model defining a body shape as a flattened sphere, with the equatorial radius and flattening provided by the user. 

   .. code-block:: cpp

       double bodyRadius = 6378.0E3;
       double bodyFlattening = 1.0 / 300.0;
       bodySettings[ "Earth" ]->shapeModelSettings = std::make_shared< OblateSphericalBodyShapeSettings >( bodyRadius, bodyFlattening ); 

.. _radiationPressureModelOptions:

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
       bodySettings[ "TestVehicle" ]->radiationPressureSettings[ sourceBody ] = std::make_shared< CannonBallRadiationPressureInterfaceSettings >(
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
       bodySettings[ "TestVehicle" ]->aerodynamicCoefficientSettings = std::make_shared< ConstantAerodynamicCoefficientSettings >( 
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