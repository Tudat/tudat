.. _walkthroughsPerturbedEarthOrbitingSatellite:

Perturbed Earth-orbiting Satellite
==================================
The example described on this page is that of Asterix, a single satellite that is orbiting the Earth. The code for this tutorial is given on Github, and is also located in your Tudat bundle at::

   tudatBundle/tudatExampleApplications/satellitePropagatorExamples/SatellitePropagatorExamples/singlePerturbedSatellitePropagator.cpp

For this example, we have the following problem statement:

    *Given the position and velocity of Asterix at a certain point in time with respect to the Earth, what will its position and velocity be after a Julian day has passed? Take into account perturbations due to gravitational forces by other celestial bodies, due to radiation pressure and due to aerodynamic forces.*

.. warning:: The example described in this page assumes that the user has read the :ref:`walkthroughsUnperturbedEarthOrbitingSatellite`. This page only describes the differences with respect such example, so please go back before proceeding.

Set Up the Environment
~~~~~~~~~~~~~~~~~~~~~~
Although this example is largely similar to the :ref:`walkthroughsUnperturbedEarthOrbitingSatellite` example, there are a number of key differences. To begin with, the environment creation includes a number of additional bodies which will provide the gravitational perturbation that affect the Asterix satellite. Such bodies are added as follows:

.. code-block:: cpp

    // Define body settings for simulation.
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Sun" );
    bodiesToCreate.push_back( "Earth" );
    bodiesToCreate.push_back( "Moon" );
    bodiesToCreate.push_back( "Mars" );
    bodiesToCreate.push_back( "Venus" );

Next, the :literal:`bodySettings` are created with a limitation on the time for which the :ref:`Spice <tudatFeaturesSpice>` ephemeris is to be interpolated. Such ephemeris is created 300 seconds ahead of the :literal:`simulationStartEpoch` and 300 seconds past the :literal:`simulationEndEpoch`:

.. code-block:: cpp

    // Create body objects.
    std::map< std::string, std::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( bodiesToCreate, simulationStartEpoch - 300.0, simulationEndEpoch + 300.0 );

In the :ref:`walkthroughsUnperturbedEarthOrbitingSatellite` example, a :class:`ConstantEphemerisSettings` is used. Thus, no limitation on the ephemeris is necessary. However, in this example additional celestial bodies are taken into account and therefore we need to provide settings for the ephemeris and rotation. The :literal:`frameOrientation` within :class:`EphemerisSettings` and the base frame orientation within :class:`RotationModelSettings` are set to the :literal:`"J2000"` frame for all bodies:

.. code-block:: cpp

    for( unsigned int i = 0; i < bodiesToCreate.size( ); i++ )
    {
        bodySettings[ bodiesToCreate.at( i ) ]->ephemerisSettings->resetFrameOrientation( "J2000" );
        bodySettings[ bodiesToCreate.at( i ) ]->rotationModelSettings->resetOriginalFrame( "J2000" );
    }
    NamedBodyMap bodyMap = createBodies( bodySettings );

Create the Vehicle
~~~~~~~~~~~~~~~~~~
A number of additional settings need to be linked to the vehicle when using additional perturbations. To begin with, the mass of the spacecraft needs to be defined:

.. code-block:: cpp

    // Create spacecraft object.
    bodyMap[ "Asterix" ] = std::make_shared< simulation_setup::Body >( );
    bodyMap[ "Asterix" ]->setConstantBodyMass( 400.0 );

We also need to set the aerodynamic coefficients of the spacecraft. These setting are stored in the :class:`AerodynamicCoefficientSettings` object. For this example, we will consider constant aerodynamic coefficients. This option is set by using the derived-class :class:`ConstantAerodynamicCoefficientSettings`. The settings for the aerodynamic coefficients are the following:

   - The reference area.
   - The aerodynamic coefficients in three directions.
   - A boolean to indicate whether the aerodynamic coefficients are defined in the aerodynamic frame (:math:`C_D`, :math:`C_S`, :math:`C_L`) or in the body frame (typically denoted as :math:`C_x`, :math:`C_y`, :math:`C_z`).
   - A boolean to define whether the aerodynamic coefficients are positive along the negative axes of the body or aerodynamic frame. 

These settings are provided in the following block of code:

.. code-block:: cpp

    // Create aerodynamic coefficient interface settings.
    double referenceArea = 4.0;
    double aerodynamicCoefficient = 1.2;
    std::shared_ptr< AerodynamicCoefficientSettings > aerodynamicCoefficientSettings =
            std::make_shared< ConstantAerodynamicCoefficientSettings >(
                referenceArea, aerodynamicCoefficient * Eigen::Vector3d::UnitX( ), 1, 1 );

    // Create and set aerodynamic coefficients object
    bodyMap[ "Asterix" ]->setAerodynamicCoefficientInterface(
                createAerodynamicCoefficientInterface( aerodynamicCoefficientSettings, "Asterix" ) );

.. tip:: Available options for :class:`AerodynamicCoefficientSettings` can be found :ref:`here <aerodynamicCoefficientOptions>`.

Next, a number of parameters necessary for the radiation pressure model are defined. This is similar to the aerodynamic coefficients as discussed above. The settings are stored in the :class:`RadiationPressureInterfaceSettings` object. This example uses a simple cannonball model. This option is set by the derived-class :class:`CannonBallRadiationPressureInterfaceSettings`. One of the assumptions made here is that Earth acts an occulting body, meaning that when Asterix enters the Earth's shadow no radiation pressure from body :literal:`"Sun"` is experienced:

.. code-block:: cpp

    // Create radiation pressure settings
    double referenceAreaRadiation = 4.0;
    double radiationPressureCoefficient = 1.2;

    std::vector< std::string > occultingBodies;
    occultingBodies.push_back( "Earth" );
    std::shared_ptr< RadiationPressureInterfaceSettings > asterixRadiationPressureSettings =
            std::make_shared< CannonBallRadiationPressureInterfaceSettings >(
                "Sun", referenceAreaRadiation, radiationPressureCoefficient, occultingBodies );

    // Create and set radiation pressure settings
    bodyMap[ "Asterix" ]->setRadiationPressureInterface(
                "Sun", createRadiationPressureInterface(
                    asterixRadiationPressureSettings, "Asterix", bodyMap ) );

.. tip:: Available options for :class:`RadiationPressureInterfaceSettings` can be found :ref:`here <radiationPressureModelOptions>`.

Set Up the Acceleration Models
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
So far we have defined the celestial bodies that will perturb the orbit of Asterix, the :class:`ArodynamicCoefficientSettings`, and  the :class:`RadiationPressureInterfaceSettings`. In summary, the Asterix spacecraft will experience the following accelerations:

   - Primary gravitational acceleration caused by Earth, according to a spherical-harmonics gravity model.
   - Perturbing gravitational acceleration caused by the Sun, the Moon, Mars and Venus.
   - Perturbing aerodynamic acceleration caused by Earth.
   - Perturbing radiation pressure acceleration caused by the Sun.

These need to be binded to the Asterix :class:`Body` object:

.. code-block:: cpp

    // Define propagation settings.
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfAsterix;
    
    accelerationsOfAsterix[ "Earth" ].push_back( std::make_shared< SphericalHarmonicAccelerationSettings >( 5, 5 ) );

    accelerationsOfAsterix[ "Sun" ].push_back( std::make_shared< AccelerationSettings >( 
                                                   basic_astrodynamics::central_gravity ) );
    accelerationsOfAsterix[ "Moon" ].push_back( std::make_shared< AccelerationSettings >(
                                                     basic_astrodynamics::central_gravity ) );
    accelerationsOfAsterix[ "Mars" ].push_back( std::make_shared< AccelerationSettings >(
                                                     basic_astrodynamics::central_gravity ) );
    accelerationsOfAsterix[ "Venus" ].push_back( std::make_shared< AccelerationSettings >(
                                                     basic_astrodynamics::central_gravity ) );
    
    accelerationsOfAsterix[ "Sun" ].push_back( std::make_shared< AccelerationSettings >(
                                                     basic_astrodynamics::cannon_ball_radiation_pressure ) );
    
    accelerationsOfAsterix[ "Earth" ].push_back( std::make_shared< AccelerationSettings >(
                                                     basic_astrodynamics::aerodynamic ) );

    accelerationMap[ "Asterix" ] = accelerationsOfAsterix;

Note that the spherical-harmonic gravitational model is implemented with the derived-class :class:`SphericalHarmonicAccelerationSettings` with inputs the degree and order of the model. Finally, :literal:`"Asterix"` is added to :literal:`bodiesToPropagate` while having :literal:`"Earth"` as the respective central body. This means that despite the inclusion of the other celestial bodies, these will not be propagated.

.. code-block:: cpp

    bodiesToPropagate.push_back( "Asterix" );
    centralBodies.push_back( "Earth" );

Results
~~~~~~~
Both the :ref:`walkthroughsUnperturbedEarthOrbitingSatellite` and :ref:`walkthroughsPerturbedEarthOrbitingSatellite` are now discussed. The orbits of both simulations can now be compared. This results are shown below. 

.. figure:: images/perturbedAndUnperturbed.png

.. tip:: Open the figure in a new tab for more detail.