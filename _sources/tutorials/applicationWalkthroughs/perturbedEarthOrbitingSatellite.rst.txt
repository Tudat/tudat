.. _walkthroughsPerturbedEarthOrbitingSatellite:

Perturbed Earth-orbiting Satellite
==================================
The example described on this page is that of Asterix, a single satellite that is orbiting the Earth. The code for this tutorial is located in your Tudat Bundle at::

   tudatBundle/tudatExampleApplications/satellitePropagatorExamples/SatellitePropagatorExamples/singlePerturbedSatellitePropagator.cpp

For this example, we have the following problem statement:

*Given the position and velocity of Asterix at a certain point in time with respect to the Earth, what will its position and velocity be after a Julian day has passed? Take into account perturbations due to gravitational forces by other celestial bodies, due to radiation pressure and due to aerodynamic forces.*

.. warning:: The example described in this page assumes that the user has read the :ref:`walkthroughsUnperturbedEarthOrbitingSatellite`. This page only describes the differences with respect to such example, so please go back before proceding.

Set up the environment
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

Next, the :literal:`bodySettings` are created with a limitation on the time for which the Spice ephemeris is to be interpolated. Such ephemeris is created 300 seconds ahead the :literal:`simulationStartEpoch` and 300 seconds past the :literal:`simulationEndEpoch`:

.. code-block:: cpp

    // Create body objects.
    std::map< std::string, boost::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( bodiesToCreate, simulationStartEpoch - 300.0, simulationEndEpoch + 300.0 );

In the :ref:`walkthroughsUnperturbedEarthOrbitingSatellite` example, a :class:`ConstantEphemerisSettings` is used thus no limitation on the ephemeris are necessary. Next, the :literal:`frameOrientation` within :class:`EphemerisSettings` and the base frame orientation within :class:`RotationModelSettings` are set to the :literal:`"J2000"` frame:

.. code-block:: cpp

    for( unsigned int i = 0; i < bodiesToCreate.size( ); i++ )
    {
        bodySettings[ bodiesToCreate.at( i ) ]->ephemerisSettings->resetFrameOrientation( "J2000" );
        bodySettings[ bodiesToCreate.at( i ) ]->rotationModelSettings->resetOriginalFrame( "J2000" );
    }

    NamedBodyMap bodyMap = createBodies( bodySettings );

Create the vehicle
~~~~~~~~~~~~~~~~~~
A number of additional settings need to be linked to the vehicle when using additional perturbations. To begin with, the mass of the spacecraft needs to be defined:

.. code-block:: cpp

    // Create spacecraft object.
    bodyMap[ "Asterix" ] = boost::make_shared< simulation_setup::Body >( );
    bodyMap[ "Asterix" ]->setConstantBodyMass( 400.0 );

Since aerodynamic forces are considered as a perturbation, we need to include a number of variables typically used when simulating vehicles in atmospheric flight:

.. code-block:: cpp

    // Create aerodynamic coefficient interface settings.
    double referenceArea = 4.0;
    double aerodynamicCoefficient = 1.2;
    boost::shared_ptr< AerodynamicCoefficientSettings > aerodynamicCoefficientSettings =
            boost::make_shared< ConstantAerodynamicCoefficientSettings >(
                referenceArea, aerodynamicCoefficient * Eigen::Vector3d::UnitX( ), 1, 1 );

    // Create and set aerodynamic coefficients object
    bodyMap[ "Asterix" ]->setAerodynamicCoefficientInterface(
                createAerodynamicCoefficientInterface( aerodynamicCoefficientSettings, "Asterix" ) );

Next, a number of parameters necessary for the cannonball radiation pressure model are defined. One of the assumptions made here is that Earth acts an occulting body, meaning that when Asterix enters the Earth's shadow no radiation pressure is experienced:

.. code-block:: cpp

    // Create radiation pressure settings
    double referenceAreaRadiation = 4.0;
    double radiationPressureCoefficient = 1.2;

    std::vector< std::string > occultingBodies;
    occultingBodies.push_back( "Earth" );
    boost::shared_ptr< RadiationPressureInterfaceSettings > asterixRadiationPressureSettings =
            boost::make_shared< CannonBallRadiationPressureInterfaceSettings >(
                "Sun", referenceAreaRadiation, radiationPressureCoefficient, occultingBodies );

    // Create and set radiation pressure settings
    bodyMap[ "Asterix" ]->setRadiationPressureInterface(
                "Sun", createRadiationPressureInterface(
                    asterixRadiationPressureSettings, "Asterix", bodyMap ) );

Set up the acceleration models
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
So far we have defined the :class:`RadiationPressureInterfaceSettings` and the celestial bodies that will perturb the orbit of Asterix. In summary, the Asterix spacecraft will experience the following accelerations:

- Primary gravitational acceleration caused by Earth, according to a spherical-harmonics gravity model.
- Perturbing gravitational acceleration caused by the Sun, the Moon, Mars and Venus.
- Perturbing radiation pressure acceleration caused by the Sun.
- Perturbing aerodynamic acceleration cause by Earth.

These needs to be binded to the Asterix :class:`Body` object:

.. code-block:: cpp

    // Define propagation settings.
    std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > accelerationsOfAsterix;

    accelerationsOfAsterix[ "Earth" ].push_back( boost::make_shared< SphericalHarmonicAccelerationSettings >( 5, 5 ) );

    accelerationsOfAsterix[ "Sun" ].push_back( boost::make_shared< AccelerationSettings >( 
                                                   basic_astrodynamics::central_gravity ) );
    accelerationsOfAsterix[ "Moon" ].push_back( boost::make_shared< AccelerationSettings >(
                                                     basic_astrodynamics::central_gravity ) );
    accelerationsOfAsterix[ "Mars" ].push_back( boost::make_shared< AccelerationSettings >(
                                                     basic_astrodynamics::central_gravity ) );
    accelerationsOfAsterix[ "Venus" ].push_back( boost::make_shared< AccelerationSettings >(
                                                     basic_astrodynamics::central_gravity ) );
    accelerationsOfAsterix[ "Sun" ].push_back( boost::make_shared< AccelerationSettings >(
                                                     basic_astrodynamics::cannon_ball_radiation_pressure ) );
    accelerationsOfAsterix[ "Earth" ].push_back( boost::make_shared< AccelerationSettings >(
                                                     basic_astrodynamics::aerodynamic ) );

    accelerationMap[  "Asterix" ] = accelerationsOfAsterix;

Finally, :literal:`"Asterix"` is added to :literal:`bodiesToPropagate` while having :literal:`"Earth"` as the respective central body. This means that despite that other celestial bodies have been included, these will not be propagated.

.. code-block:: cpp

    bodiesToPropagate.push_back( "Asterix" );
    centralBodies.push_back( "Earth" );
