.. _earthOrbiterBasicStateEstimation:

Orbit Determination and Parameter Estimation (Basic)
====================================================

In all previous tutorials, we were only concerned with the propagation of orbits, and the analysis of the numerical results. In this tutorial, we will show how to simulate tracking observables, and use these observations to estimate the state of a spacecraft, as well as a variety of physical parameters of the environment. The estimation framework in Tudat has a broad variety of features, and we only cover the basic aspects here. A tutorial of more extensive features is given on the page :ref:`earthOrbiterStateEstimation`.

The code for this tutorial is given on Github, and is also located in your tudat bundle at::

   tudatBundle/tudatExampleApplications/satellitePropagatorExamples/SatellitePropagatorExamples/earthOrbiterBasicStateEstimation.cpp

As you can see in the code, setting up the environment is done in a similar manner as the previous tutorials: a set of celestial bodies is created, as well as a vehicle, which is endowed with radiation pressure and aerodynamic properties. 

The first modification is that we change the Earth rotation model

.. code-block:: cpp

   bodySettings[ "Earth" ]->rotationModelSettings = std::make_shared< SimpleRotationModelSettings >(
               "ECLIPJ2000", "IAU_Earth", spice_interface::computeRotationQuaternionBetweenFrames(
                   "ECLIPJ2000", "IAU_Earth", initialEphemerisTime ),
               initialEphemerisTime, 2.0 * mathematical_constants::PI / physical_constants::JULIAN_DAY );

We take this extra step, so that we can estimate the rotational properties of the Earth, as they are described by this model (fixed rotation axis and rotation rate). Details of the architecture of the parameter estimation is given on the page on :ref:`parameterArchitecture`.

Secondly, after creating the bodies, we add a number of ground stations on our body Earth. Ground stations serve as reference points from which observations can be performed. Their body-fixed state is, in this case, defined in geodetic coordinates. See the page :ref:`groundStationCreation` for more details:

.. code-block:: cpp

   // Create ground stations from geodetic positions.
   std::vector< std::string > groundStationNames;
   groundStationNames.push_back( "Station1" );
   groundStationNames.push_back( "Station2" );
   groundStationNames.push_back( "Station3" );

   createGroundStation( bodyMap.at( "Earth" ), "Station1",
                        ( Eigen::Vector3d( ) << 0.0, 1.25, 0.0 ).finished( ), geodetic_position );
   createGroundStation( bodyMap.at( "Earth" ), "Station2",
                        ( Eigen::Vector3d( ) << 0.0, -1.55, 2.0 ).finished( ), geodetic_position );
   createGroundStation( bodyMap.at( "Earth" ), "Station3",
                        ( Eigen::Vector3d( ) << 0.0, 0.8, 4.0 ).finished( ), geodetic_position );

The subsequent creation of acceleration, propagation and integration settings is done in the same manner as previous tutorials.

Defining Observation Settings
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
In Tudat, an observation model is referred to a number of link ends (see :ref:`linkEndSetup`). For a one-way range model, for instance, a receiver and transmitter are required, both of which are termed a 'link end'. An :literal:`enum` termed :literal:`LinkEndType` is available that lists all the possible kinds of link ends. A full list of observation models, as well as the link ends that they require, is given on page :ref:`observationTypes`. A set of link ends used for a given observable are stored in a :literal:`LinkEnds` type, which is a :literal:`typedef` for :literal:`std::map< LinkEndType, std::pair< std::string, std::string > >`. As you can see, the map value is a pair of strings. The first entry of the pair is the body on which the link end is placed, whereas the second entry is the reference point on this body (typically the ground station). In the case where the second entry is empty (as is often the case for spacecraft), the body's center of mass is used.

Here, we want to create a set of link ends that use each of the ground stations as a receiver, and the spacecraft as a transmitter, as well as vice-versa. We do this by:

.. code-block:: cpp

   // Create list of link ends in which station is receiver and in which station is transmitter (with spacecraft other link end).
   std::vector< LinkEnds > stationReceiverLinkEnds;
   std::vector< LinkEnds > stationTransmitterLinkEnds;

   for( unsigned int i = 0; i < groundStationNames.size( ); i++ )
   {
       LinkEnds linkEnds;
       linkEnds[ transmitter ] = std::make_pair( "Earth", groundStationNames.at( i ) );
       linkEnds[ receiver ] = std::make_pair( "Vehicle", "" );
       stationTransmitterLinkEnds.push_back( linkEnds );

       linkEnds.clear( );
       linkEnds[ receiver ] = std::make_pair( "Earth", groundStationNames.at( i ) );
       linkEnds[ transmitter ] = std::make_pair( "Vehicle", "" );
       stationReceiverLinkEnds.push_back( linkEnds );
   }

For instance, :literal:`stationReceiverLinkEnds.at( 1 )` will now denote a set of link ends where the spacecraft is the transmitter, and the ground station termed :literal:`"Station2"` is the receiver. 

Next, we need to define which link ends are to be used for which observable. We do this somewhat arbitrarily, and define:

.. code-block:: cpp

   // Define (arbitrarily) link ends to be used for 1-way range, 1-way doppler and angular position observables
   std::map< ObservableType, std::vector< LinkEnds > > linkEndsPerObservable;
   linkEndsPerObservable[ one_way_range ].push_back( stationReceiverLinkEnds[ 0 ] );
   linkEndsPerObservable[ one_way_range ].push_back( stationTransmitterLinkEnds[ 0 ] );
   linkEndsPerObservable[ one_way_range ].push_back( stationReceiverLinkEnds[ 1 ] );

   linkEndsPerObservable[ one_way_doppler ].push_back( stationReceiverLinkEnds[ 1 ] );
   linkEndsPerObservable[ one_way_doppler ].push_back( stationTransmitterLinkEnds[ 2 ] );

   linkEndsPerObservable[ angular_position ].push_back( stationReceiverLinkEnds[ 2 ] );
   linkEndsPerObservable[ angular_position ].push_back( stationTransmitterLinkEnds[ 1 ] );

Where you can see that the :literal:`ObservableType` enumeration denotes which types of observations are considered. Here, we limit ourselves to 1-way range, 1-way Doppler and angular position observables.

Now that we have defined which link ends are used for which observables, we can start adding more properties to the observation models. This is done by using the :class:`ObservationSettings` class. This class is discussed in more detail on the page :ref:`observationSettings`. For this tutorial, we restrict ourselves to simple observation models (which do not require any information in addition to their type) and we do not use observation biases or light-time corrections.

The resulting code to create settings for the observation models then becomes:

.. code-block:: cpp

   observation_models::ObservationSettingsMap observationSettingsMap;
   for( std::map< ObservableType, std::vector< LinkEnds > >::iterator linkEndIterator = linkEndsPerObservable.begin( );
        linkEndIterator != linkEndsPerObservable.end( ); linkEndIterator++ )
   {
       ObservableType currentObservable = linkEndIterator->first;

       std::vector< LinkEnds > currentLinkEndsList = linkEndIterator->second;
       for( unsigned int i = 0; i < currentLinkEndsList.size( ); i++ )
       {
           // Define settings for observable, no light-time corrections, and biases for selected 1-way range links
           observationSettingsMap.insert( std::make_pair( currentLinkEndsList.at( i ),
                                                          std::make_shared< ObservationSettings >( currentObservable ) ) );
       }
   }

Where we have defined a map :class:`ObservationSettingsMap` (a typedef for :literal:`std::multimap< LinkEnds, std::shared_ptr< ObservationSettings > >`) that contains all the settings necessary to create the observation models.

Defining Estimation Settings 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We have now defined the settings for the observation models that are to be used. Next are the settings for the parameters that are to be estimated. In this tutorial, we use only a limited set of parameters, namely:

   - The spacecraft initial state :math:`x_{0}`, where use only a single arc to estimate its dynamics.
   - A constant radiation pressure coefficient :math:`C_{r}` of the spacecraft (assuming a cannonball radiation pressure model)
   - A constant aerodynamic drag coefficient :math:`C_{D}` of the spacecraft
   - Spherical harmonic cosine coefficients at degree 2, and orders 0 to 2 (so :math:`C_{20}, C_{21}, C_{22}`)
   - Spherical harmonic sine coefficients at degree 2, and orders 1 to 2 (so :math:`S_{21}, S_{22}`)

Defining the settings for these parameters is done by:

.. code-block:: cpp

   // Define list of parameters to estimate.
   std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;
   parameterNames.push_back( std::make_shared< InitialTranslationalStateEstimatableParameterSettings< double > >(
                                 "Vehicle", systemInitialState, "Earth" ) );
   parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Vehicle", radiation_pressure_coefficient ) );
   parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Vehicle", constant_drag_coefficient ) );
   parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                                 2, 0, 2, 2, "Earth", spherical_harmonics_cosine_coefficient_block ) );
   parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                                 2, 1, 2, 2, "Earth", spherical_harmonics_sine_coefficient_block ) );

details on the set up of the parameters can be found on the page :ref:`parameterSettingCreation`. The general idea behind the settings may be familiar: they are similar to the acceleration settings. Some parameters (:math:`C_{r}` and :math:`C_{D}`) require no information in addition to the type of parameter and associated bodies and are created using the :class:`EstimatableParameterSettings` base class. The other parameters require additional information, and have a dedicated derived class.

Now, the actual objects that are used in the simulation are created by:

.. code-block:: cpp

   // Create parameters
   std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate =
           createParametersToEstimate( parameterNames, bodyMap );

   // Print identifiers and indices of parameters to terminal.
   printEstimatableParameterEntries( parametersToEstimate );

Where the second part (i.e., :literal:`printEstimatableParameterEntries`) is optional, and produces a list of the estimated parameters to your console. The output should be something like: ::

   Parameter start index, Parameter definition
   0, translational state of (Vehicle).
   6, radiation pressure coefficient of (Vehicle).
   7, constant drag coefficient of (Vehicle).
   8, cosine spherical harmonic coefficient block of (Earth), Minimum D/O: (2, 0), Maximum D/O: (2, 2). 
  11, sine spherical harmonic coefficient block of (Earth), Minimum D/O: (2, 1), Maximum D/O: (2, 2). 

which provides information on the physical meaning of the entries of the parameter vector (note that the order is not necessarilly the same as in the :literal:`parameterNames` list). Here, the initial state starts at index 0, the radiation pressure at index 6, etc.

Initializing Dynamics, Observation Models and Partial Derivatives
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

All required objects that compute the dynamics, variational equations, observation models and observation partials are created by the following line of code:

.. code-block:: cpp

   // Create orbit determination object (propagate orbit, create observation models)
   OrbitDeterminationManager< double, double > orbitDeterminationManager =
           OrbitDeterminationManager< double, double >(
               bodyMap, parametersToEstimate, observationSettingsMap,
               integratorSettings, propagatorSettings );

The :class:`OrbitDeterminationManager` object that is created will automatically propagate the dynamics (accordinng to the :literal:`integratorSettings` and :literal:`propagatorSettings`), as well as the associated variational equations (according to, in this case, the same :literal:`propagatorSettings`). Observation models are created using :literal:`observationSettingsMap`, as well as the associated models for the observation partial derivatives. More details can be found in :ref:`estimationObjectCreation`.

Simulating Observations
~~~~~~~~~~~~~~~~~~~~~~~

The tutorial is concerned with using *simulated* data to perform the estimation. Here, we discuss how to generate simulated observations. First, we start by defining the times at which we want to simulate observations:

.. code-block:: cpp

   // Define time of first observation
   double observationTimeStart = initialEphemerisTime + 1000.0;

   // Define time between two observations
   double  observationInterval = 20.0;

   // Simulate observations for 3 days
   std::vector< double > baseTimeList;
   for( unsigned int i = 0; i < 3; i++ )
   {
       // Simulate 500 observations per day (observationInterval apart)
       for( unsigned int j = 0; j < 500; j++ )
       {
           baseTimeList.push_back( observationTimeStart + static_cast< double >( i ) * 86400.0 +
                                   static_cast< double >( j ) * observationInterval );
       }
   }

So, we start simulating at :math:`t=t_{0}+1000` (with :math:`t_{0}` the start time of the simulation), and then simulate 500 simulations 20 seconds apart at the start of each of the 3 days in the simulations. This list of time is then stored in the :literal:`baseTimeList` vector.

In general, observation times will be different for each link end/observable type. Here, however, we take a simpler approach and use the same observation time for each link:

.. code-block:: cpp

   // Create measureement simulation input
   std::map< ObservableType, std::map< LinkEnds, std::pair< std::vector< double >, LinkEndType > > > measurementSimulationInput;
   for( std::map< ObservableType, std::vector< LinkEnds > >::iterator linkEndIterator = linkEndsPerObservable.begin( );
        linkEndIterator != linkEndsPerObservable.end( ); linkEndIterator++ )
   {
       // Define observable type and link ends
       ObservableType currentObservable = linkEndIterator->first;
       std::vector< LinkEnds > currentLinkEndsList = linkEndIterator->second;

       // Define observation times and reference link ends
       for( unsigned int i = 0; i < currentLinkEndsList.size( ); i++ )
       {
           measurementSimulationInput[ currentObservable ][ currentLinkEndsList.at( i ) ] =
                   std::make_pair( baseTimeList, receiver );
       }
   }

This code iterates over all observable types and link ends (which were stored in the :literal:`linkEndsPerObservable` variable), and then populate the :literal:`measurementSimulationInput` map. This map contains a list of observation times for each link ends/observable. Note that the input to the map is :literal:`std::make_pair( baseTimeList, receiver )`, not only :literal:`baseTimeList`. The :literal:`receiver` identifier denotes that the observation time is valid at reception of the signal (not at its transmission).

Simulating the observations is then done as:

.. code-block:: cpp

   // Set typedefs for POD input (observation types, observation link ends, observation values, associated times with
   // reference link ends.
   typedef Eigen::Matrix< double, Eigen::Dynamic, 1 > ObservationVectorType;
   typedef std::map< LinkEnds, std::pair< ObservationVectorType, std::pair< std::vector< double >, LinkEndType > > >
           SingleObservablePodInputType;
   typedef std::map< ObservableType, SingleObservablePodInputType > PodInputDataType;

   // Simulate observations
   PodInputDataType observationsAndTimes = simulateObservations< double, double >(
               measurementSimulationInput, orbitDeterminationManager.getObservationSimulators( ) );

In this simulation, we have completely neglected the fact that the spacecraft may not be visible from the ground station from which the observation is taken. Tudat has the capabilities to prune the observations with this (and other) checks, but this is discussed in a later tutorial (and in more detail on the page on :ref:`observationViability`). The simulated data type :literal:`observationsAndTimes` now contains a simulated observables, along with information on the associated observable type, link ends and times.
 
Performing the Estimation
~~~~~~~~~~~~~~~~~~~~~~~~~

Now that we have our simulated data and our estimation objects all ready to go, we can perform the actual simulated estimation. Since this is a simulated scenario without noise, we first need to perturb our parameter vector a bit, otherwise the postfit residuals will all be exactly zero even on the first iteration. This we do by:

.. code-block:: cpp
   
   // Perturb parameter estimate
   Eigen::Matrix< double, Eigen::Dynamic, 1 > initialParameterEstimate =
           parametersToEstimate->template getFullParameterValues< double >( );
   Eigen::Matrix< double, Eigen::Dynamic, 1 > truthParameters = initialParameterEstimate;
   Eigen::Matrix< double, Eigen::Dynamic, 1 > parameterPerturbation =
           Eigen::Matrix< double, Eigen::Dynamic, 1 >::Zero( truthParameters.rows( ) );
   parameterPerturbation.segment( 0, 3 ) = Eigen::Vector3d::Constant( 10.0 );
   parameterPerturbation.segment( 3, 3 ) = Eigen::Vector3d::Constant( 1.0E-2 );
   parameterPerturbation( 6 ) = 0.01;
   parameterPerturbation( 7 ) = 0.01;
   initialParameterEstimate += parameterPerturbation;

In which we perturb the initial position and velocity by 10 m and 0.01 m/s, respectively. Both :math:`C_{r}` and  :math:`C_{D}`, we perturb by 0.01. Note that only the parameters of the model are changed, so that the estimation should converge (to within its numerical capabilities) to the original parameter set.

We define the input to the estimation with the :class:`PodInput` class:

.. code-block:: cpp

   // Define estimation input
   std::shared_ptr< PodInput< double, double > > podInput =
           std::make_shared< PodInput< double, double > >(
               observationsAndTimes, initialParameterEstimate.rows( ),
               Eigen::MatrixXd::Zero( truthParameters.rows( ), truthParameters.rows( ) ),
               initialParameterEstimate - truthParameters );
   podInput->defineEstimationSettings( true, true, false, true );

where details on the input (and the :literal:`defineEstimationSettings` function) is given here: :ref:`estimationInput`. Additionally, since we are using different observables we must set their weights explicitly (they are all set as 1 if we don't). This we do by:

.. code-block:: cpp

   // Define observation weights (constant per observable type)
   std::map< observation_models::ObservableType, double > weightPerObservable;
   weightPerObservable[ one_way_range ] = 1.0 / ( 1.0 * 1.0 );
   weightPerObservable[ angular_position ] = 1.0 / ( 1.0E-5 * 1.0E-5 );
   weightPerObservable[ one_way_doppler ] = 1.0 / ( 1.0E-11 * 1.0E-11 );
   podInput->setConstantPerObservableWeightsMatrix( weightPerObservable );

which sets the expected range data precision at 1.0 m, the angular position data at 10 mrad, and the Doppler data at :math:`10^{-11}` (Doppler data is range-rate uncertainty non-dimensionalized by the speed of light).

The estimation is then performed by:

.. code-block:: cpp

   // Perform estimation
   std::shared_ptr< PodOutput< double > > podOutput = orbitDeterminationManager.estimateParameters(
               podInput, std::make_shared< EstimationConvergenceChecker >( 4 ) );

Where the index 4 indicates that the estimation will perform 4 iterations. The estimation should produce output similar to the following: ::

   Calculating residuals and partials 13500
   Parameter update   -6.08378    -18.9676    -7.76344 -0.00696149  -0.0186691 -0.00883959   -0.165673  -0.0139745 -1.4458e-09 9.05548e-08 3.77696e-08 1.44452e-07   8.374e-10
   Current residual: 2502.59
   Warning, tabulated ephemeris is being reset using data at different precision
   Calculating residuals and partials 13500
   Parameter update    -3.91642      8.96792     -2.23662  -0.00303872   0.00866954  -0.00116045      0.15567   0.00397531  1.44563e-09 -9.05592e-08  -3.7772e-08 -1.44462e-07 -8.40113e-10
   Current residual: 1.40276
   Warning, tabulated ephemeris is being reset using data at different precision
   Calculating residuals and partials 13500
   Parameter update 0.000194697 -0.000314532  5.82188e-05  2.08074e-07 -4.85731e-07  4.99801e-08   4.5677e-06 -2.21992e-06  2.18631e-13   4.6825e-12  2.48334e-12  9.83626e-12  2.55367e-12
   Current residual: 0.000221693
   Warning, tabulated ephemeris is being reset using data at different precision
   Calculating residuals and partials 13500
   Parameter update 5.18789e-06 -9.23697e-06 -2.86618e-06  4.95883e-09  4.12459e-09  2.61187e-09 -2.04626e-06  1.16347e-06 -5.01255e-14 -2.42908e-13 -1.31398e-13  9.26677e-14  1.81624e-13
   Current residual: 0.000496538
   Maximum number of iterations reached
   Final residual: 0.000221693

which shows the estimation progress. Clearly, no improvements are made in the final iteration, so that only 3 iterations would have been needed. The example also prints the true and formal estimation errors where, in this case, the true error is much smaller, as the observations are noise-free (and not noisy at the level presumed by the weights). The performance is in this case only limited by the numerical precision, causing the initial state to be estimated at the 0.01 mm level. The simulation also prints various quantities to files, which are then processed by the MATLAB function provided.

.. A few examples of results are given below: 