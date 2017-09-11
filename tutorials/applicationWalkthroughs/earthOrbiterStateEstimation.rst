.. _earthOrbiterStateEstimation:

Orbit Determination and Parameter Estimation
============================================

In all previous tutorials, we were only concerned with the propagation of orbits, and the analysis of the numerical results. In this tutorial, we will show how to simulate tracking observables, and use these observations to estimate the state of a spacecraft, as well as a variety of physical parameters of the environment.

 The code for this tutorial is given on Github, and is also located in your tudat bundle at::

    tudatBundle/tudatExampleApplications/satellitePropagatorExamples/SatellitePropagatorExamples/earthOrbiterStateEstimation.cpp

As you can see, setting up the environment is done in a similar manner as the previous tutorials: a set of celestial bodies is created, as well as a vehicle, which is endowed with radiation pressure and aerodynamic properties. 

The first modification is that we change the Earth rotation model

.. code-block:: cpp

   bodySettings[ "Earth" ]->rotationModelSettings = boost::make_shared< SimpleRotationModelSettings >(
                    "ECLIPJ2000", "IAU_Earth",
                    spice_interface::computeRotationQuaternionBetweenFrames(
                        "ECLIPJ2000", "IAU_Earth", initialEphemerisTime ),
                    initialEphemerisTime, 2.0 * mathematical_constants::PI /
                    ( physical_constants::JULIAN_DAY + 40.0 * 60.0

Which we do so that we can estimate the rotational properties of the Earth, as they are described by this model (fixed rotation axis and rotation rate).

Secondly, after creating the bodies, we now create a number of ground stations on our body Earth. Ground stations serve as reference points from which observations can be performed. Their body-fixed state is, in this case, defined in geodetic coordinates. See THIS PAGE for more details:

.. code-block:: cpp

    // Create ground stations from geodetic positions.
    std::vector< std::string > groundStationNames;
    groundStationNames.push_back( "Station1" );
    groundStationNames.push_back( "Station2" );
    groundStationNames.push_back( "Station3" );

    createGroundStation( bodyMap.at( "Earth" ), "Station1",
                         ( Eigen::Vector3d( ) << 0.0, 0.35, 0.0 ).finished( ), geodetic_position );
    createGroundStation( bodyMap.at( "Earth" ), "Station2",
                         ( Eigen::Vector3d( ) << 0.0, -0.55, 2.0 ).finished( ), geodetic_position );
    createGroundStation( bodyMap.at( "Earth" ), "Station3",
                         ( Eigen::Vector3d( ) << 0.0, 0.05, 4.0 ).finished( ), geodetic_position );
                         
The subsequent creation of acceleration, propagation and integration settings is done in the same manner as previous tutorials.

Defining Observation Settings
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
In Tudat, an observation model is referred to a number of link ends. For a one-way range model, for instance, a receiver and transmitter are required, both of which are termed a 'link end'. An :literal:`enum` termed :literal:`LinkEndType` is available that lists all the possible kinds of link ends. A full list of observation models, as well as the link ends that they require, is given on THIS PAGE. A set of link ends used for a given observable are stored in a :literal:`LinkEnds` type, which is a :literal:`typedef` for :literal:`std::map< LinkEndType, std::pair< std::string, std::string > >`. As you can see, the map value is a pair of strings. The first entry of the string is the body on which the link end is placed, the second entry the reference point on this body (typically the ground station). In the case where teh second entry is empty (as is often the case for spacecraft), the body's center of mass is used.

Here, we want to create a set of link ends that use each of the ground stations as a receiver, and the spacecraft as a transmitter, as well as vice versa. We do this by:

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

For instance, :literal:`stationReceiverLinkEnds.at( 1 )` will now denote a set of link ends where the spacecraft is the transmitter, and ground station 1 is the receiver. 

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

Where you can see that the :literal:`ObservableType` denotes a type of observation. Here, we limit ourselves to 1-way range, 1-way Doppler and angular position observables.

Now that we've defined which link ends are used for which observables, we can start adding more properties to the observation models. This is done by using the :literal:`ObservationSettings` class. This class is discussed in more detail on THIS PAGE. For this tutorial, we restrict ourselves to simple observation models (which do not require any information in addition to their type) and we do not use observation biases or light-time corrections.

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
            observationSettingsMap.insert(
                        std::make_pair( currentLinkEndsList.at( i ),
                                        boost::make_shared< ObservationSettings >(
                                            currentObservable ) );
                    //, boost::shared_ptr< LightTimeCorrectionSettings >( ),
                    //                        biasSettings ) ) );
        }
    }
    
Where we have defined a map :literal:`ObservationSettingsMap` that contains all the settings necessary to create the observation models.

Defining Estimation Settings 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~