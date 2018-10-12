.. _linkEndSetup:

Setting up the link ends
========================

In Tudat, an observation model is referred to a number of link ends. For a one-way range model, for instance, a receiver and transmitter are required, both of which are termed a 'link end'. An :literal:`enum` termed :literal:`LinkEndType` is available that lists all the possible kinds of link ends. A full list of observation models, as well as the link ends that they require, is given on the page on :ref:`observationModelSetup`.

.. _groundStationCreation:

Ground Station Creation
~~~~~~~~~~~~~~~~~~~~~~~

Often, you will need to define the positions of ground stations on celestial bodies to/from which observations are made. Presently, the position of a ground station is fixed in the body-fixed frame of the body on which it is located, but modifications to allow a time-varying position (due to tides, continental drift, *etc.*) are planned.

With the present setup, you must provide only the name and position (as an :literal:`Eigen::Vector3d`) of the ground station to create it. The position may be defined in several coordinate systems, defined by an enum of :literal:`PositionElementTypes`:

* :literal:`cartesian_position` For this type, the three entries of the position are simply the *x*, *y* and *z* components in the body-fixed frame.
* :literal:`spherical_position` For this type, the three entries of the position are the distance from the body's center of mass, and the geocentric latitude and longitude on the body (in that order).
* :literal:`geodetic_position` For this type, the three entries of the position are the altitude (w.r.t. the body's shape model), the geodetic latitude and the longitude (in that order). 

Now, creating a ground station is done by:

.. code-block:: cpp

    std::shared_ptr< Body > earth = ...;// Define object for Earth

    // Define name of station
    std::string stationName = "Graz"; 
    
    // Position of ground station
    Eigen::Vector3d stationPosition;
    stationPosition << 4194511.7, 1162789.7, 4647362.5;
    PositionElementTypes positionType = cartesian_position;
     
    createGroundStation( earth, stationName, stationPosition, positionType );

Which will create a ground station named :literal:`"Graz"` on the body :literal:`earth` at the given position.

Creating a Set of Link Ends
~~~~~~~~~~~~~~~~~~~~~~~~~~~

A set of link ends used for a given observable are stored in a :literal:`LinkEnds` type, which is a :literal:`typedef` for :literal:`std::map< LinkEndType, std::pair< std::string, std::string > >`. As you can see, the map value is a pair of strings. The first entry of the string is the body on which the link end is placed, the second entry the reference point on this body (typically the ground station). In the case where the second entry is empty (as is often the case for spacecraft), the body's center of mass is used. An example of defining link ends is given below:

.. code-block:: cpp

    LinkEnds oneWayLinkEnds;
    oneWayLinkEnds[ transmitter ] = std::make_pair( "Earth", "Graz" );
    oneWayLinkEnds[ receiver ] = std::make_pair( "LRO", "" );
    
This defines a link for which the ground station termed Graz on the boy called Earth acts as transmitter, and the body called LRO is used as the receiver (in this case placed at the body's center of mass).

An example of link-ends for a two-way link from Graz to LRO and back to Graz is:

.. code-block:: cpp

    LinkEnds twoWayLinkEnds;
    twoWayLinkEnds[ transmitter ] = std::make_pair( "Earth", "Graz" );
    twoWayLinkEnds[ reflector1 ] = std::make_pair( "LRO", "" );
    twoWayLinkEnds[ receiver ] = std::make_pair( "Earth", "Graz" );

Where the Graz station now acts as both transmitter and receiver. Similarly, the receiver may be different from the transmitter (in what is typically called a three-way observable in Deep Space tracking ), so:

.. code-block:: cpp

    LinkEnds twoWayLinkEnds;
    twoWayLinkEnds[ transmitter ] = std::make_pair( "Earth", "Graz" );
    twoWayLinkEnds[ reflector1 ] = std::make_pair( "LRO", "" );
    twoWayLinkEnds[ receiver ] = std::make_pair( "Earth", "Matera" );
    
where the signal is transmitter by Graz station, retransmitter or reflected by LRO, and then received by the Matera station.