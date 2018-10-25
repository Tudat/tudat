.. _interplanetaryTrajectoryDesign:

MGA Trajectory Design
=====================

This tutorial will focus on the trajectory design tool that is included in Tudat. This tool allows the definition of a multiple gravity assist (MGA) trajectory using a patched conics approach.

The code for this tutorial is given on Github, and is also located in your Tudat bundle at: ::

   tudatBundle/tudatExampleApplications/satellitePropagatorExamples/SatellitePropagatorExamples/interplanetaryTrajectoryDesign.cpp

The first trajectory that will be calculated is that of Cassini using only gravity assists and no deep space maneuvers (DSM). The code starts by defining the different leg types:

.. code-block:: cpp

   // Specify required parameters
   // Specify the number of legs and type of legs.
   int numberOfLegs = 6;
   std::vector< int > legTypeVector;
   legTypeVector.resize( numberOfLegs );
   legTypeVector[ 0 ] = mga_Departure;
   legTypeVector[ 1 ] = mga_Swingby;
   legTypeVector[ 2 ] = mga_Swingby;
   legTypeVector[ 3 ] = mga_Swingby;
   legTypeVector[ 4 ] = mga_Swingby;
   legTypeVector[ 5 ] = capture;

The leg types can be found in the :literal:`legTypes` enum, located in the following file: ::

   tudatBundle/tudat/Tudat/Astrodynamics/TrajectoryDesign/trajectory.h

A trajectory always has to start with a departure leg, and has to end with a capture leg. 
The following step is creating a vector with pointers to the ephemerides of the planets that the vehicle will visit:

.. code-block:: cpp

   // Create the ephemeris vector.
   std::vector< ephemerides::EphemerisPointer >
           ephemerisVector( numberOfLegs );
   ephemerisVector[ 0 ] = std::make_shared< ephemerides::ApproximatePlanetPositions >(
               ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::earthMoonBarycenter );
   ephemerisVector[ 1 ] = std::make_shared< ephemerides::ApproximatePlanetPositions >(
               ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::venus );
   ephemerisVector[ 2 ] = std::make_shared< ephemerides::ApproximatePlanetPositions >(
               ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::venus );
   ephemerisVector[ 3 ] = std::make_shared< ephemerides::ApproximatePlanetPositions >(
               ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::earthMoonBarycenter );
   ephemerisVector[ 4 ] = std::make_shared< ephemerides::ApproximatePlanetPositions >(
               ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::jupiter );
   ephemerisVector[ 5 ] = std::make_shared< ephemerides::ApproximatePlanetPositions >(
               ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::saturn );

The variable has to be of type :literal:`std::vector< ephemerides::EphemerisPointer >`, but the type of ephemeris used is up to the user. This vector also determines the specific planets that will be visited along the trajectory, and the order of the visits.
Next, a vector of gravitational parameters are given for the planets.

.. code-block:: cpp

   // Create gravitational parameter vector
   Eigen::VectorXd gravitationalParameterVector( numberOfLegs );
   gravitationalParameterVector << 3.9860119e14, 3.24860e14, 3.24860e14, 3.9860119e14, 1.267e17, 3.79e16;

All the different legs in the trajectory need different characteristics to be defined by the user. For swingby only legs, the trajectory variables vector looks as follows:

.. code-block:: cpp

   // Create variable vector.
   Eigen::VectorXd variableVector( numberOfLegs + 1 );
   variableVector << -789.8117, 158.302027105278, 449.385873819743, 54.7489684339665,
           1024.36205846918, 4552.30796805542, 1/*dummy*/;
   variableVector *= physical_constants::JULIAN_DAY;

The first variable is the departure date (in MJD2000), after that, each number is the time of flight (in days) of the corresponding leg. The time of flight of the capture leg is not important for further calculations, thus it can be given a dummy value. These variables are all then converted to seconds (since MJD2000 for the departure time).

Finally, some extra variables that characterise the orbits of the vehicle in between planets and around them are given:

.. code-block:: cpp

   // Create departure and capture variables.
   Eigen::VectorXd semiMajorAxes( 2 ), eccentricities( 2 );
   semiMajorAxes << std::numeric_limits< double >::infinity( ), 1.0895e8 / 0.02;
   eccentricities << 0., 0.98;

   // Sun gravitational parameter
   const double sunGravitationalParameter = 1.32712428e20;

   // Create minimum pericenter radii vector
   Eigen::VectorXd minimumPericenterRadii( numberOfLegs );
   minimumPericenterRadii << 6778000.0, 6351800.0, 6351800.0, 6778000.0, 600000000.0, 600000000.0;

First the semi-major axes and eccentricities of the departure and target planet are given, then the central body gravitational parameter (the central body meaning the main gravitational influence when the vehicle is not in the sphere of influence of any planet) is given, and finally the minimum distance between a planet's centre and the vehicle during a swing-by is given.
Once all these variables are defined, they can be used as input to the trajectory class, which will the be able to calculate various properties of the trajectory:

.. code-block:: cpp

   // Create the trajectory problem.
   Trajectory Cassini1( numberOfLegs, legTypeVector, ephemerisVector,
                        gravitationalParameterVector, variableVector, sunGravitationalParameter,
                        minimumPericenterRadii, semiMajorAxes, eccentricities );

   // Vectors for the specific maneuvers and the total delta v
   std::vector< Eigen::Vector3d > positionVector;
   std::vector< double > timeVector;
   std::vector< double > deltaVVector;
   double resultingDeltaV;

   // Calculate the orbits
   Cassini1.calculateTrajectory( resultingDeltaV );
   Cassini1.maneuvers( positionVector, timeVector, deltaVVector );

In this example, the total delta V needed in the trajectory is calculated using :literal:`calculateTrajectory( resultingDeltaV )`, where :literal:`resultingDeltaV` will contain the final value, and the individual delta V contributions, and time and positions of these contributions are given by :literal:`maneuvers( positionVector, timeVector, deltaVVector )`.

The second example is the trajectory of Messenger and shows how DSM's can be added into the trajectory. The first thing that changes is the leg types that are defined:

.. code-block:: cpp

   // Specify required parameters
   // Specify the number of legs and type of legs.
   numberOfLegs = 5;
   legTypeVector.resize( numberOfLegs );
   legTypeVector[ 0 ] = mga1DsmVelocity_Departure;
   legTypeVector[ 1 ] = mga1DsmVelocity_Swingby;
   legTypeVector[ 2 ] = mga1DsmVelocity_Swingby;
   legTypeVector[ 3 ] = mga1DsmVelocity_Swingby;
   legTypeVector[ 4 ] = capture;

The velocity in :literal:`mga1DsmVelocity` stands for the way the DSM is calculated, more information on this can be found `here <https://repository.tudelft.nl/islandora/object/uuid%3A02468c77-5c64-4df8-9a24-1ed7ad9d1408?collection=education>`_. 
The other part of the code that is different is the trajectory variables that need to be given. It looks as follows:

.. code-block:: cpp

   // Add the time of flight and start epoch, which are in JD.
   variableVector << 1171.64503236 * physical_constants::JULIAN_DAY,
           399.999999715 * physical_constants::JULIAN_DAY,
           178.372255301 * physical_constants::JULIAN_DAY,
           299.223139512 * physical_constants::JULIAN_DAY,
           180.510754824 * physical_constants::JULIAN_DAY,
           1, // The capture time is irrelevant for the final leg.
           // Add the additional variables.
           0.234594654679, 1408.99421278, 0.37992647165 * 2 * 3.14159265358979,
           std::acos(  2 * 0.498004040298 - 1. ) - 3.14159265358979 / 2, // 1st leg.
           0.0964769387134, 1.35077257078, 1.80629232251 * 6.378e6, 0.0, // 2nd leg.
           0.829948744508, 1.09554368115, 3.04129845698 * 6.052e6, 0.0, // 3rd leg.
           0.317174785637, 1.34317576594, 1.10000000891 * 6.052e6, 0.0; // 4th leg.

The first part is the same as the MGA without DSMs, departure time and time of flights. However, the DSMs need to be defined by the user, which is done here. For each leg there are four variables, for the departure leg:

	1: the time of flight fraction at which the DSM is performed. 
	2: the hyperbolic excess velocity magnitude for the start. 
	3: the in-plane angle for the hyperbolic excess velocity. 
	4: the out-of-plane angle for the hyperbolic excess velocity.
  	
For the swing-by leg:

	1: the time of flight fraction at which the DSM is performed. 
	2: the rotation angle of the GA. 
	3: the pericenter radius of the GA.
   4: the :math:`\Delta V` added for the powered GA.

This is done for each leg containing a DSM. The calculation of the final values of the trajectory is done in the same manner as before, but now the :literal:`maneuvers( positionVector, timeVector, deltaVVector )` also contains the DSMs.

Application Output
~~~~~~~~~~~~~~~~~~
All the output of the trajectory is handled as follows:

.. code-block:: cpp

   // Define vectors to calculate intermediate points
   std::vector< Eigen::Vector3d > interPositionVectorMessenger;
   std::vector< double > interTimeVectorMessenger;

   // Calculate intermediate points and write to file
   std::string outputFileTraj = tudat_applications::getOutputPath( ) + "messengerTrajectory.dat";
   Messenger.intermediatePoints( 1000.0 , interPositionVectorMessenger, interTimeVectorMessenger );
   writeTrajectoryToFile( interPositionVectorMessenger, interTimeVectorMessenger, outputFileTraj );

   // Define vectors to calculate intermediate points
   std::vector< Eigen::Vector3d > manPositionVectorMessenger;
   std::vector< double > manTimeVectorMessenger;
   std::vector< double > manDeltaVVectorMessenger;

   // Calculate maneuvers and write to file
   std::string outputFileMan = tudat_applications::getOutputPath( ) + "messengerManeuvers.dat";
   Messenger.maneuvers( manPositionVectorMessenger, manTimeVectorMessenger, manDeltaVVectorMessenger );
   writeTrajectoryToFile( manPositionVectorMessenger, manTimeVectorMessenger, outputFileMan );

   // Calculate trajectories of the planets and output to file
   std::vector< Eigen::Vector3d > positionVectorEarth;
   std::vector< double > timeVectorEarth;
   std::vector< Eigen::Vector3d > positionVectorVenus;
   std::vector< double > timeVectorVenus;
   std::vector< Eigen::Vector3d > positionVectorMercury;
   std::vector< double > timeVectorMercury;

   // Earth
   returnSingleRevolutionPlanetTrajectory(
               std::make_shared< ephemerides::ApproximatePlanetPositions >(
                   ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::earthMoonBarycenter ),
               sunGravitationalParameter,
               1171.64503236,
               1000.0,
               positionVectorEarth,
               timeVectorEarth );

   // Venus
   returnSingleRevolutionPlanetTrajectory(
               std::make_shared< ephemerides::ApproximatePlanetPositions >(
                   ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::venus ),
               sunGravitationalParameter,
               1171.64503236,
               1000.0,
               positionVectorVenus,
               timeVectorVenus );

   // Mercury
   returnSingleRevolutionPlanetTrajectory(
               std::make_shared< ephemerides::ApproximatePlanetPositions >(
                   ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::mercury ),
               sunGravitationalParameter,
               1171.64503236,
               1000.0,
               positionVectorMercury,
               timeVectorMercury );

   std::string outputFilePlanetE = tudat_applications::getOutputPath(  ) + "earthTrajectory.dat";
   writeTrajectoryToFile( positionVectorEarth, timeVectorEarth, outputFilePlanetE );

   std::string outputFilePlanetV = tudat_applications::getOutputPath(  ) + "venusTrajectory.dat";
   writeTrajectoryToFile( positionVectorVenus, timeVectorVenus, outputFilePlanetV );

   std::string outputFilePlanetM = tudat_applications::getOutputPath(  ) + "mercuryTrajectory.dat";
   writeTrajectoryToFile( positionVectorMercury, timeVectorMercury, outputFilePlanetM );

These lines of code use several functions that can produce output from the calculated trajectory. First, intermediate points are calculated along the trajectory using :literal:`intermediatePoints`. The produced position and time vector are then passed to :literal:`writeTrajectoryToFile` to produce a file that can be read by an external plotting program. The same is done afterwards, but instead of producing some intermediate positions, now only the manuevers are produced. These are then also written to a file. Finally, for every planet the trajectory is outputted to a file to be able to see where the trajectory encounters the planet. The final figure is shown below.

.. figure:: images/interpTrajFig.png