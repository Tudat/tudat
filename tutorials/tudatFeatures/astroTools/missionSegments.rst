.. _tudatFeaturesMissionSegments:

Mission Segments
===========================

Lambert Targeter
~~~~~~~~~~~~~~~~
A Lambert targeter is an algorithm which can solve Lambert's problem. The input into this problem is the initial position vector, the final position vector, and the time-of-flight (TOF) for the peoblem that needs to be solved. Several Lambert Targeting algorithms, developed over the last few years, are included in Tudat. Each of these algorithm are a seperate class that is derived from one base class: :literal:`lambertTargeter`.

.. class:: lambertTargeter

This is the base class that is used to derive the different algorithms to solve Lambert's problem. Each class that is derived from this base class contains a method called :literal:`execute()`. This method is called in the constructor of each class and performs the specific algorithm. Other methods like: :literal:`getInertialVelocityAtDeparture()` can then be called to get the results of the algorithm. 

There are several algorithms available in Tudat. Which algorithms are available are discussed below, however, the specifics of these algorithms are not explained here. 

.. class:: lambertTargeterIzzo

This class uses Izzo's algorithm to solve the Lambert problem. This class can be called as follows:

.. code-block:: cpp
   
      LambertTargeterIzzo( const Eigen::Vector3d& aCartesianPositionAtDeparture,
                         const Eigen::Vector3d& aCartesianPositionAtArrival,
                         const double aTimeOfFlight,
                         const double aGravitationalParameter,
                         const bool isRetrograde = false,
                         const double convergenceTolerance = 1e-9,
                         const int maximumNumberOfIterations = 50 )


where:

	- :literal:`aCartesianPositionAtDeparture`:
  
	   The cartesian position of the body at departure, given in meters.

	- :literal:`aCartesianPositionAtArrival`:

	   The cartesian position of the body at arrival, given in meters.

	- :literal:`aTimeOfFlight`:

	   The time-of-flight between departure and arrival, given in seconds.

	- :literal:`aGravitationalParameter`:

	   The gravitational parameter of the main body.

	- :literal:`isRetrograde`:

	   Boolean variable that determines if the orbital motion is retrograde or not. 

	- :literal:`convergenceTolerance`:

	   :literal:`double` that gives the tolerance of the root finding algorithm.

	- :literal:`maximumNumberOfIterations`:

	   :literal:`int` that gives the maximum number of iterations of the root finding algorithm.

.. class:: lambertTargeterGooding

This class uses Gooding's algorithm to solve the Lambert problem. This class can be called as follows:

.. code-block:: cpp
   
      LambertTargeterGooding( const Eigen::Vector3d& aCartesianPositionAtDeparture,
                            const Eigen::Vector3d& aCartesianPositionAtArrival,
                            const double aTimeOfFlight,
                            const double aGravitationalParameter,
                            root_finders::RootFinderPointer aRootFinder = 
                            root_finders::RootFinderPointer( ) );


where:

	- :literal:`aCartesianPositionAtDeparture`:
  
	   The cartesian position of the body at departure, given in meters.

	- :literal:`aCartesianPositionAtArrival`:

	   The cartesian position of the body at arrival, given in meters.

	- :literal:`aTimeOfFlight`:

	   The time-of-flight between departure and arrival, given in seconds.

	- :literal:`aGravitationalParameter`:

	   The gravitational parameter of the main body.

	- :literal:`aRootFinder`:

           A pointer of type: :literal:`root_finders::RootFinderPointer` that points to a root finding algorithm.
	   
.. class:: ZeroRevolutionlambertTargeterIzzo

This class uses Izzo's zero revolution algorithm to solve the Lambert problem. This class can be called as follows:

.. code-block:: cpp
   
      ZeroRevolutionLambertTargeterIzzo( const Eigen::Vector3d& aCartesianPositionAtDeparture,
                                       const Eigen::Vector3d& aCartesianPositionAtArrival,
                                       const double aTimeOfFlight,
                                       const double aGravitationalParameter,
                                       const bool aIsRetrograde = false,
                                       const double aConvergenceTolerance = 1.0e-9,
                                       const int aMaximumNumberOfIterations = 50 )


where:

	where:

	- :literal:`aCartesianPositionAtDeparture`:
  
	   The cartesian position of the body at departure, given in meters.

	- :literal:`aCartesianPositionAtArrival`:

	   The cartesian position of the body at arrival, given in meters.

	- :literal:`aTimeOfFlight`:

	   The time-of-flight between departure and arrival, given in seconds.

	- :literal:`aGravitationalParameter`:

	   The gravitational parameter of the main body.

	- :literal:`isRetrograde`:

	   Boolean variable that determines if the orbital motion is retrograde or not. 

	- :literal:`convergenceTolerance`:

	   :literal:`double` that gives the tolerance of the root finding algorithm.

	- :literal:`maximumNumberOfIterations`:

	   :literal:`int` that gives the maximum number of iterations of the root finding algorithm.

.. class:: MultiRevolutionlambertTargeterIzzo

This class uses Izzo's multi-revolution algorithm to solve the Lambert problem. This class is not directly derived from the :literal:`lambertTargeter` class, but it is derived from the :literal:`ZeroRevolutionlambertTargeterIzzo` class. This class can be called as follows:

.. code-block:: cpp
   
      MultiRevolutionLambertTargeterIzzo( const Eigen::Vector3d& aCartesianPositionAtDeparture,
                                        const Eigen::Vector3d& aCartesianPositionAtArrival,
                                        const double aTimeOfFlight,
                                        const double aGravitationalParameter,
                                        const int aNumberOfRevolutions = 0,
                                        const bool aIsRightBranch = false,
                                        const bool aIsRetrograde = false,
                                        const double aConvergenceTolerance = 1e-9,
                                        const int aMaximumNumberOfIterations = 50 )


where:

	where:

	- :literal:`aCartesianPositionAtDeparture`:
  
	   The cartesian position of the body at departure, given in meters.

	- :literal:`aCartesianPositionAtArrival`:

	   The cartesian position of the body at arrival, given in meters.

	- :literal:`aTimeOfFlight`:

	   The time-of-flight between departure and arrival, given in seconds.

	- :literal:`aGravitationalParameter`:

	   The gravitational parameter of the main body.

        - :literal:`aNumberOfRevolutions`:

	   Required number of revolutions in the problem solution. 

        - :literal:`aIsRightBranch`:

	   A boolean flag to indicate whether the right or left branch (corresponding to a low or high energy transfer arc) should be used in the solution . 

	- :literal:`aIsRetrograde`:

	   Boolean variable that determines if the orbital motion is retrograde or not. 

	- :literal:`convergenceTolerance`:

	   :literal:`double` that gives the tolerance of the root finding algorithm.

	- :literal:`maximumNumberOfIterations`:

	   :literal:`int` that gives the maximum number of iterations of the root finding algorithm.
	   

	
