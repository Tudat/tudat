.. _tudatFeaturesTrajectoryDesign:

Trajectory Design
===========================

This code allows the user to build up a trajectory for an interplanetary spacecraft using gravity assists (GA) and deep space maneuvers (DSM). The code is based off the work done `here <https://repository.tudelft.nl/islandora/object/uuid%3A02468c77-5c64-4df8-9a24-1ed7ad9d1408?collection=education>`_, users who wish to have more detail on the subject and/or would like to see several applications of this code are referred to this `thesis <https://repository.tudelft.nl/islandora/object/uuid%3A02468c77-5c64-4df8-9a24-1ed7ad9d1408?collection=education>`_.

The user can build a trajectory by defining the departure planet, the sequence of planets visited, the leg types of the mission, the defining variables of those legs, and the final capture planet. This code works especially well with the pagmo optimization code to generate optimal trajectories. 


Trajectory Models
~~~~~~~~~~~~~~~~~ 
A trajectory can be of three different types: Multiple Gravity Assist (MGA) trajectory, Multiple Gravity Assist with 1 Deep Space Maneuver per leg using Velocity Formulation (MGA-1DSM-VF) trajectory, and Multiple Gravity Assist with 1 Deep Space Maneuver per leg using Position Formulation (MGA-1DSM-PF) trajectory. The differences between the VM and PM can be found `here <https://repository.tudelft.nl/islandora/object/uuid%3A02468c77-5c64-4df8-9a24-1ed7ad9d1408?collection=education>`_. Each trajectory model is build up from multiple legs, which can be selected by the user. There are three different leg types that can be selected:

	- :literal:`Departure leg`:
  
	   This leg starts at a parking orbit around the departure planet and ends when the sphere of influence of the target planet is reached.

	- :literal:`Swing-by leg`:

	   A swing-by leg starts at the beginning of the swing-by maneuver (at the edge of the sphere of influence of the swing-by planet), and ends at the sphere of influence of the target planet.

	- :literal:`Capture leg`:

	   Capture legs are used at the end of a trajectory and consists of entering the sphere of influence of the capture planet and can last until the parking orbit around the capture planet is
	   left again.


Each leg has the option to return the :math:`\Delta` V of the leg (:literal:`CalculateLeg`), the location and size of all the maneuvers (:literal:`Maneuvers`), and the position of the spacecraft at specified intervals during the trajectory (:literal:`IntermediatePoints`). The code works by concatenating these legs and thereby forming a trajectory, for which several properties (time of flight, :math:`\Delta` V, etc.) can be calculated. The next section will explain how a user can build the interplanetary trajectory using the :literal:`Trajectory` class.
	
Trajectory Calculation
~~~~~~~~~~~~~~~~~~~~~~
The interplanetary trajectory is made using the :literal:`Trajectory` class:

.. code-block:: cpp
   
      Trajectory( const double numberOfLegs,
                  const std::vector< int >& legTypeVector,
                  std::vector< ephemerides::EphemerisPointer >& 
			ephemerisVector,
                  const Eigen::VectorXd& gravitationalParameterVector,
                  const Eigen::VectorXd& trajectoryVariableVector,
		  const double centralBodyGravitationalParameter,
                  const Eigen::VectorXd& minimumPericenterRadiiVector,
                  const Eigen::VectorXd& semiMajorAxesVector,
		  const Eigen::VectorXd& eccentricityVector, )


where the inputs are: 

	- :literal:`numberOfLegs`:
  
	   The amount of legs in the trajectory.

	- :literal:`legTypeVector`:

	   A vector containing all the types of the different legs (in order), which are defined by the :literal:`legTypes` enum located in the :literal:`trajectory.h` file.

	- :literal:`ephemerisVector`:

	   A vector containing pointers to ephemeris objects of the planets that are visited (in order).

	- :literal:`gravitationalParameterVector`:

	   A list of the gravitational parameters of the planets that are defined in :literal:`ephemerisVector`.

	- :literal:`trajectoryVariableVector`:

	   Vector of trajectory variables for which the structure will be explained below.

	- :literal:`centralBodyGravitationalParameter`:

	   The gravitational parameter of the central body of the trajectory (for most cases the Sun, but could be a planet for trajectories visiting the moons of a planet).

	- :literal:`minimumPericenterRadiiVector`:

	   A vector containing the closest the spacecraft can approach the planet during a swing-by.

	- :literal:`semiMajorAxesVector`:

	   A vector containing the semi-major axes of the parking orbits around the departure and capture planet.

	- :literal:`eccentricityVector`:

	   A vector containing the eccentricities of the parking orbits around the departure and capture planet.

The :literal:`trajectoryVariableVector` contains the variables that define the different legs. The first entry should be the departure time, then the time of flights of each of the legs should be entered (in order) and finally the additional variables can be entered for each leg. The additional variables are different for each trajectory model and are defined below:

	- :literal:`MGA`:
  
	   There are no additional variables needed for the MGA trajectory model.

	- :literal:`MGA-1DSM-VF`:

	   For the departure leg:
	   1: the time of flight fraction at which the DSM is performed. 2: the hyperbolic excess velocity magnitude for the start. 3: the in-plane angle for the hyperbolic excess velocity. 4:
	   the out-of-plane angle for the hyperbolic excess velocity.
  	
	   For the swing-by leg:
	   1: the time of flight fraction at which the DSM is performed. 2: the rotation angle of the GA. 3: the pericenter radius of the GA. 4: the :math:`\Delta` V added for the powered GA.

	- :literal:`MGA-1DSM-PF`:

	   For the departure and swing-by legs:
	   1: the time of flight fraction at which the DSM is performed. 2: the dimesnionless radius of the DSM (position of the DSM wrt the central body divided by the departure planets position 
	   wrt the central body. 3: the in-plane angle for DSM. 4: the out-of-plane angle for the DSM.

The trajectory variable thus is structured as follows: (departure time, time of flights for all the legs, additional variables for each leg). 

The trajectory class has functions to calculate the complete :math:`\Delta` V needed for the trajectory, called: :literal:`calculateTrajectory( double& totalDeltaV )`, a function that returns the location and :math:`\Delta` V for all the maneuvers, called: :literal:`maneuvers( positionVector, timeVector, deltaVVector)`, and a function to return intermediate points along a trajectory, called: :literal:`intermediatePoints( maxTimeStep, positionVector, timeVector )`.
  	
	   

