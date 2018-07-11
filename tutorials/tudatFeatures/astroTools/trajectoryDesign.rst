.. _tudatFeaturesTrajectoryDesign:

Trajectory Design
===========================

This code allows the user to build up a trajectory for an interplanetary spacecraft using gravity assists (GA) and deep space maneuvers (DSM). The code is based off the work done `here <https://repository.tudelft.nl/islandora/object/uuid%3A02468c77-5c64-4df8-9a24-1ed7ad9d1408?collection=education>`_, users who wish to have more detail on the subject and/or would like to see several applications of this code are referred to this `thesis <https://repository.tudelft.nl/islandora/object/uuid%3A02468c77-5c64-4df8-9a24-1ed7ad9d1408?collection=education>`_.

The user can build a trajectory by defining the departure planet, the sequence of planets visited, the leg types of the mission, the defining variables of those legs, and the final capture planet. This code works especially well with the pagmo optimization code to generate optimal trajectories. 

The trajectories calculated here are based on several assumptions. They are as follows:

        - Patched-conics is used in the complete calculation of the trajectory. No forces except for the gravity of the central body and the thrust applied by the vehicle are considered. 

	- The planetary sequence is fixed beforehand.

	- Only one revolution per transfer is considered.

	- The transfer direction are all counter-clockwise.

	- If a DSM is applied, it can only be done once every transfer. 

If the trajectory design code is used, these assumptions need to be taken into account.


Trajectory Models
~~~~~~~~~~~~~~~~~ 
Each trajectory is built up from different legs. Each leg can be either a departure leg, swing-by leg, or a capture leg. Each trajectory has to stat with a departure leg and end with a capture leg. The departure and swing-by legscan be divided into 3 different types. These 3 different types are: Multiple Gravity Assist (MGA), Multiple Gravity Assist with 1 Deep Space Maneuver per leg using Velocity Formulation (MGA-1DSM-VF), and Multiple Gravity Assist with 1 Deep Space Maneuver per leg using Position Formulation (MGA-1DSM-PF). The differences between the VM and PM can be found `here <https://repository.tudelft.nl/islandora/object/uuid%3A02468c77-5c64-4df8-9a24-1ed7ad9d1408?collection=education>`_. Each combination of leg and type is discussed below  . 

	- :literal:`Departure leg MGA`:
  
	   This leg starts at a parking orbit around the departure planet and ends when the sphere of influence of the target planet is reached. No DSM is performed during the transfer.

	- :literal:`Departure leg MGA-1DSM-PF`:
  
	   This leg starts at a parking orbit around the departure planet and ends when the sphere of influence of the target planet is reached. One DSM is performed, and calculated using the position formulation.

	- :literal:`Departure leg MGA-1DSM-VF`:
  
	   This leg starts at a parking orbit around the departure planet and ends when the sphere of influence of the target planet is reached. One DSM is performed, and calculated using the velocity formulation.

	- :literal:`Swing-by leg MGA`:

	   A swing-by leg starts at the beginning of the swing-by maneuver (at the edge of the sphere of influence of the swing-by planet), and ends at the sphere of influence of the target planet. After the GA is performed at the swing-by planet, no additional DSM is performed.

	- :literal:`Swing-by leg MGA-1DSM-PF`:

	   A swing-by leg starts at the beginning of the swing-by maneuver (at the edge of the sphere of influence of the swing-by planet), and ends at the sphere of influence of the target planet. After the GA is performed at the swing-by planet, an additional DSM is performed. The location, magnitude, and direction are defined using the position formulation.

	- :literal:`Swing-by leg MGA-1DSM-VF`:

	   A swing-by leg starts at the beginning of the swing-by maneuver (at the edge of the sphere of influence of the swing-by planet), and ends at the sphere of influence of the target planet. After the GA is performed at the swing-by planet, an additional DSM is performed. The location, magnitude, and direction are defined using the velocity formulation.

	- :literal:`Capture leg`:

	   Capture legs are used at the end of a trajectory and consists of entering the sphere of influence of the capture planet and can last until the parking orbit around the capture planet is
	   left again.


Each leg has the option to return the :math:`\Delta` V of the leg (:literal:`CalculateLeg`), the location and size of all the maneuvers (:literal:`Maneuvers`), and the position of the spacecraft at specified intervals during the trajectory (:literal:`IntermediatePoints`). The code works by concatenating these legs and thereby forming a trajectory, for which several properties (time of flight, :math:`\Delta` V, etc.) can be calculated. The next section will explain how a user can build the interplanetary trajectory using the :literal:`Trajectory` class.
	
Trajectory Calculation
~~~~~~~~~~~~~~~~~~~~~~
The :class:`trajectory` class takes the settings for the full trajectory, and calculates the total :math:`\Delta` V, the time-of-flight, the location of the vehicle over time, and the time and location of the different maneuvers. The class is defined as follows:

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

	   A vector containing pointers to ephemeris objects of the planets that are visited (in order). See :ref:`ephemerisModel` for more information.

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
	   	- the time of flight fraction at which the DSM is performed. 
		- the hyperbolic excess velocity magnitude for the start. 
		- the in-plane angle for the hyperbolic excess velocity. 
		- the out-of-plane angle for the hyperbolic excess velocity.  	
	   For the swing-by leg:
	   	- the time of flight fraction at which the DSM is performed. 
		- the rotation angle of the GA. 
		- the pericenter radius of the GA. 
		- the :math:`\Delta` V added for the powered GA.

	- :literal:`MGA-1DSM-PF`:

	   For the departure and swing-by legs:
	   	- the time of flight fraction at which the DSM is performed. 
		- the dimesnionless radius of the DSM (position of the DSM wrt the central body divided by the departure planets position wrt the central body. 
		- the in-plane angle for DSM. 
		- the out-of-plane angle for the DSM.

The trajectory variable thus is structured as follows: (departure time, time of flights for all the legs, additional variables for each leg). 

The trajectory class contains several functions that can be used to examine the complete trajectory, they are listed below:

	- :literal:`void Trajectory::intermediatePoints( )`
	  
	  this function returns the position of the vehicle at specific times, which are defined by the maximum time step. The vectors containing the intermediate points can then be used for the :literal:`writeTrajectoryToFile` function.

	- :literal:`void Trajectory::calculateTrajectory( )`

	  this function returns the total :math:`\Delta` V of the trajectory.

	- :literal:`void Trajectory::maneuvers( )`

	  the maneuvers function returns the position and time of all the manuevers executed during the trajectory. These can be passed to the :literal:`writeManeuversToFile` function.

	- :literal:`void Trajectory::planetaryOrbits( )`

	  returns a single revolution of the planets that are encountered during the trajectory, to be plotted from the file generated by the :literal:`writeTrajectoryToFile` function.

	
  	
	   

