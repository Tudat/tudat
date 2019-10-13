.. _tudatFeaturesSimsFlanagan:

Sims-Flanagan Method
====================

The Sims-Flanagan method is the most commonly used direct method. It consists in dividing the trajectory into two separate parts: a forward propagation subleg from departure to half of the time of flight, and a backward propagation subleg from arrival to half of the time of flight. Those two sublegs are themselves divided into a given number of segments, so that the low-thrust trajectory is modelled as a succession of impulsive shots, which are applied at the middle of each trajectory segment. The direction and magnitude of the impulsive shots are left as free parameters, but are constrained to ensure they do not overshoot the maximum allowable thrust magnitude. The Sims-Flanagan method thus transforms the two-boundary problem into a parametrised optimisation problem, which is defined as follows:

	- **Objective**
		Minimisation of the deltaV required by the trajectory.

	- **Constraints**
		- The state vectors obtained with respectively forward propagation from departure, and backward propagation from arrival should match at half of the time of flight.
		- The thrust magnitude shall not exceed its maximum value (provided as input by the user).

	- **Design parameters**
		So-called thrust throttles: direction and magnitude of the impulsive shots applied at the middle of each trajectory segment.
	

Sims-Flanagan optimisation problem
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Sims-Flanagan method has been implemented in the class :literal:`SimsFlanagan`, which is derived from the base class :literal:`LowThrustLeg`. 

.. class:: SimsFlanagan

This class creates the Sims-Flanagan problem and runs the associated optimisation problem. It inherits from the base class :literal:`LowThrustLeg` and is defined as follows:

.. code-block:: cpp
	
	SimsFlanagan(
            const Eigen::Vector6d& stateAtDeparture,
            const Eigen::Vector6d& stateAtArrival,
            const double maximumThrust,
            const std::function< double ( const double ) > specificImpulseFunction,
            const int numberSegments,
            const double timeOfFlight,
            simulation_setup::NamedBodyMap bodyMap,
            const std::string bodyToPropagate,
            const std::string centralBody,
            pagmo::algorithm optimisationAlgorithm,
            const int numberOfGenerations,
            const int numberOfIndividualsPerPopulation,
            const double relativeToleranceConstraints = 1.0e-6,
            std::pair< std::function< Eigen::Vector3d( const double ) >, double > initialGuessThrustModel = std::make_pair( nullptr, 0.0 ) )
			
where the input parameters are:
	
	- :literal:`stateAtDeparture`
		State of the spacecraft at departure.
		
	- :literal:`stateAtArrival`
		State of the spacecraft at arrival.
		
	- :literal:`maximumThrust`
		Maximum thrust magnitude per impulsive shot.
		
	- :literal:`specificImpulseFunction`
		Function returning the specific impulse value at a given time.
		
	- :literal:`numberSegments`
		Total number of segments into which the trajectory leg is subdivided (this includes both the forward and backward propagation legs into which the whole trajectory is divided into).
		
	- :literal:`timeOfFlight`
		Time of flight required for the trajectory.		
		
	- :literal:`bodyMap`
		Map of pointers to :literal:`Body` objects involved in the low-thrust trajectory.
		
	- :literal:`bodyToPropagate`
		Name of the spacecraft to be propagated.
		
	- :literal:`centralBody`
		Name of the central body of the low-thrust trajectory.

	- :literal:`optimisationSettings`
		Settings for the Sims-Flanagan optimisation (includes algorithm to be used, number of generations, number of individuals per population, tolerance with respect to contraints, and possibly initial guess)			
					
		
The :literal:`SimsFlanagan` class is directly derived from the base class :literal:`LowThrustLeg`, and all the methods contained in that base class are thus available from any :literal:`SimsFlanagan` object (see :ref:`tudatFeaturesLowThrustTrajectory` for more details). This includes, among others, the methods allowing the user to retrieve the trajectory, mass, thrust, and thrust acceleration history along the trajectory.
		
The constructor of the class :literal:`SimsFlanagan` automatically calls the method :literal:`performOptimisation`, which tries to find the solution to the Sims-Flanagan parametrised optimisation problem, using the optimisation settings provided by the user. This method creates an object of the class :literal:`SimsFlanaganProblem`, which defines the Sims-Flanagan optimisation problem so that it is compatible with the PAGMO library. The :literal:`performOptimisation` function then runs the optimisation algorithm and saves the identified optimum in the :literal:`simsFlanaganLeg_` private variable, which can be retrieved from the :literal:`getOptimalSimsFlanaganLeg` method.


Sims-Flanagan trajectory model
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Trying the solve the Sims-Flanagan parametrised optimisation problem, as presented above, implies considering many different parametrised trajectories with various thrust profiles, and trying to converge towards the best performing one (in terms of constraints satisfaction, and deltaV budget). Each individual parametrised trajectory whose performance is to be assessed with respect to the trajectory design optimisation objective is created as an object of the following class.

.. class:: SimsFlanaganModel

This class models the low-thrust trajectory as described by the Sims-Flanagan method: forward and backward propagation legs from departure and arrival respectively, matching at half of the time of flight, and low-thrust profile discretised as a succession of impulsive shots applied along the trajectory. This class does not solve the parametrised optimisation problem, but can simply propagate the trajectory defined by a set of user-defined impulsive shots. The class is defined as:

.. code-block::cpp
	
	SimsFlanaganLeg( const Eigen::Vector6d& stateAtDeparture,
                     const Eigen::Vector6d& stateAtArrival,
                     const double maximumThrust,
                     const std::function< double ( const double ) > specificImpulseFunction,
                     const double timeOfFlight,
                     simulation_setup::NamedBodyMap& bodyMap,
                     std::vector< Eigen::Vector3d >& throttles,
                     const std::string bodyToPropagate,
                     const std::string centralBody )
	
The input parameters of this class constructor are:
	
	- :literal:`stateAtDeparture`
		State of the spacecraft at departure.
			
	- :literal:`stateAtArrival`
		State of the spacecraft at arrival.
			
	- :literal:`maximumThrust`
		Maximum thrust magnitude per impulsive shot.
		
	- :literal:`specificImpulseFunction`
		Function returning the specific impulse value at a given time.
		
	- :literal:`timeOfFlight`
		Time of flight required for the trajectory.
		
	- :literal:`bodyMap`
		Map of pointers to :literal:`Body` objects involved in the low-thrust trajectory.
		
	- :literal:`throttles`
		Vector containing the thrust vectors (normalised with respect to the maximum thrust value) for each of the impulsive shots.
		
	- :literal:`bodyToPropagate`
		Name of the spacecraft to be propagated.
			
	- :literal:`centralBody`
		Name of the central body of the low-thrust trajectory.

.. _tudatFeaturesSimsFlanaganInitialGuessFromShaping:

Using shape-based trajectory as an initial guess
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

It can be rather difficult to reach convergence when trying to solve the parametrised Sims-Flanagan optimisation problem. This is mostly due to the fact that the constant thrust vector of each of the trajectory segment can be arbitrarily chosen, so that the parameters search space is extremely large. This is why it is recommended to use a good initial guess as a starting point for the Sims-Flanagan trajectory design method. The shaping methods which have been implemented in Tudat are good candidates to compute a preliminary trajectory design in an efficient way. As the design parameters for the Sims-Flanagan method are the thrust vectors of the trajectory segments, the initial guess must be provided as a set of n thrust vectors, constant in magnitude and directions (n being the number of segments of the Sims-Flanagan trajectory). One can use the function :literal:`getInitialGuessFunctionFromShaping` to approximate a thrust profile delivered by a shaping method by a set of n successive, constant thrust vectors.

.. code-block:: cpp

    std::function< Eigen::Vector3d( const double ) > getInitialGuessFunctionFromShaping(
        std::shared_ptr< shape_based_methods::ShapeBasedMethodLeg > shapeBasedLeg,
        const int numberSegmentsSimsFlanagan,
        const double timeOfFlight,
        std::function< double( const double ) > specificImpulseFunction,
        std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings )

The input parameters of this function are the following ones:

	- :literal:`shapeBasedLeg`
		Pointer to a :literal:`ShapeBasedObject` from which the thrust profile of the shape-based trajectory used as initial guess is to be retrieved.  

	- :literal:`numberSegmentsSimsFlanagan`
		:literal:`int` object defining the number of segments into which the trajectory is to be subdivided when using the Sims-Flanagan method.

	- :literal:`timeOfFlight`
		Expected time-of-flight for the trajectory.

	- :literal:`specificImpulseFunction`
		Function returning the specific impulse of the spacecraft as a function of time.

	- :literal:`integratorSettings`
		Integrator settings to be used to retrieve the thrust profile of the shape-based trajectory.

Below is an example of how a thrust profile derived from a shape-based trajectory is approximated to a set of succcessive, constant thrust vectors which can be used as an initial guess for the Sims-Flanagan thrust throttles. It is based on an Earth-Mars transfer, and the shaping method used to get the rough preliminary design is hodographic shaping. 

.. figure:: images/initialGuessSimsFlanagan.png





