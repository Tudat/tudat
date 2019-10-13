.. _tudatFeaturesLowThrustTrajectory:

Low-Thrust Trajectories
=======================

This code allows to design a low-thrust trajectory, starting at a given departure state to reach a specified final state, in a given time-of-flight. The user can define the departure and arrival states as corresponding to the positions and velocities of existing celestial bodies, or they can be arbitrarily defined. 

The current implementation of low-thrust trajectories does not include gravity assists, so that they are referred to as low-thrust trajectory "legs". Various low-thrust trajectory models have been implemented:
	
	- **Direct methods:** 
	They rely on discretisation to turn the problem into a parameter optimisation problem. One of the most commonly used direct method is the so-called Sims-Flanagan method, developed in (ADD REFERENCE). Direct methods are easy to implement, but are computationally demanding as they require optimisation of a large set of free parameters, and they are found hard to converge without a good initial guess.
	

	- **Hybrid methods:** 
	Hybrid methods offer an interesting alternative to direct and indirect methods. Indirect methods are also a common approach to deal with low-thrust trajectory problems. They are based on the optimal control theory (costates) and require an analytical solution to the two-boundary problem, to ensure that the boundary and optimality conditions are fulfilled. The inherent complexity of the analytical derivations often leads to the need for simplifying the problem, which limits the accuracy of the solution. Modifying the dynamical model also implies to re-derive the analytical solution. Hybrid methods aim at combining the advantages of both direct and indirect methods, and limiting their drawbacks. The hybrid method developped by (ADD REFERENCE) and further improved by (ADD REFERENCE) has been implemented in TUDAT. 

	
	- **Shape-based methods:** 
	Shaping methods assume a certain shape for the low-thrust trajectory, and tune the parameters defining the trajectory shape so that the boundary conditions can be fulfilled. Their main advantage lies in their analytical formulation, which makes them extremely computationnally effective when compared to direct, indirect, or even hybrid methods. This allows for exploring a larger design space, reducing the chances to converge to a local optimum (which is one of the risks of using direct/indirect/hybrid methods), and thus to provide a good preliminary trajectory design. On the other hand, their accuracy is limited, due to the assumed trajectory shape and to the simplified dynamical model required to solve the problem analytically.


The low-thrust trajectory methods currenly available are:

	- Shape-based methods
		- Hodographic shaping
		- Spherical shaping
	- Sims-Flanagan method
	- Hybrid method


Low-Thrust trajectory leg
~~~~~~~~~~~~~~~~~~~~~~~~~

.. class:: LowThrustLeg

Base class for low-thrust trajectory leg. This class is defined as follows:

.. code-block:: cpp

	LowThrustLeg( const Eigen::Vector6d& stateAtDeparture,
                  const Eigen::Vector6d& stateAtArrival,
                  const double timeOfFlight,
                  simulation_setup::NamedBodyMap& bodyMap,
                  const std::string bodyToPropagate,
                  const std::string centralBody )
				  
where the input parameters are the following ones:
	
		- :literal:`stateAtDeparture`
			State of the spacecraft at departure (one of the boundary conditions to be fulfilled).
			
		- :literal:`stateAtArrival`
			State of the spacecraft at arrival (one of the boundary conditions to be fulfilled).
			
		- :literal:`timeOfFlight`
			Time of flight specified for the targeted low-thrust trajectory (one of the boundary conditions to be fulfilled).
			
		- :literal:`bodyMap`
			Map of pointers to body objects involved in the low-thrust trajectory.
			
		- :literal:`bodyToPropagate`
			Name of the spacecraft to be propagated.
			
		- :literal:`centralBody`
			Name of the central body of the low-thrust trajectory.

Each class derived from this base class contains the following methods:

		- :literal:`convertTimeToIndependentVariable`
			Converts current time to the corresponding independent variable value, depending on what kind of low-thrust trajectory model is used.
		
		- :literal:`convertIndependentVariableToTime`
			Converts current independent variable value to corresponding time, depending on what kind of low-thrust trajectory model is used.
		
		- :literal:`computeCurrentThrustAccelerationDirection`
			Returns direction of the current low-thrust acceleration vector.
		
		- :literal:`computeCurrentThrustAccelerationMagnitude`
			Returns magnitude of the current acceleration vector.
			
		- :literal:`computeCurrentStateVector`
			Returns current state vector along the low-thrust trajectory.
		
		- :literal:`getTrajectory`
			Returns trajectory map (by reference), filled with state history of the spacecraft at a given set of epochs, provided in an input vector.
			
		- :literal:`computeCurrentMass`
			Returns current mass of the spacecraft.
		
		- :literal:`getMassProfile`
			Returns map (by reference), filled with mass history of the spacecraft, at a given set of epochs provided in an input vector.
			
		- :literal:`getLowThrustAccelerationModel`
			Returns the thrust acceleration model corresponding to the designed trajectory.
		
		- :literal:`retrieveLowThrustAccelerationMap`
			Returns accelerations map corresponding to the low-thrust trajectory. This includes thrust acceleration and gravitational acceleration exerted by the central body of the trajectory.
		
		- :literal:`computeCurrentThrust`
			Returns current thrust vector.
		
		- :literal:`getThrustProfile`
			Returns map (by reference) filled with thrust history along the low-thrust trajectory for a set of epochs provided in an input vector.
		
		- :literal:`computeCurrentThrustAcceleration`
			Returns current thrust acceleration vector.
		
		- :literal:`getThrustAccelerationProfile`
			Returns map (by reference) filled with thrust acceleration history along the low-thrust trajectory for a set of epochs provided in an input vector.
		
		- :literal:`computeDeltaV`
			Returns deltaV associated with the low-thrust trajectory.
		
		- :literal:`computeSemiAnalyticalAndFullPropagation`
			Computes the analytical or semi-analytical low-thrust trajectory, and propagates the associated fully perturbed problem. The propagation starts at half of the time-of-flight and is performed backward until departure, and forward until arrival. Once the numerical propagation is over, the semi-analytical results are computed for the same set of epochs, to make direct comparison possible and assess the quality of the semi-analytical method. Both full propagation and semi-analytical results maps, as well as the dependent variables history maps are returned by reference.
		
		- :literal:`createLowThrustPropagatorSettings`
			Returns pair of appropriate propagator settings for the backward and forward propagations of the fully perturbed problem, to be used as input for the computeSemiAnalyticalAndFullPropagation method described above. This function returns a multi-type propagator settings, including settings to propagate the translational state and the mass of the spacecraft.
		
		- :literal:`createLowThrustTranslationalStatePropagatorSettings`
			Returns pair of appropriate translational state propagator settings for the backward and forward propagations of the fully perturbed problem, to be used as input for the computeSemiAnalyticalAndFullPropagation method described above. 


Setting up a low-thrust trajectory
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A :literal:`LowThrustTrajectoryObject` can be created using the settings class :literal:`LowThrustLegSettings`, using the function :literal:`createLowThrustLeg`. 

.. code-block:: cpp

	createLowThrustLeg(
		const std::shared_ptr< LowThrustLegSettings >& lowThrustLegSettings,
	        const Eigen::Vector6d& stateAtDeparture,
	        const Eigen::Vector6d& stateAtArrival,
	        const double& timeOfFlight,
	        simulation_setup::NamedBodyMap& bodyMap,
	        const std::string& bodyToPropagate,
	        const std::string& centralBody )

Except for the :literal:`LowThrustLegSettings` object, this function takes the departure and arrival states, as well as the required time-of-flight, the names of the spacecraft and of the central body of the trajectory, and the body map defining the trajectory environment. Using this :literal:`createLowThrustLeg` function allows the user for switching easily from one trajectory type to another, while still addressing the same design problem.

.. class:: LowThrustLegSettings

This is the base class to create :literal:`LowThrustTrajectoryObject`. The low-thrust trajectory is constructed from the settings classes derived from this base class.

.. class:: HodographicShapingLegSettings

This class defines the settings to construct a :literal:`HodographicShaping` object, which will provide a preliminary, hodographically shaped trajectory design. For more details about the hodographic shaping method, the reader is referred to (ADD REFERENCE). The definition of :literal:`HodographicShapingLegSettings` requires the user to specify the base functions to be used for the trajectory shaping, and the values of the free parameters,  if any.

.. code-block:: cpp
	
	HodographicShapingLegSettings(
            const int numberOfRevolutions,
            const double centralBodyGravitationalParameter,
            std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > >& radialVelocityFunctionComponents,
            std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > >& normalVelocityFunctionComponents,
            std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > >& axialVelocityFunctionComponents,
            const Eigen::VectorXd freeCoefficientsRadialVelocityFunction, 
            const Eigen::VectorXd freeCoefficientsNormalVelocityFunction,
            const Eigen::VectorXd freeCoefficientsAxialVelocityFunction )

.. class:: SphericalShapingLegSettings

This is the settings class for :literal:`SphericalShaping` object, which provides a spherically shaped trajectory (spherical shaping is described in more details in ADD REFERENCE). This shaping method only has one parameter which is not directly inferred from the satisfaction of departure and arrival boundary conditions. The value of this parameter is tuned until the targeted time-of-flight can be achieved. The :literal:`SphericalShapingLegSettings` settings class thus requires to specify the initial value for this free parameter, along with a :literal:`rootFinderSettings` object, to be used to find the free parameter value which will match the required time-of-flight.

.. code-block :: cpp

	SphericalShapingLegSettings(
            const int numberOfRevolutions,
            const double centralBodyGravitationalParameter,
            const double initialValueFreeCoefficient,
            const std::shared_ptr< root_finders::RootFinderSettings >& rootFinderSettings,
            const std::pair< double, double > boundsFreeCoefficient = std::make_pair( TUDAT_NAN, TUDAT_NAN ) )

.. class:: SimsFlanaganLegSettings

This is the settings class for :literal:`SimsFlanaganLeg` object. Among other inputs, it requires the user to provide the number of segments into which the trajectory is to be subdivided, according to the Sims-Flanagan parametrisation method (more details here (ADD LINK) ). An :literal:`OptimisationSettings` object is also to be provided to solve the Sims-Flanagan parametrised optimisation problem.	

.. code-block:: cpp

	SimsFlanaganLegSettings(
            const double maximumThrust,
            std::function< double( const double ) > specificImpulseFunction,
            const int numberOfSegments,
            const std::string centralBody,
            std::shared_ptr< transfer_trajectories::OptimisationSettings > optimisationSettings )

.. class:: HybridMethodLegSettings

This class defines the settings to construct a :literal:`HybridMethodLeg` object (more details about the hybrid method can be found in ADD REFERENCE). Among the different inputs of this settings class, an :literal:`OptimisationSettings` object must be provided to define the way the inherent hybrid method optimisation problem will be tackled.

.. code-block:: cpp
	
	HybridMethodLegSettings(
            const double maximumThrust,
            const double specificImpulse,
            const std::string centralBody,
            std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
            std::shared_ptr< transfer_trajectories::OptimisationSettings > optimisationSettings )
	

Optimising a low-thrust trajectory
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A pre-defined optimisation problem has been implemented, to allow the user to address a simple trajectory design optimisation problem, with limited coding effort. It tries to identify the best time-of-flight and/or departure date in order to minimize the deltaV required by the trajectory. 

The optimisation problem is defined in the class :literal:`TrajectoryOptimisationProblem`, which is implemented so that it is compatible with the PAGMO library.

.. class:: TrajectoryOptimisationProblem

This is the low-thrust trajectory optimisation class, defined as follows:

.. code-block:: cpp

	TrajectoryOptimisationProblem(
            simulation_setup::NamedBodyMap bodyMap,
            const std::string bodyToPropagate,
            const std::string centralBody,
            std::function< Eigen::Vector6d( const double ) > departureStateFunction,
            std::function< Eigen::Vector6d( const double ) > arrivalStateFunction,
            std::pair< double, double > departureTimeBounds,
            std::pair< double, double > timeOfFlightBounds,
            const std::shared_ptr< transfer_trajectories::LowThrustLegSettings >& lowThrustLegSettings )

The input parameters are the following ones:

	- :literal:`bodyMap`
		Map of pointers to :literal:`Body` objects defining the trajectory environment.

	- :literal:`bodyToPropagate`
		Name of the spacecraft to be propagated.

	- :literal:`centralBody`
		Name of the central body of the trajectory.

	- :literal:`departureStateFunction`
		Function returning the state vector at departure, as a function of the departure date.

	- :literal:`arrivalStateFunction`
		Function returning the state vector at arrival, as a function of the arrival date (defined as departure date + time-of-flight)

	- :literal:`departureTimeBounds`
		:literal:`pair` object containing the lower and upper bounds for the departure date of the trajectory.

	- :literal:`timeOfFlightBounds`
		:literal:`pair` object containing the lower and upper bounds for the time-of-flight of the trajectory.

	- :literal:`LowThrustLegSettings`
		Settings for the low-thrust trajectory to be designed.

The :literal:`fitness` function creates the relevant :literal:`LowThrustTrajectoryLeg` object out of the provided low-thrust leg settings. It then calculates the corresponding trajectory, and returns the associated deltaV.

The :literal:`get_bounds` function simply returns the time-of-flight and departure bounds which are provided as inputs of the :literal:`TrajectoryOptimisationProblem` constructor.


