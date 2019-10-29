.. _tudatFeaturesHybridMethod:

Hybrid Method
=============

The hybrid method implemented in TUDAT follows from the work done in the theses by Boudestijn (2014) and Jimenez-Lluva (2018). It aims at combining the advantages of direct and indirect methods, while limiting their respective drawbacks.

It makes use of the optimal control theory to reduce the set of free parameters compared to direct methods. As opposed to traditional indirect methods, the costates are not derived from analytically solving the two-point-boundary-value-problem. A linear interpolation is assumed between the initial and final values of the Modified Equinoctial Elements (MEE) costates, and those initial and final costate values are left as the only free parameters of the problem. The thrust profile along the low-thrust trajectory is defined  from the MEE costates guidance and magnitude models. 

The trajectory design problem is turned into a simplified optimisation problem, defined as follows:

	- **Objective**
		Minimisation of the deltaV required by the trajectory.
	
	- **Constraints**
		The state vector obtained at the end of the propagation should match with the targeted state at arrival.
	
	- **Design parameters**
		Initial and final values of the MEE costates.
	
	
Hybrid method optimisation problem
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The hybrid method has been implemented in the class :literal:`HybridMethod`, which inherits from the base class :literal:`LowThrustLeg`. 

.. class:: HybridMethod
	
This class defines the trajectory design problem as addressed by the hybrid method and runs the associated optimisation problem. The constructor is defined as follows:

.. code-block:: cpp
	
	HybridMethod(
		Eigen::Vector6d& stateAtDeparture,
		const Eigen::Vector6d& stateAtArrival,
		const double maximumThrust,
		const double specificImpulse,
		const double timeOfFlight,
		simulation_setup::NamedBodyMap& bodyMap,
		const std::string bodyToPropagate,
		const std::string centralBody,
		std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
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
		
	- :literal:`specificImpulse`
		Function returning the specific impulse value at a given time.
		
	- :literal:`timeOfFlight`
		Time of flight required for the trajectory.		
		
	- :literal:`bodyMap`
		Map of pointers to :literal:`Body` objects involved in the low-thrust trajectory.
		
	- :literal:`bodyToPropagate`
		Name of the spacecraft to be propagated.
		
	- :literal:`centralBody`
		Name of the central body of the low-thrust trajectory.
			
	- :literal:`integratorSettings`
		Integrator settings to propagate the spacecraft trajectory.

	- :literal:`optimisationSettings`
		Settings for the Sims-Flanagan optimisation (includes algorithm to be used, number of generations, number of individuals per population, tolerance with respect to contraints, and possibly initial guess)
			
					
The :literal:`HybridMethod` class is directly derived from the base class :literal:`LowThrustLeg`, and all the methods contained in that base class are thus available from any :literal:`HybridMethod` object (see :ref:`tudatFeaturesLowThrustTrajectory` for more details). This includes, among others, the methods allowing the user to retrieve the trajectory, mass, thrust, and thrust acceleration history along the trajectory.
		
The constructor of the class :literal:`HybridMethod` automatically calls the method :literal:`performOptimisation`, which tries to find the solution to the optimisation problem as defined by the hybrid method, using the optimisation settings provided by the user. This method creates the structure :literal:`HybridMethodProblem`, which defines the corresponding optimisation problem. It thus saves the identified optimum in the :literal:`hybridMethodLeg_` private variable, which can be retrieved from the (TO BE COMPLETED) method.


Hybrid method trajectory model
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Addressing a trajectory design problem with the hybrid method described above requires to run an optimisation algorithm to find the best trajectory, i.e. a trajectory which fulfills the departure and arrival boundary conditions, as well as the required time-of-flight, for the lowest deltaV. This requires the assessment of numerous trajectories to finally identify the best performing one. Each of the individual trajectories parsed by the optimisation algorithm is obtained out of a given set of initial and final MEE costates, and is defined as an object of the class described below. 

.. class:: HybridMethodModel

This class models the low-thrust trajectory as described by the hybrid method. It simply propagates the trajectory assuming a simplified model with thrust and central body gravitational accelerations only. The thrust acceleration is derived from the optimal control theory, using a vector of costates functions which return the value of each costate as a function of time. The costate functions are directly derived from the linear interpolation between their initial and final values, which are provided by the user. This class does not solve the optimisation problem, but it defines the low-thrust trajectory corresponding to a given set of initial and final costate values. The class is defined as:
	
.. code-block:: cpp

	HybridMethodLeg( const Eigen::Vector6d& stateAtDeparture,
                     const Eigen::Vector6d& stateAtArrival,
                     const Eigen::VectorXd& initialCoStates,
                     const Eigen::VectorXd& finalCoStates,
                     const double maximumThrust,
                     const double specificImpulse,
                     const double timeOfFlight,
                     simulation_setup::NamedBodyMap& bodyMap,
                     const std::string bodyToPropagate,
                     const std::string centralBody,
                     std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings )
	
The input parameters of this class constructor are:
	
	- :literal:`stateAtDeparture`
		State of the spacecraft at departure.
			
	- :literal:`stateAtArrival`
		State of the spacecraft at arrival.
		
	- :literal:`initialCoStates`
		Vector containing the values of each of the MEE costates at departure.
		
	- :literal:`finalCoStates`
		Vector containing the values of each of the MEE costates at arrival.
			
	- :literal:`maximumThrust`
		Maximum thrust magnitude. The thrust model used in the hybrid method based on the optimal control theory is defined as a so-called "bang-bang" thrust model: the magnitude of the thrust vector is equal to either 0 or the maximum thrust value.
		
	- :literal:`specificImpulse`
		Specific Impulse value. The current implementation of the hybrid method does not allow for time-varying specific impulse. 
		
	- :literal:`timeOfFlight`
		Time of flight required for the trajectory.
		
	- :literal:`bodyMap`
		Map of pointers to :literal:`Body` objects involved in the low-thrust trajectory.
		
	- :literal:`bodyToPropagate`
		Name of the spacecraft to be propagated.
			
	- :literal:`centralBody`
		Name of the central body of the low-thrust trajectory.
			
	- :literal:`integratorSettings`
		Integrator settings to be used to propagate the spacecraft trajectory.
		
		







