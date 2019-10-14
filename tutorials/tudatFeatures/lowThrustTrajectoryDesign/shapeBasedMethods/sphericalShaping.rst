.. _tudatFeaturesSphericalShaping:

Spherical shaping
=================

The spherical shape-based method shapes the position of the spacecraft using spherical coordinates, and uses the azimuth angle as independent variable. This shaping method assumes an inverse function for the radial distance, and a polynomial function for the elevation angle. It only requires two shaping functions, as one of the spherical coordinates is used as independent variable. As opposed to hodographic shaping, the current implementation of the spherical shaping method does not let the choice of base functions free. The ones proposed in Novak, 2012 (Methods and tools for preliminary low thrust mission analysis), where this shaping method has been initially developed, are used directly. The weighting coefficients of those base functions are calculated so that the boundary conditions of the low-thrust trajectory are fulfilled. The spherical shaping method has no free parameter (on the contrary to the hodographic shaping method), so that the shaped trajectory is entirely defined by its boundary conditions and time-of-flight and cannot be optimised in terms of deltaV budget.

One of the main drawbacks of this shaping method is that it assumes the thrust is always tangential to the orbit. Consequently, there is no guarantee this shape-based method always leads to a feasible trajectory, depending on the boundary conditions and required time-of-flight. This is especially the case for trajectories implying large inclination changes.

Moreover, because this method shapes the position of the spacecraft and not its velocity, matching the required time of flight is not inherently ensured by the implementation of the shape-based method itself (as it is the case for hodographic shaping for instance). An iterative process is required to match the expected time-of-flight. The coefficient associated to the second base function of the inverse polynomial radial distance function (referred to as the "free coefficient" is the following) is thus not determined by the departure and arrival boundary conditions, but its value is tuned until the expected time-of-flight is matched. 

Setting up a spherically shaped trajectory
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A spherically shaped trajectory is defined with the :literal:`SphericalShaping` class. It is a derived class from :literal:`ShapeBasedMethodLeg`. A corresponding settings class (:literal:`SphericalShapingLegSettings`) has been implemented to facilitate the creation of :literal:`SphericalShaping` objects (see :ref:`tudatFeaturesSetUpLowThrustTrajectory` for more details on the :literal:`SphericalShapingLegSettings`).

.. class:: SphericalShaping

The SphericalShaping class is defined as follows:

.. code-block:: cpp
   
      SphericalShaping(Eigen::Vector6d initialState,
                     Eigen::Vector6d finalState,
                     double requiredTimeOfFlight,
                     int numberOfRevolutions,
                     simulation_setup::NamedBodyMap& bodyMap,
                     const std::string bodyToPropagate,
                     const std::string centralBody,
                     double initialValueFreeCoefficient,
                     std::shared_ptr< root_finders::RootFinderSettings >& rootFinderSettings,
                     const double lowerBoundFreeCoefficient = TUDAT_NAN,
                     const double upperBoundFreeCoefficient = TUDAT_NAN,
                     std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings
                        = std::shared_ptr< numerical_integrators::IntegratorSettings< > >( ) )


where the inputs are:

	- :literal:`initialState`
		State of the spacecraft at departure.

	- :literal:`finalState`
		State of the spacecraft at arrival.

	- :literal:`requiredTimeOfFlight`
		Time of flight required for the shape-based trajectory.

	- :literal:`numberOfRevolutions`
		Required number of revolutions before the spacecraft reaches its final state.

	- :literal:`bodyMap`
		Map of pointers to :literal:`Body` objects involved in the low-thrust trajectory.

	- :literal:`bodyToPropagate`	
		Name of the spacecraft to be propagated.

	- :literal:`centralBody`
		Name of the central body of the low-thrust trajectory.

	- :literal:`initialValueFreeCoefficient`
		Initial value for the free coefficient, to be tuned until the required time-of-flight is matched.

	- :literal:`rootFinderSettings`
		Settings that define the root finder algorithm, used to find the proper free coefficient value for which the required time-of-flight would be achieved.

	- :literal:`lowerBoundFreeCoefficient`
		Lower bound for the free coefficient search.

	- :literal:`upperBoundFreeCoefficient`
		Upper bound for the free coefficient search.

	- :literal:`integratorSettings`
		Integrator settings (empty by default), used to propagate the spacecraft mass or thrust profiles, or to numerically propagate the fully perturbed trajectory (as a means to assess the quality of the analytical shaped-based preliminary design).



.. note::

	The creation of a spherically shaped trajectory requires much less inputs from the user than hodographic shaping does. The shaping functions are indeed fixed in spherical shaping, while they are user-defined in hodogaphic shaping. Morever, spherical shaping does not allow for free parameter, as opposed to hodographic shaping. This makes the definition of a spherically shaped trajectory more straightforward (no need to set up the base functions) from the user's point of view. The downside of it is that spherical is less flexible than hodographic shaping, and cannot be optimised to get a lower deltaV. 

	   

	
