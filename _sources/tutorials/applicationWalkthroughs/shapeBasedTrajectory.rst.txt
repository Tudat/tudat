.. _walkthroughsShapeBasedTrajectory:

Shape-Based Low-Thrust Trajectories
===================================
The example described in this tutorial presents the use of shape-based methods to design low-thrust trajectories. The code for this tutorial is given on Github, and is also located in your Tudat bundle at:

   tudatBundle/tudatExampleApplications/libraryExamples/satellitePropagatorExamples/shapeBasedTrajectoryDesign.cpp

This tutorial presents how a preliminary low-thrust trajectory can be obtained using various shape-based methods. The optimisation of the trajectory will not be addressed here, but this tutorial rather focuses on how to design a feasible trajectory for an Earth-Mars transfer which respects boundary conditions at departure and arrival, as well as the specified value of the time-of-flight.

Shaping methods use a simplified, unperturbed model to compute shaped trajectories. This tutorial also presents the numerical propagation of the fully perturbed problem, and compares the results with the analytical trajectory derived from the shape-based methods.

Set up the problem
~~~~~~~~~~~~~~~~~~
The boundary conditions of the problem, required Cartesian state and time at departure and arrival (or in this case the time-of-flight), are first defined  for the targeted trajectory. In addition, the expected number of revolutions are specified. This describes the trajectory design problem that is to be tackled with shape-based methods.

.. code-block:: cpp

    int numberOfRevolutions = 1;
    double julianDateAtDeparture = 8174.5 * physical_constants::JULIAN_DAY;
    double  timeOfFlight = 580.0 * physical_constants::JULIAN_DAY;

    // Ephemeris departure body.
    ephemerides::EphemerisPointer pointerToDepartureBodyEphemeris = std::make_shared< ephemerides::ApproximatePlanetPositions>(
                ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::earthMoonBarycenter );

    // Ephemeris arrival body.
    ephemerides::EphemerisPointer pointerToArrivalBodyEphemeris = std::make_shared< ephemerides::ApproximatePlanetPositions >(
                ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::mars );

    // Retrieve cartesian state at departure and arrival.
    Eigen::Vector6d cartesianStateAtDeparture = pointerToDepartureBodyEphemeris->getCartesianState( julianDateAtDeparture );
    Eigen::Vector6d cartesianStateAtArrival = pointerToArrivalBodyEphemeris->getCartesianState( julianDateAtDeparture + timeOfFlight );

	
The environment within which the trajectory of the spacecraft is to be calculated is set up in a similar manner as in the previous tutorials: the required celestial bodies are created (from default body settings), along with a body object representing our vehicle. The ephemeris of the Sun (central body of the trajectory) is defined so that it stays at the Solar System Barycenter. 


Set up the accelerations (unperturbed and perturbed cases)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
Two different sets of accelerations are defined here:

	- The first one corresponds to the unperturbed problem: only the two accelerations already accounted for in the shape-based trajectory design (namely central body gravitational acceleration and low-thrust acceleration) are exerted on the spacecraft.

	.. code-block:: cpp

	    // Define acceleration map for the simplified problem.
	    // (empty map as central gravity and low-thrust accelerations are already included in the shape-based methods)
	    basic_astrodynamics::AccelerationMap perturbingAccelerationsMapSimplifiedProblem;


	- In the second case, additional perturbations are taken into account:

		- Gravitational attraction exerted by the Earth, Mars, and Jupiter.
		- Solar radiation pressure

	.. code-block:: cpp

	    // Define acceleration map for the fully perturbed problem.

	    // Define propagation settings.
	    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationSettingsPerturbedProblem;
	    accelerationSettingsPerturbedProblem[ "Earth" ].push_back( std::make_shared< AccelerationSettings >(
		                                             basic_astrodynamics::central_gravity ) );
	    accelerationSettingsPerturbedProblem[ "Mars" ].push_back( std::make_shared< AccelerationSettings >(
		                                             basic_astrodynamics::central_gravity ) );
	    accelerationSettingsPerturbedProblem[ "Jupiter" ].push_back( std::make_shared< AccelerationSettings >(
		                                             basic_astrodynamics::central_gravity ) );
	    accelerationSettingsPerturbedProblem[ "Sun" ].push_back( std::make_shared< AccelerationSettings >(
		                                             basic_astrodynamics::cannon_ball_radiation_pressure ) );


The definition of additional perturbations will later aim at quantifying the effects of the simplifying assumptions used in the shape-based trajectory design, when propagating the fully perturbed trajectory of the spacecraft.


.. warning::
	
	The set of accelerations defined above are later used to define appropriate :literal:`propagatorSettings` object for the numerical propagation. The function :literal:`computeSemiAnalyticalAndFullPropagation` which performs the propagation of the fully perturbed problem takes a :literal:`propagatorSettings` object as input, but the set of accelerations used to define them must contain **perturbing** accelerations only. So gravitational acceleration exerted by the central body and thrust acceleration should not be considered here, as they are already taken into account by the shaping method.


Set up hodographic shaping
~~~~~~~~~~~~~~~~~~~~~~~~~~

Using hodographic shaping to design a low-thrust trajectory requires the definition of three different shaping functions, one for each of the cylindrical velocity components. The shaping functions are defined as a combination of so-called base functions. As described in the hodographic shaping documentation (:ref:`tudatFeaturesHodographicShaping`), at least three base functions must be defined for each velocity component to ensure the boundary conditions are fulfilled.

The frequency of any trigonometric-like base functions and the scale factor used are defined as follows (recommended values depend on the time-of-flight):

.. code-block:: cpp

	double frequency = 2.0 * mathematical_constants::PI / timeOfFlight;
   	double scaleFactor = 1.0 / timeOfFlight;

The settings for each of the radial velocity component base functions are then defined:

.. code-block:: cpp

    // Create base function settings for the components of the radial velocity composite function.
    std::shared_ptr< BaseFunctionHodographicShapingSettings > firstRadialVelocityBaseFunctionSettings =
            std::make_shared< BaseFunctionHodographicShapingSettings >( );
    std::shared_ptr< BaseFunctionHodographicShapingSettings > secondRadialVelocityBaseFunctionSettings =
            std::make_shared< PowerFunctionHodographicShapingSettings >( 1.0, scaleFactor );
    std::shared_ptr< BaseFunctionHodographicShapingSettings > thirdRadialVelocityBaseFunctionSettings =
            std::make_shared< PowerFunctionHodographicShapingSettings >( 2.0, scaleFactor );


The function :literal:`createBaseFunctionHodographicShaping` can be called to create the corresponding base functions from the above-defined settings. The base functions defined that way are pushed to a :literal:`std::vector< std::shared_ptr< BaseFunctionHodographicShaping > >` object, which is latter used as an input parameter to create the :literal:`HodographicShaping` object.

.. code-block:: cpp
	
    // Create components of the radial velocity composite function.
    std::vector< std::shared_ptr< BaseFunctionHodographicShaping > > radialVelocityFunctionComponents;
    radialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( constant, firstRadialVelocityBaseFunctionSettings ) );
    radialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPower, secondRadialVelocityBaseFunctionSettings ) );
    radialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPower, thirdRadialVelocityBaseFunctionSettings ) );


A similar process is repeated for the normal and axial components of the spacecraft cylindrical velocity.

.. code-block:: cpp

    // Create base function settings for the components of the normal velocity composite function.
    std::shared_ptr< BaseFunctionHodographicShapingSettings > firstNormalVelocityBaseFunctionSettings =
            std::make_shared< BaseFunctionHodographicShapingSettings >( );
    std::shared_ptr< BaseFunctionHodographicShapingSettings > secondNormalVelocityBaseFunctionSettings =
            std::make_shared< PowerFunctionHodographicShapingSettings >( 1.0, scaleFactor );
    std::shared_ptr< BaseFunctionHodographicShapingSettings > thirdNormalVelocityBaseFunctionSettings =
            std::make_shared< PowerFunctionHodographicShapingSettings >( 2.0, scaleFactor );

    // Create components of the normal velocity composite function.
    std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > > normalVelocityFunctionComponents;
    normalVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( constant, firstNormalVelocityBaseFunctionSettings ) );
    normalVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPower, secondNormalVelocityBaseFunctionSettings ) );
    normalVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPower, thirdNormalVelocityBaseFunctionSettings ) );

    // Create base function settings for the components of the axial velocity composite function.
    std::shared_ptr< BaseFunctionHodographicShapingSettings > firstAxialVelocityBaseFunctionSettings =
            std::make_shared< TrigonometricFunctionHodographicShapingSettings >( ( numberOfRevolutions + 0.5 ) * frequency );
    std::shared_ptr< BaseFunctionHodographicShapingSettings > secondAxialVelocityBaseFunctionSettings =
            std::make_shared< PowerTimesTrigonometricFunctionHodographicShapingSettings >
            ( 3.0, ( numberOfRevolutions + 0.5 ) * frequency, scaleFactor );
    std::shared_ptr< BaseFunctionHodographicShapingSettings > thirdAxialVelocityBaseFunctionSettings =
            std::make_shared< PowerTimesTrigonometricFunctionHodographicShapingSettings >(
                3.0, ( numberOfRevolutions + 0.5 ) * frequency, scaleFactor );

    // Create components for the axial velocity composite function.
    std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > > axialVelocityFunctionComponents;
    axialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( cosine, firstAxialVelocityBaseFunctionSettings ) );
    axialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPowerCosine, secondAxialVelocityBaseFunctionSettings ) );
    axialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPowerSine, thirdAxialVelocityBaseFunctionSettings ) );


In hodographic shaping, the values of the free coefficients, if any, should be provided as inputs of the :literal:`HodographicShaping` object constructor. Here, only three base functions are provided per velocity component, which corresponds to the minimum required to satisfy the boundary conditions. So the free coefficients vectors are here just empty vectors. 

.. code-block:: cpp

    // Initialize free coefficients vector for radial velocity function (empty here, only 3 base functions).
    Eigen::VectorXd freeCoefficientsRadialVelocityFunction = Eigen::VectorXd::Zero( 0 );

    // Initialize free coefficients vector for normal velocity function (empty here, only 3 base functions).
    Eigen::VectorXd freeCoefficientsNormalVelocityFunction = Eigen::VectorXd::Zero( 0 );

    // Initialize free coefficients vector for axial velocity function (empty here, only 3 base functions).
    Eigen::VectorXd freeCoefficientsAxialVelocityFunction = Eigen::VectorXd::Zero( 0 );

Finally, the :literal:`HodographicShaping` object can be created:

.. code-block:: cpp

    // Create hodographic-shaping object with defined velocity functions and boundary conditions.
    shape_based_methods::HodographicShaping hodographicShaping(
                cartesianStateAtDeparture, cartesianStateAtArrival, timeOfFlight, 1, bodyMap, "Borzi", "Sun",
                radialVelocityFunctionComponents, normalVelocityFunctionComponents, axialVelocityFunctionComponents,
                freeCoefficientsRadialVelocityFunction, freeCoefficientsNormalVelocityFunction, freeCoefficientsAxialVelocityFunction,
                integratorSettings );


Set up spherical shaping
~~~~~~~~~~~~~~~~~~~~~~~~

The definition of a spherically shaped trajectory is much more straightforward than that of a hodographically shaped one. This is mostly due to the fact that the base functions used to map the spherical position of the spacecraft in spherical shaping are fixed, while those used in hodographic shaping have to be selected by the user. Also, there is no free parameters in spherical shaping, so the shaping function is pre-defined and cannot be tuned by the user.   

.. code-block:: cpp

    // Define root finder settings (used to update the value of the free coefficient, so that it matches the required time of flight).
    std::shared_ptr< root_finders::RootFinderSettings > rootFinderSettings =
            std::make_shared< root_finders::RootFinderSettings >( root_finders::bisection_root_finder, 1.0e-6, 30 );

    // Compute shaped trajectory.
    shape_based_methods::SphericalShaping sphericalShaping = shape_based_methods::SphericalShaping(
                cartesianStateAtDeparture, cartesianStateAtArrival, timeOfFlight,
                numberOfRevolutions, bodyMap, "Borzi", "Sun", 0.000703,
                rootFinderSettings, 1.0e-6, 1.0e-1, integratorSettings );


Retrieve trajectory, mass, thrust, and thrust acceleration profiles
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

From any :literal:`ShapeBasedMethodLeg` object, the trajectory of the spacecraft, along with the corresponding mass, thrust, and thrust acceleration profiles can be retrieved. This is done here for both the hodographically and spherically shaped trajectories. The code used for hodographic shaping is reproduced below (and it is done for spherical shaping in exactly the same way).

.. code-block:: cpp

    // Hodographic shaping
    std::vector< double > epochsVectorHodographicShaping;
    for ( std::map< double, Eigen::Vector6d >::iterator itr = hodographicShapingAnalyticalResults.begin( ) ;
          itr != hodographicShapingAnalyticalResults.end( ) ; itr++ )
    {
        epochsVectorHodographicShaping.push_back( itr->first );
    }

    std::map< double, Eigen::VectorXd > hodographicShapingMassProfile;
    std::map< double, Eigen::VectorXd > hodographicShapingThrustProfile;
    std::map< double, Eigen::VectorXd > hodographicShapingThrustAccelerationProfile;

    hodographicShaping.getMassProfile(
                epochsVectorHodographicShaping, hodographicShapingMassProfile, specificImpulseFunction, integratorSettings );
    hodographicShaping.getThrustProfile(
                epochsVectorHodographicShaping, hodographicShapingThrustProfile, specificImpulseFunction, integratorSettings );
    hodographicShaping.getThrustAccelerationProfile(
                epochsVectorHodographicShaping, hodographicShapingThrustAccelerationProfile, specificImpulseFunction, integratorSettings );

The plot below presents the Earth-Mars trajectories obtained with both hodographic (red) and spherical (blue) shaping methods. The associated thrust acceleration, thrust, and mass profiles are plotted too.

.. figure:: images/shapeBasedProfiles.png

Numerically propagate the unperturbed problem
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The unperturbed problem (with central body gravitational acceleration and spacecraft low-thrust acceleration as defined by the shaping method) is propagated numerically. To this end, the method :literal:`computeSemiAnalyticalAndFullPropagation` of the :literal:`ShapeBasedMethodLeg` is called and provides the analytical, shape-based trajectory, and the result of the numerically propagated trajectory.

.. code-block:: cpp

    std::map< double, Eigen::VectorXd > hodographicShapingFullPropagationResults;
    std::map< double, Eigen::Vector6d > hodographicShapingAnalyticalResults;
    std::map< double, Eigen::VectorXd > hodographicShapingDependentVariablesHistory;

    // Create propagator settings for hodographic shaping.
    std::pair< std::shared_ptr< propagators::PropagatorSettings< double > >, std::shared_ptr< propagators::PropagatorSettings< double > > >
            hodographicShapingPropagatorSettings = hodographicShaping.createLowThrustPropagatorSettings(
                specificImpulseFunction, perturbingAccelerationsMapSimplifiedProblem, integratorSettings, dependentVariablesToSave );

    // Compute shaped trajectory and propagated trajectory.
    hodographicShaping.computeSemiAnalyticalAndFullPropagation(
                integratorSettings, hodographicShapingPropagatorSettings, hodographicShapingFullPropagationResults,
                hodographicShapingAnalyticalResults, hodographicShapingDependentVariablesHistory );


    std::map< double, Eigen::VectorXd > sphericalShapingFullPropagationResults;
    std::map< double, Eigen::Vector6d > sphericalShapingAnalyticalResults;
    std::map< double, Eigen::VectorXd > sphericalShapingDependentVariablesHistory;

    // Create propagator settings for spherical shaping.
    std::pair< std::shared_ptr< propagators::PropagatorSettings< double > >, std::shared_ptr< propagators::PropagatorSettings< double > > >
            sphericalShapingPropagatorSettings = sphericalShaping.createLowThrustPropagatorSettings(
                specificImpulseFunction, perturbingAccelerationsMapSimplifiedProblem, integratorSettings, dependentVariablesToSave );

    // Compute shaped trajectory and propagated trajectory.
    sphericalShaping.computeSemiAnalyticalAndFullPropagation(
                integratorSettings, sphericalShapingPropagatorSettings, sphericalShapingFullPropagationResults,
                sphericalShapingAnalyticalResults, sphericalShapingDependentVariablesHistory );
	

Numerically propagate the perturbed problem
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The same is done using a different, more complete set of perturbing accelerations. The only difference is that the set of :literal:`PropagatorSettings` is defined differently, using the more complete set of perturbing accelerations that have been defined previously (gravitational attractions from Earth, Mars and Jupiter, and solar radiation pressure). Otherwise, the code is strickly the same as the one used to propagate the unperturbed problem. 

Results
~~~~~~~

The application output should look like:

.. code-block:: cpp

	Starting C:\tudatBundle\tudatExampleApplications\satellitePropagatorExamples\bin\applications\application_ShapeBasedTrajectoryDesign.exe...
	Dependent variables being saved, output vectors contain: 
	Vector entry, Vector contents
	Dependent variables being saved, output vectors contain: 
	Vector entry, Vector contents
	Dependent variables being saved, output vectors contain: 
	Vector entry, Vector contents
	Dependent variables being saved, output vectors contain: 
	Vector entry, Vector contents
	Dependent variables being saved, output vectors contain: 
	Vector entry, Vector contents
	Dependent variables being saved, output vectors contain: 
	Vector entry, Vector contents
	Dependent variables being saved, output vectors contain: 
	Vector entry, Vector contents
	Dependent variables being saved, output vectors contain: 
	Vector entry, Vector contents
	deltaV hodographic shaping: 21051.4

	deltaV spherical shaping: 5698.9

	C:/tudatBundle/tudatExampleApplications/satellitePropagatorExamples/bin/applications/application_ShapeBasedTrajectoryDesign.exe exited with code 0

The results of the numerical propagation (for both the unperturbed and perturbed cases) obtained with the two different shape-based methods are presented in the plot below. The difference in position between the analytical solution (so shaped trajectory, by definition computed under simplifying assumptions) and the full propagation numerical solution is plotted. The fully perturbed trajectory is propagated from half of the time-of-flight, backwards until departure, and forwards until arrival. This explains why the difference between analytical and numerical solutions is always zero in the middle of the trajectory, and grows larger when getting closer to either departure or arrival, as the effects of the perturbing accelerations keep propagating and adding up to each other. 

.. figure:: images/analyticalVsPropagationShapingMethods.png

In the unperturbed case, the analytical and numerical solutions are extremely similar. This is in agreement with the fact that no additional perturbing accelerations are considered in the numerical propagation of the problem, so that the observed differences are only due to integration errors. The fact that the difference between analytical and numerical results is (significantly) higher for spherical shaping can be explained by its independent variable which is not time, but azimuth angle. This requires an additional step to convert time to azimuth angle (and the other way around), compared to hodographic shaping. This conversion makes use of an interpolator, and the higher differences are due to interpolation errors (which unfortunately build up along the propagation).

The perturbed case shows larger differences between the shape-based and the propagated trajectories, because of the perturbing accelerations which are introduced in the numerical propagation.



