.. _walkthroughspropTargeting:

Propagation Targeting
=====================================
The example described in this tutorial will cover the following problem statement:

   A satellite is going to be launched on an elliptical orbit with min and max altitudes respectively 180 km and 40000 km with inclination i=0 (and argument of perigee omega = 0). A fixed target is set on the same plane at an altitude of 35000 km, fixed latitude of 30 deg. What is the value of the RAAN for which the satellite achieves a minimum approach distance from the target?

The code for this tutorial is given here on Github, and is also located in your tudat bundle at:

   tudatBundle/tudatExampleApplications/libraryExamples/PaGMOEx/propagationTargetingExample.cpp

This tutorial will use the same concepts described in previous examples, however for this example, the Tudat propagation is implemented in the problem. Thus this example will show how the PAGMO library and Tudat can work together. It is necessary that the previous examples are understood to be able to go through this example.

Set up the problem
~~~~~~~~~~~~~~~~~~~~~~
The problem is set-up as follows:

.. code-block:: cpp

    //Load spice kernels
    tudat::spice_interface::loadStandardSpiceKernels( );

    // Create object to compute the problem fitness; no perturbations
    double altitudeOfPerigee = 180000.0;
    double altitudeOfApogee = 40000000.0;
    double altitudeOfTarget = 35000000.0;
    double longitudeOfTarget = 30.0; // In degrees
    problem prob{PropagationTargetingProblem( altitudeOfPerigee, altitudeOfApogee, altitudeOfTarget,
                  longitudeOfTarget, false ) };

First, the Spice (ephemeris) kernels are loaded in to get an accurate ephemeris of all the bodies included in this example. Then the input to the problem is defined: the perigee and apogee of the targeter, and the altitude and longitude of the target. Then the problem is object is made in the same manner as was done in previous examples, using the :literal:`PropagationTargetingProblem(...)` as the UDP.

.. warning:: Don't load the spice kernels inside the UDP! This will cause the kernels to be loaded each time the UDP is called, which will slow the application down, and result in an error by spice when a certain limit is reached.

The UDP contains an empty constructor and a normal constructor. In the normal constructor, the final argument, :literal:`useExtendedDynamics`, is a boolean variable that determines if either only the Earth gravitational influence is used (:literal:`false`, default) or if also the lunar and solar gravitational influence are implememted. The :literal:`get_bounds()` is implemented in the same way as before, but the fitness function is implemented differently, which will be discussed later. To reduce the computation time of the optimization example, it is important that a part of the code is put in the constructor. 

First, in the constructor, the orbit of the targeter is defined:

.. code-block:: cpp

    // Definition of the orbit
    earthRadius_ = spice_interface::getAverageRadius( "Earth" );
    radiusOfPerigee_ =  earthRadius_ + altitudeOfPerigee_;
    radiusOfApogee_ = earthRadius_ + altitudeOfApogee_;
    earthGravitationalParameter_ = spice_interface::getBodyGravitationalParameter( "Earth" );
    semiMajorAxis_ = (radiusOfApogee_ + radiusOfPerigee_)/2.0;

The next part of the constructor is the set up of the integration, the bodies, and the environment. This is done in a similar manner as in the :ref:`walkthroughsUnperturbedEarthOrbitingSatellite`. Thus if a part of this tutorial is not clear, the reader is referred to this example (and the following examples on that page).

The following code shows the initialization of the integration, the bodies, and the environment:

.. code-block:: cpp

    //Integration time: half a orbit
    simulationStartEpoch_ = 0.0;
    simulationEndEpoch_ = 1.2 * mathematical_constants::PI *
            std::sqrt(pow(semiMajorAxis_,3)/earthGravitationalParameter_);

    // Create the body Earth from Spice interface
    std::map< std::string, std::shared_ptr< BodySettings > > bodySettings;
    if( useExtendedDynamics_ )
    {
        bodySettings =
                getDefaultBodySettings( {"Earth", "Moon", "Sun"}, simulationStartEpoch_ - 3600.0, simulationEndEpoch_ + 3600.0 );
        bodySettings[ "Moon" ]->rotationModelSettings->resetOriginalFrame( "J2000" );
        bodySettings[ "Moon" ]->ephemerisSettings->resetFrameOrientation( "J2000" );
        bodySettings[ "Sun" ]->rotationModelSettings->resetOriginalFrame( "J2000" );
        bodySettings[ "Sun" ]->ephemerisSettings->resetFrameOrientation( "J2000" );
    }
    else
    {
        bodySettings =
                getDefaultBodySettings( {"Earth"} );
    }
    bodySettings[ "Earth" ]->ephemerisSettings = std::make_shared< simulation_setup::ConstantEphemerisSettings >(
                Eigen::Vector6d::Zero( ), "SSB", "J2000" );
    bodySettings[ "Earth" ]->atmosphereSettings = NULL;
    bodySettings[ "Earth" ]->shapeModelSettings = NULL;

    bodySettings[ "Earth" ]->rotationModelSettings->resetOriginalFrame( "J2000" );
    bodySettings[ "Earth" ]->ephemerisSettings->resetFrameOrientation( "J2000" );


    //Create bodyMap and add the satellite as an empty body
    bodyMap_ = simulation_setup::createBodies( bodySettings );

It is important to realize why this is done in the constructor and not in the fitness function of the UDP. If this was put into the fitness function, the :literal:`bodyMap` would be created each time the fitness function would be called. The :literal:`simulation_setup::createBodies()` function will store the information every time it is called, and thus the RAM usage on your computer will increase over time, until there is no more RAM left and the program is terminated, or slowed down considerably.
 
After the constructor is setup, the next step is to define the target orbit using the input values of the UDP in the fitness function, and set-up the initial conditions for the satellite:

.. code-block:: cpp

    //Define position of the target at 35000 km from Earth at 30 deg latitude
    Eigen::Vector3d target;
    target[0] = (earthRadius_ + altitudeOfTarget_) * cos(longitudeOfTarget_*mathematical_constants::PI/180);
    target[1] = (earthRadius_ + altitudeOfTarget_) * sin(longitudeOfTarget_*mathematical_constants::PI/180);
    target[2] = 0.0;

    //Define initial position of satellite at the perigee
    Eigen::Vector6d initialKeplerElements;
    initialKeplerElements[ semiMajorAxisIndex ] = semiMajorAxis_;
    initialKeplerElements[ eccentricityIndex ] = (radiusOfApogee_ - radiusOfPerigee_)/(radiusOfApogee_ + radiusOfPerigee_);
    initialKeplerElements[ inclinationIndex ] = 35.0 * mathematical_constants::PI/180.0;
    initialKeplerElements[ argumentOfPeriapsisIndex ] = x[0] * mathematical_constants::PI/180.0;
    initialKeplerElements[ longitudeOfAscendingNodeIndex ] = x[1] * mathematical_constants::PI/180.0;
    initialKeplerElements[ trueAnomalyIndex ] = 0.0;

    const Eigen::Vector6d systemInitialState = convertKeplerianToCartesianElements(
    initialKeplerElements, earthGravitationalParameter_ );


The only acceleration models that are implemented are gravitational of nature. The propagator setting use the Cowell method and a RK4 integrator:

.. code-block:: cpp

    //Setup simulation. Simple Keplerian orbit (only central-gravity of Earth)
    std::vector< std::string > bodiesToPropagate = { "Satellite" };
    std::vector< std::string > centralBodies = { "Earth" };
    SelectedAccelerationMap accelerationMap;
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfSatellite;
    if( useExtendedDynamics_ )
    {
        accelerationsOfSatellite[ "Earth" ].push_back( std::make_shared< SphericalHarmonicAccelerationSettings >(
                                                           2, 2 ) );
        accelerationsOfSatellite[ "Moon" ].push_back( std::make_shared< AccelerationSettings >(
                                                          point_mass_gravity ) );
        accelerationsOfSatellite[ "Sun" ].push_back( std::make_shared< AccelerationSettings >(
                                                         point_mass_gravity ) );
        accelerationMap[ "Satellite" ] = accelerationsOfSatellite;
    }
    else
    {
        accelerationsOfSatellite[ "Earth" ].push_back( std::make_shared< AccelerationSettings >(
                                                           point_mass_gravity ) );
        accelerationMap[ "Satellite" ] = accelerationsOfSatellite;
    }
    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodyMap_, accelerationMap, bodiesToPropagate, centralBodies );

    //Setup propagator (cowell) and integrator (RK4 fixed stepsize)
    std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings< double > >
            ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState, simulationEndEpoch_ );
    std::shared_ptr< IntegratorSettings< > > integratorSettings =
            std::make_shared< IntegratorSettings< > >
            ( rungeKutta4, simulationStartEpoch_, fixedStepSize );

    //Start simulation
    SingleArcDynamicsSimulator< > dynamicsSimulator(
            bodyMap_, integratorSettings, propagatorSettings, true, false, false );


After the simulation is defined, it is time to actually define the optimization part of this example: 

.. code-block:: cpp

    //Retrieve results
    std::map< double, Eigen::VectorXd > integrationResult =
            dynamicsSimulator.getEquationsOfMotionNumericalSolution( );


    //Find minimum distance from target
    Eigen::Vector3d separationFromTarget =
            Eigen::Quaterniond( Eigen::AngleAxisd(
                                    -earthRotationRate * integrationResult.begin( )->first, Eigen::Vector3d::UnitZ( ) ) ) *
            ( integrationResult.begin( )->second.segment( 0, 3 ) )- target;

    double bestDistanceFromTarget = separationFromTarget.norm( );
    double timeForBestDistanceFromTarget = integrationResult.begin( )->first;

    for( std::map< double, Eigen::VectorXd >::const_iterator stateIterator = integrationResult.begin( );
         stateIterator != integrationResult.end( ); stateIterator++ )
    {
        separationFromTarget = Eigen::Quaterniond( Eigen::AngleAxisd(
                                                       -earthRotationRate * stateIterator->first, Eigen::Vector3d::UnitZ( ) ) ) *
                stateIterator->second.segment( 0, 3 ) - target;
        const double distanceFromTarget = separationFromTarget.norm( );

        if( distanceFromTarget < bestDistanceFromTarget )
        {
            bestDistanceFromTarget = distanceFromTarget;
            timeForBestDistanceFromTarget = stateIterator->first;
        }

    }

    std::vector< double > output = {bestDistanceFromTarget} ;

    return output;

First, the results of the simulation are stored in a variable: :literal:`integrationResult`. Then the smalles distance and time of the smallest distance between the targeter and the target are intialized. These variables are needed in the for-loop that follows to determine the best results of the simulation. This result is then stored and used as the fitness value. 

Selecting the Algorithm
~~~~~~~~~~~~~~~~~~~~~~~~
In this example, the de1220 algorithm is selected to optimize the trajectory:

.. code-block:: cpp

    // Instantiate a pagmo algorithm
    algorithm algo{de1220( )};

A grid-search is also performed, however, this is only done to compare the resuolts and it shall thus not be discussed here. 

Building the Island
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The island is built in the same way as in :ref:`walkthroughsHimmelblau`:

.. code-block:: cpp

        // Create an island with 128 individuals
        pagmo::population::size_type populationSize = 128;
        island isl{algo, prob, populationSize};


Perform the Optimization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Finally, the optimization is performed in the same manner as in :ref:`walkthroughsHimmelblau`:

.. code-block:: cpp

        // Evolve for 25 generations
        for( int i = 0; i < 25; i++ )
        {
            isl.evolve();
            while( isl.status()!=pagmo::evolve_status::idle )
                isl.wait();

            // Write current iteration results to file
            printPopulationToFile( isl.get_population( ).get_x( ), "targetingPropagation_" + std::to_string( i ) + "_" + std::to_string( i ) , false );
            printPopulationToFile( isl.get_population( ).get_f( ), "targetingPropagation_" + std::to_string( i ) + "_" + std::to_string( i ) , true );

            std::cout<<i<<std::endl;
        }


Perturbed Example
~~~~~~~~~~~~~~~~~
In the second half of the example, the code is repeated, but with extended dynamics. This optimization problem is almost the same as the previous problem, but there are some differences due to the fact that the convergence of this problem is non-trivial. 

The first difference is the fact that when the UDP is initialized, a boolean variable in the constructor is set to :literal:`true`, namely: :literal:`useExtendedDynamics`. This variable will make sure that the third-body lunar and solar gravitational perturbation is included in the dynamics. The second difference is that the initial population is set to the final population of the previous problem. This will help with the convergence of this problem, and is shown below:

.. code-block:: cpp

        // Create an empty population for perturbed problem
        population population_pert = population( prob_pert, 0 );

        // Retrieve population of unperturbed problem, and instantiate population of perturbed problem
        std::vector<vector_double> original_population = isl.get_population( ).get_x( );
        for( int k = 0; k < populationSize; k++ )
        {
            population_pert.push_back( original_population.at( k ) );
        }

The most important thing to remember here is the fact that the population needs to be set in a for-loop, not just in one line. The rest of the example is fairly straightforward and is left to the reader to understand.


Results
~~~~~~~
The output of the application should look as follows (specific numbers could be different):

.. code-block:: cpp

        Starting ...\tudatBundle.git\tudatExampleApplications\libraryExamples\bin\applications\application_PagmoPropagationTargetingExample.exe...

	Grid search 0
	Grid search 1
	Grid search 2
	...
	...
	...
	Grid search 97
	Grid search 98
	Grid search 99
	0
	1
	2
	...
	...
	...
	22
	23
	24
	Grid search 0
	Grid search 1
	Grid search 2
	...
	...
	...
	Grid search 97
	Grid search 98
	Grid search 99
	0
	1
	2
	3
	.../tudatBundle.git/tudatExampleApplications/libraryExamples/bin/applications/application_PagmoPropagationTargetingExample.exe exited with code 0



