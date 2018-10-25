.. _walkthroughsEarthMars:

(multi-objective) Earth Mars Transfer
=====================================
The examples described in this tutorial will cover a direct transfer from Earth to Mars using a Lambert targeter. The code for this tutorial is given on Github, and is also located in your Tudat bundle at:

   tudatBundle/tudatExampleApplications/libraryExamples/PaGMOEx/earthMarsTransferExample.cpp

and

   tudatBundle/tudatExampleApplications/libraryExamples/PaGMOEx/multiObjectiveEarthMarsTransferExample.cpp

This tutorial will  cover two examples, one with a single objective, and one with multiple objectives. It will also show how the bounds of the decision vector can be defined in a more intuitive manner, and how Pagmo can be integrated into Tudat. It is necessary that the previous example, the :ref:`walkthroughsHimmelblau`, is understood to be able to go through this example.

Set up the problem
~~~~~~~~~~~~~~~~~~~~~~
As mentioned earlier, in this example the definition of the bounds will be done in a more intuitive manner, as is shown below:

.. code-block:: cpp

    // We have two decision variables each with a lower and upper bound, create a vector of vectors that will contain these.
    std::vector< std::vector< double > > bounds( 2, std::vector< double >( 2, 0.0 ) );

    // Define bounds: Search between 2020 and 2025 for flight duration between 200 and 1000 days.
    bounds[ 0 ][ 0 ] = 2458849.5;
    bounds[ 1 ][ 0 ] = 2460676.5;
    bounds[ 0 ][ 1 ] = 200;
    bounds[ 1 ][ 1 ] = 1000;

In this example a two-dimensional vector is defined, wich is initialized to be of size 2-by-2, with 0.0 as initial values. The decision variable in this example contains two variables: the launch date and the flight duration. The first index of the vector represents the lower (0) and the upper (1) bound, the second index represents the launch date (0) and the flight duration (1). This vector can then be passed to the UDP, as is done in the following part of the code:

.. code-block:: cpp

    // Create object to compute the problem fitness
    problem prob{EarthMarsTransfer( bounds )};

and for the multi-objective example:

.. code-block:: cpp

    // Create object to compute the problem fitness
    problem prob{EarthMarsTransfer( bounds, true )};

where :literal:`EarthMarsTransfer( bounds )` is the UDP. The true statement enables the UDP to not only use the delta-V as an objective, but also the trip time. 

This UDP contains an empty constructor, and a constructor that initializes the bounds, just as in :ref:`walkthroughsHimmelblau`. Afterwards two methods are added:

.. code-block:: cpp

    //! Descriptive name of the problem
    std::string EarthMarsTransfer::get_name() const {
       return "Multi-revolution Lambert Earth-Mars transfer trajectory";
    }

    //! Get bounds
    std::pair<std::vector<double>, std::vector<double> > EarthMarsTransfer::get_bounds() const {
       return { problemBounds_[0], problemBounds_[1] };
    }

They are in charge of producing an descriptive output, and getting the bounds (which is a mandatory method). After this the fitness function is defined. This function integrates several tudat features into the problem, thus it is important that these features are first understood. This will not be done in this tutorial, readers who do not undderstand these features are referred to: :ref:`walkthroughsIndex` or :ref:`tudatFeaturesIndex` to get a better understanding of them. The fitness function looks as follows:

.. code-block:: cpp

    //! Implementation of the fitness function (return delta-v)
    std::vector<double> EarthMarsTransfer::fitness( const std::vector<double> &xv ) const{

       using tudat::mission_segments::MultiRevolutionLambertTargeterIzzo;

       std::vector<double> f;

       // Gravitational parameter of the Sun
       double mu = 1.32712440018e+20;

       // Set initial and final position as those of Earth and Mars at
       // departure and arrival respectively.

       StateType initialState = getPlanetPosition( xv[0], "Earth");

       StateType finalState   = getPlanetPosition( xv[0] + xv[1], "Mars" );

       MultiRevolutionLambertTargeterIzzo lambertTargeter( initialState.segment(0,3),
           finalState.segment(0,3), xv[1]*86400, mu );

       double deltaV = std::numeric_limits<double>::infinity();

       unsigned int maxrev = lambertTargeter.getMaximumNumberOfRevolutions( );

       // Go through all multi-revolution solutions and select the one
       // with the lowest delta-V
       for( unsigned int i = 0; i <= maxrev; ++i){
	   lambertTargeter.computeForRevolutionsAndBranch( i, false );
	   deltaV = std::min( deltaV, ( initialState.segment(3,3)
		   - lambertTargeter.getInertialVelocityAtDeparture( )).norm() +
	       + ( finalState.segment(3,3)
		   - lambertTargeter.getInertialVelocityAtArrival( )).norm());
       }

       f.push_back(deltaV);

       if( useTripTime_ )
       {
           f.push_back( xv[1] );
       }

       return f;
   }

After the necessary features and variables are defined, the first step in the fitness function is to get the initial state and the final state of the spacecraft. The initial state is the state of the Earth and the final state is that of Mars. These states are retrieved using a method in the same class called: :literal`getPlanetPosition()`, which doesn't contain any special Pagmo features, thus is not discussed in this tutorial. The following three lines initialize the lambert targeter using the intitial and final positions, and the decision variables. Finally, in the for-loop, the solutions of the Lambert targeter are evaluated and the one with the lowest delta-v is selected. For the multi-objective optimizer, the variable :literal:`useTripTime_` is set to true, thus the trip time for the selected trajectory is added to the fitness vector. 

Selecting the Algorithm
~~~~~~~~~~~~~~~~~~~~~~~~
In this example, 8 optimizers are compared. These optimizers are selected using a for-loop:

.. code-block:: cpp

    for( int j = 0; j < 8; j++ )
    {
        // Retrieve algorothm
        int algorithmIndex = j;
        algorithm algo{getAlgorithm( algorithmIndex )};
        ...
        ...
        ...

and for the multi-objective example

.. code-block:: cpp

    for( int j = 0; j < 8; j++ )
    {
        // Retrieve MO algorithm
        algorithm algo{getMultiObjectiveAlgorithm( j )};
        ...
        ...
        ...
              
where the :literal:`getAlgorithm()` function is used, just as in the previous example. For the multi-objective example, a similar function is used: :literal:`getMultiObjectiveAlgorithm()`, which works the same as :literal:`getAlgorithm()`, but only selects algorithms suitable for multi-objective problems. 

Building the Island
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The island is built in the same way as in :ref:`walkthroughsHimmelblau`:

.. code-block:: cpp

        // Create an island with 1024 individuals
        island isl{algo, prob, 1024};


Perform the Optimization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Finally, the optimization is performed in the same manner as in :ref:`walkthroughsHimmelblau`:

.. code-block:: cpp

                    ...
                    ...
                    ...
                    // Evolve for 100 generations
                    for( int i = 0 ; i < 100; i++ )
                    {
                         isl.evolve();
                         while( isl.status()!=pagmo::evolve_status::idle )
                         isl.wait();

                         // Write current iteration results to file
                         printPopulationToFile( isl.get_population( ).get_x( ), "earthMarsLambert_" + std::to_string( j ) + "_" + std::to_string( i ) , false );
                         printPopulationToFile( isl.get_population( ).get_f( ), "earthMarsLambert_" + std::to_string( j ) + "_" + std::to_string( i ) , true );

                          std::cout<<i<<" "<<algorithmIndex<<std::endl;
                     }
                 }


Results
~~~~~~~
The single-objective application output should look like this: 

.. code-block:: cpp

	Starting ...\tudatBundle\tudatExampleApplications\libraryExamples\bin\applications\application_PagmoEarthMarsTransferExample.exe...

	Grid search 0
	Grid search 1
	Grid search 2
	...
	...
	...
	Grid search 997
	Grid search 998
	Grid search 999
	0 0
	1 0
	2 0
	3 0
	...
	...
	...
	96 7
	97 7
	98 7
	99 7

	.../tudatBundle/tudatExampleApplications/libraryExamples/bin/applications/application_PagmoEarthMarsTransferExample.exe exited with code 0

For the multi-objective example, the following output plot can be seen here (there are several figures that all represent different generations, the 10th generation is shown here):

.. figure:: images/porkchop10.png


