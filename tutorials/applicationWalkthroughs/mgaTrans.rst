.. _walkthroughsmga:

Multiple Gravity Assist Transfer
=====================================
The example described in this tutorial will cover a transfer using multiple gravity assists. The code for this tutorial is given here on Github, and is also located in your tudat bundle at:

   tudatBundle/tudatExampleApplications/libraryExamples/PaGMOEx/mgaTransferExample.cpp

This tutorial will use the same concepts described in previous examples, but now applied to a complex problem that may require a high number of generations to find the optima. It is necessary that the previous examples are understood to be able to go through this example.

Set up the problem
~~~~~~~~~~~~~~~~~~~~~~
The gravity assist are done in the following order: Eart, Venus, Earth, Earth, Jupiter (EVEEJ). The decision variables in this example are the start date, and the duration of the legs. The bounds are defined as before:

.. code-block:: cpp

    // We have five decision variables each with a lower and upper
    // bound, create a vector of vectors that will contain these.
    int numberOfParameters = 5;
    std::vector< std::vector< double > > bounds( 2, std::vector< double >( numberOfParameters, 0.0 ) );

    // Define search bounds: first parameter is start date, following parameters are leg durations
    bounds[ 0 ][ 0 ] = 7304.5; //MJD2000
    bounds[ 1 ][ 0 ] = 7304.5 + 10 * 365; //MJD2000
    bounds[ 0 ][ 1 ] = 200;
    bounds[ 1 ][ 1 ] = 500;
    bounds[ 0 ][ 2 ] = 50;
    bounds[ 1 ][ 2 ] = 300;
    bounds[ 0 ][ 3 ] = 50;
    bounds[ 1 ][ 3 ] = 300;
    bounds[ 0 ][ 4 ] = 500;
    bounds[ 1 ][ 4 ] = 3000;

The decision variable is set up in the same manner as the last tutorial: the first index will determine if it is the lower or upper bound, and the second index determines the decision variable. For the input of the UDP, a vector containing the planets need to be passed to the constructor, this vector is set-up as follows:

.. code-block:: cpp

    // Define the problem: EVEEJ flyby sequence
    std::vector< int > flybySequence;
    flybySequence.push_back( 3 );
    flybySequence.push_back( 2 );
    flybySequence.push_back( 3 );
    flybySequence.push_back( 3 );
    flybySequence.push_back( 5 );

The index here represents the planet where a flyby will be performed.

This vector, plus the bounds, can then be passed to the UDP, as is done in the following part of the code:

.. code-block:: cpp

    // Create object to compute the problem fitness
    problem prob{ MultipleGravityAssist( bounds, flybySequence, true ) };


where :literal:`MultipleGravityAssist( bounds, flybySequence, true )` is the UDP. The true statement enables the UDP to also use the trip time as an objective for the optimalization. 

This UDP contains an empty constructor, and a constructor that initializes the problem. The first part of the constructor initializes the leg types of the trajectory:

.. code-block:: cpp

    // Specify required parameters
    // Specify the number of legs and type of legs.
    numberOfLegs_ = flybySequence.size( );
    legTypeVector_.resize( numberOfLegs_ );
    legTypeVector_[ 0 ] = mga_Departure;
    legTypeVector_[ numberOfLegs_ - 1 ] = capture;

    for(int i = 1; i < numberOfLegs_ - 1; i++){
        legTypeVector_[ i ] = mga_Swingby;
    }


The trajectory starts with a departure from Earth, and ends with a capture at Jupiter. All the legs between the departure and capture are swing-by legs. This application will not go into detail on the trajectory design setup, if this is not understood the reader is referred to the :ref:`interplanetaryTrajectoryDesign` application tutorial.

After the trajectory is setup, the parameters for the selected order of planets is intialized using a large switch statement:

.. code-block:: cpp

    
    // Create the ephemeris, gravitational parameter, and minimum pericentre vector.
    ephemerisVector_.resize( numberOfLegs_ );
    gravitationalParameterVector_.resize( numberOfLegs_ );
    minimumPericenterRadii_.resize( numberOfLegs_ );
    for(int i = 0; i < numberOfLegs_; i++)
    {
        switch(flybySequence[ i ])
        {
        case( 1 ):
            ephemerisVector_[ i ] = boost::make_shared< ephemerides::ApproximatePlanetPositions >
                    ( ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::mercury );
            gravitationalParameterVector_[ i ] = 2.2032E13;
            minimumPericenterRadii_[ i ] = 2639.7E3;
            break;
        case( 2 ):
            ephemerisVector_[ i ] = boost::make_shared< ephemerides::ApproximatePlanetPositions >
                    ( ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::venus );
            gravitationalParameterVector_[ i ] = 3.24859E14;
            minimumPericenterRadii_[ i ] = 6251.8E3;
            break;
	...
	...
	...
        case( 8 ):
            ephemerisVector_[ i ] = boost::make_shared< ephemerides::ApproximatePlanetPositions >
                    ( ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::neptune );
            gravitationalParameterVector_[ i ] = 6.836529E15;
            minimumPericenterRadii_[ i ] = 25000.0E3;
            break;
        case( 9 ):
            ephemerisVector_[ i ] = boost::make_shared< ephemerides::ApproximatePlanetPositions >
                    ( ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::pluto );
            gravitationalParameterVector_[ i ] = 8.71E11;
            minimumPericenterRadii_[ i ] = 1395.0E3;
            break;
        default:
            std::cerr<<"Planet in flyby sequence is not defined.";
        }
    }

    // Create departure and capture variables.
    semiMajorAxes_.resize( 2 );
    eccentricities_.resize( 2 );
    semiMajorAxes_ << std::numeric_limits< double >::infinity( ), 1.0895e8 / 0.02;
    eccentricities_ << 0., 0.98;

Now, the fitness function needs to be set up:

.. code-block:: cpp

	//! Implementation of the fitness function (return delta-v)
	std::vector<double> MultipleGravityAssist::fitness( const std::vector<double> &xv ) const
	{
	    // Sun gravitational parameter
	    const double sunGravitationalParameter = 1.32712428e20;

	    // Create variable vector.
	    Eigen::VectorXd variableVector ( numberOfLegs_ + 1 );

	    double TOF = 0;
	    for(int i = 0; i < numberOfLegs_ ; i++){
		variableVector[ i ] = xv[ i ];
		if( i > 0 ){
		    TOF += xv[i];
		}
	    }
	    variableVector[ numberOfLegs_ ] = 1;//dummy
	    variableVector *= physical_constants::JULIAN_DAY;

	    // Create the trajectory problem.
	    Trajectory mgaTraj( numberOfLegs_, legTypeVector_, ephemerisVector_,
		                  gravitationalParameterVector_, variableVector, sunGravitationalParameter,
		                  minimumPericenterRadii_, semiMajorAxes_, eccentricities_ );

	    // Start the deltaV vector.
	    double resultingDeltaV;
	    mgaTraj.calculateTrajectory( resultingDeltaV );

	    if (std::isnan(resultingDeltaV))
	    {
		resultingDeltaV = 1.0E10;
	    }

	    if ( useTripTime_ ){
		return { resultingDeltaV, TOF };
	    }
	    else {
		return { resultingDeltaV };
	    }

	}




The :literal:`fitness( const std::vector<double> &xv )` method is made in a similar way as the previous examples. The decision variables are the departure time and the time-of-flight, thus they need to be entered into the trajectory design code to be able to calculate the delta V. A check is also made, that if the final Delta V was not able to be calculated, it will give a large penalty. 

Selecting the Algorithm
~~~~~~~~~~~~~~~~~~~~~~~~
In this example, the nsga2 algorithm is selected to optimize the trajectory:

.. code-block:: cpp

    // Select NSGA2 algorithm for problem
    algorithm algo{nsga2( )};


Building the Island
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The island is built in the same way as in :ref:`walkthroughsHimmelblau`:

.. code-block:: cpp

        // Create an island with 1000 individuals
        island isl{algo, prob, 1000 };


Perform the Optimization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Finally, the optimization is performed in the same manner as in :ref:`walkthroughsHimmelblau`:

.. code-block:: cpp

        // Evolve for 512 generations
        for( int i = 0 ; i < 512; i++ )
        {
            isl.evolve();
            while( isl.status()!=pagmo::evolve_status::idle )
                isl.wait();

            // Write current iteration results to file
            printPopulationToFile( isl.get_population( ).get_x( ), "mo_mga_EVEEJ_" + std::to_string( i ), false );
            printPopulationToFile( isl.get_population( ).get_f( ), "mo_mga_EVEEJ_" + std::to_string( i ), true );
            std::cout<<i<<std::endl;
        }


Results
~~~~~~~
The application output should look like this: 

.. code-block:: cpp

	Starting ...\tudatBundle.git\tudatExampleApplications\libraryExamples\bin\applications\application_PagmoMgaTransferExample.exe...

        0
        1
        2
        ...
        ...
        ...
        509
        510
        511

        .../tudatBundle.git/tudatExampleApplications/libraryExamples/bin/applications/application_PagmoMgaTransferExample.exe exited with code 0



