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
    bounds[ 0 ][ 0 ] = 2458849.5;
    bounds[ 1 ][ 0 ] = 2458849.5 + 20 * 365;
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

This UDP contains an empty constructor, and a constructor that initializes the bounds, and the problem:

.. code-block:: cpp

    MultipleGravityAssist::MultipleGravityAssist( std::vector< std::vector< double > > &bounds,
                                              std::vector< int > flybySequence,
                                              const bool useTripTime ) :
    problemBounds_( bounds ), useTripTime_( useTripTime ) {

       mgaObject_.type = total_DV_orbit_insertion;
       mgaObject_.sequence = flybySequence;
       mgaObject_.rev_flag.resize( flybySequence.size( ) );

       for( unsigned int i = 0; i < flybySequence.size( ); i++ )
       {
           mgaObject_.rev_flag[ i ] = 0;
       }

       mgaObject_.Isp = 300.0;
       mgaObject_.mass = 8000.0;
       mgaObject_.DVlaunch = 0.0;

       mgaObject_.rp = 76.0E3;
       mgaObject_.e = 0.9;
   }

Where :literal:`mgaObject` is an object of type :literal:`mgaproblem`, which is defined in the AstroToolbox folder. 


.. code-block:: cpp

    //! Get bounds
    std::pair<std::vector<double>, std::vector<double> > MultipleGravityAssist::get_bounds() const {

        return { problemBounds_[0], problemBounds_[1] };
    }

    //! Implementation of the fitness function (return delta-v)
    std::vector<double> MultipleGravityAssist::fitness( const std::vector<double> &xv ) const{

        std::vector<double> rp;
        std::vector<double> DV;
        double obj_funct;
        int result = MGA( xv, mgaObject_, rp, DV, obj_funct );
        if( !useTripTime_ )
        {
            return { obj_funct };
        }
        else
        {
            double tof = 0.0;
            for( unsigned int i = 1; i < problemBounds_.at( 0 ).size( ); i++ )
            {
                tof += xv[ i ];
            }
            return { obj_funct, tof };

        }
    }

The code above implements the two necessary methods: :literal:`get_bounds()` and :literal:`fitness( const std::vector<double> &xv )`. The :literal:`get_bounds()` function is implemented in the same way as previous examples. The :literal:`fitness( const std::vector<double> &xv )` method is relatively simple, however it uses an external function called :literal:`MGA`, which calculates the delta-v for the proposed gravity assist trajectory. 

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



