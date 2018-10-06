.. _walkthroughsHimmelblau:

Himmelblau Optimization
====================================
The example described on this page is that of finding a minimum of the Himmelblau function. The Himmelblau function is often used to test optimization algorithms, thus it is a good start to understand the Pagmo 2 optimization library. The code for this tutorial is given on Github, and is also located in your tudat bundle at:

   tudatBundle/tudatExampleApplications/libraryExamples/PaGMOEx/himmelblauOptimization.cpp

The Himmelblau function is given by the following equation:

.. math:: f(x,y) = (x^2 + y - 11)^2 + (x - y^2 - 7)^2

and it has four local minima, all with the same functional value of 0.

The following concepts will be discussed in this tutorial:

    - Setting up a problem.
    - Choosing an optimization algorithm and its parameters.
    - Building an island for the population.

The Pagmo 2 library has several options that are not covered in this tutorial (and the upcoming tutorials). For more in depth information on the Pagmo 2 library, users are referred to the `documentation <https://esa.github.io/pagmo2/>`_ of this library. This tutorial, and the following tutorials, will also not go in depth of the theory behind the optimization algorithms. These tutorials will only contain information on how to implement optimization algorithms into your code.

Set up the problem
~~~~~~~~~~~~~~~~~~~~~~
In the Himmelblau optimization code, the set up of the problem is a very short piece of code:

.. code-block:: cpp

   //Set seed for reproducible results
   pagmo::random_device::set_seed( 12345 );

   // Define Himmelblau function (range [-5,5]x[-5,5])
   pagmo::problem prob{ HimmelblauFunction( -5, 5, -5, 5) };

The first line sets the seed for the problem. As most algorithms need random number generators, a seed can be set explicitly to be able to produce the same set of random numbers. The second line is the most important part of setting up a problem. The :literal:`pagmo::problem` class represents a generic evolutionary problem that tries to find a decision vector, :math:`\vec{x}`, bounded by an upper and lower boundary, to minimize the problem: :math:`f(\vec{x})`, subject to several (in-)equality constraints. In this case, a user defined problem (UDP) is set-up: the :literal:`HimmelblauFunction`. The header file of this UDP is found in the problems folder, located in the same directory as the :literal:`himmelblauOptimization.cpp` file. In this header file the following can be seen:

.. code-block:: cpp

    // Empty constructor
    // Without an empty constructor the problem is not accepted
    // as a multithreading type
    HimmelblauFunction( ){ }

    //Actual constructor allowing the user to define the boundaries
    HimmelblauFunction( const double x_min, const double x_max, const double y_min,
            const double y_max ) :
        x_min_( x_min ), x_max_( x_max ), y_min_( y_min ), y_max_( y_max )
    { }

The **empty constructor** used in the example above will avoid multi-threading identity related problems associated with the ``pthread`` library. The actual constructor below the empty constructor allows the user to define the bounds of the problem. For the Himmelblau function the decision vector has a size of two, thus there are four entries in the constructor (upper and lower bounds for each decision variable). A :literal:`pagmo::problem` class has two mandatory methods:

    - :literal:`std::vector< double > fitness( const std::vector< double > &x ) const`
    - :literal:`std::pair< std::vector< double >, std::vector< double > > get_bounds( ) const`

Which can be seen in the following code:

.. code-block:: cpp

    // Mandatory, computes the fitness, i.e. the Himmelblau's function
    std::vector< double > fitness( const std::vector< double > &x ) const
    {

        std::vector< double > return_value;

        return_value.push_back( pow( x[0]*x[0] + x[1] - 11.0, 2.0 )
                + pow( x[0] + x[1]*x[1] - 7.0, 2.0 ) );

        return return_value;

    }

    // Mandatory, returns the box-bounds
    std::pair< std::vector< double >, std::vector< double > > get_bounds( ) const
    {

        std::pair< std::vector< double >, std::vector< double > > box_bounds;

        box_bounds.first.push_back( x_min_ );
        box_bounds.first.push_back( y_min_ );

        box_bounds.second.push_back( x_max_ );
        box_bounds.second.push_back( y_max_ );

        return box_bounds;

     }

The first method defines the function that will be minimized (the fitness function), which in this case is the Himmelblau function. Its input is the decision vector and its output is the value of the fitness function. The second function is the other mandatory method that returns the bounds of the decision vector. The :literal:`std::pair` type has a :literal:`std::pair.first` and :literal:`std::pair.second` option which define the lower and upper bounds respectively. 

This header file is sufficient for the definition of this problem. In future tutorials, other options will be explored.

Selecting the Algorithm
~~~~~~~~~~~~~~~~~~~~~~~~
The :literal:`pagmo::algorithm` class is used to select an algorithm that solves the optimization problem. In the Himmelblau example it is used as follows:

.. code-block:: cpp

        // Solve using DE algorithm
        pagmo::algorithm algo{ pagmo::de( ) };

This creates an object, :literal:`algo`, that contains the :literal:`pagmo::de( )` algorithm. The algorithm can either be a pre-defined algorithm, that can be found in this `list <https://esa.github.io/pagmo2/docs/algorithm_list.html>`_, or it can be defined using an user defined algorithm (UDA). The second option will not be discussed in this tutorial. In this example, a differential evolutionary algorithm is chosen to solve the problem. 

Building the Island
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The next part concerns the :literal:`pagmo::island` class, and it looks as follows:

.. code-block:: cpp

        // Create island with 1000 individuals
        pagmo::island isl = pagmo::island{ algo, prob, 1000 };

This piece of code creates and object :literal:`isl` of type :literal:`pagmo::island`, which manages the evolution of the population, here defined to be 1000 random samples, using the :literal:`algo` object to solve the problem defined by the :literal:`prob` object. 

Perform the Optimization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The final part of the code is the actual optimization. This process is shown below:

.. code-block:: cpp

        // Evolve for 1000 generations
        for( int i = 1; i <= 100; i++ )
        {
            isl.evolve( );
            while( isl.status()!=pagmo::evolve_status::idle )
                isl.wait();

            printPopulationToFile( isl.get_population( ).get_x( ), "himmelblau_" + std::to_string( i ) , false );
            printPopulationToFile( isl.get_population( ).get_f( ), "himmelblau_" +  std::to_string( i ) , true );

            // Print current optimum to console
            std::cout << "Minimum: " <<i<<" "<<std::setprecision( 16 ) <<"f= "<< isl.get_population().champion_f()[0] <<", x="<<
                         isl.get_population().champion_x()[0] <<" y="<<isl.get_population().champion_x()[1] <<std::endl;
        
        }

The optimization is done using methods from the :literal:`pagmo::island` class. In a for-loop that runs for the desired amount of generations (here 1000), the population on the island is first evolved using the :literal:`evolve()` method. The evolve method will generate a new population on the island using the algorithm specified in :literal:`algo`. After this, the status of the island is checked before iterating. Finally, the decision vectors and the fitness values are then written to a file and the best member of the population (the champion) is also displayed. 


Results
~~~~~~~
The output of the program should look similar to the output below:

.. code-block:: cpp

        Starting ...\tudatBundle\tudatExampleApplications\libraryExamples\bin\applications\application_PagmoHimmelblauOptimization.exe...
        Minimum: 1 f= 0.0399354434859056, x=-2.786081426696553 y=3.10439194057476
        Minimum: 2 f= 0.0399354434859056, x=-2.786081426696553 y=3.10439194057476
        Minimum: 3 f= 0.0290540463514649, x=3.029539573873633 y=1.972445487248521
        Minimum: 4 f= 0.0290540463514649, x=3.029539573873633 y=1.972445487248521
        Minimum: 5 f= 0.0290540463514649, x=3.029539573873633 y=1.972445487248521
        ...
        ...
        ...
        Minimum: 95 f= 4.282411893194119e-010, x=3.0000021972896 y=1.999994663573049
        Minimum: 96 f= 4.282411893194119e-010, x=3.0000021972896 y=1.999994663573049
        Minimum: 97 f= 4.282411893194119e-010, x=3.0000021972896 y=1.999994663573049
        Minimum: 98 f= 4.282411893194119e-010, x=3.0000021972896 y=1.999994663573049
        Minimum: 99 f= 4.282411893194119e-010, x=3.0000021972896 y=1.999994663573049
        Minimum: 100 f= 4.282411893194119e-010, x=3.0000021972896 y=1.999994663573049
        .../tudatBundle/tudatExampleApplications/libraryExamples/bin/applications/application_PagmoHimmelblauOptimization.exe exited with code 0





