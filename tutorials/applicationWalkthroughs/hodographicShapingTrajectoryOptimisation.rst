.. _walkthroughsHodographicShapingOptimisation:

Shaping Methods Trajectory Optimisation
=======================================
The example described in this tutorial presents the use of shape-based methods to find the best trajectory to transfer from Earth to Mars. The code for this tutorial is given on Github, and is also located in your Tudat bundle at:

   tudatBundle/tudatExampleApplications/libraryExamples/PaGMOEx/lowThrustTrajectoryExample.cpp

This tutorial presents the optimisation of a shape-based trajectory, using hodographic shaping. As described in more details here (ADD LINK), this shaping method allows the user to introduce some extra degrees of freedom in the problem, by adding more base functions than what is required to satisfy the boundary conditions. The weighting coefficients associated with those additional base functions are then free parameters of the problem and can be tuned to minimize the deltaV required by the shaped trajectory. 

In this tutorial, we determine the best departure date and time-of-flight for an Earth-Mars transferm using hodographic shaping to design the trajectory. We first use the lowest-order solution (no additional base functions, so no degree of freedom in the design problem). Then, we optimise the shaped trajectory by introducing some free parameters (high-order solution) while focusing on a reduced search space. 

Set up the trajectory design problem  FOR THE GLORY OF THE BORZI!!!
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The environment within which the trajectory is to be computed is defined similarly to what is done in the previous tutorials. The central body (Sun) is created using default settings and its position is fixed to the origin of the inertial reference frame. A vehicle body is also created, and its initial mass is specified.

Grid search for the lowest-order solution
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Two nested for loops are used to parse the departure date and time-of-flight ranges for hodographically shaped Earth-Mars transfers. The departure dates range from 7304 MJD (Modified Julian Date) to 10225 MJD, while the time-of-flight are constrained between 500 and 2000 days.

.. code-block:: cpp

    std::pair< double, double > departureTimeBounds = std::make_pair( 7304.5 * physical_constants::JULIAN_DAY, 10225.5 * physical_constants::JULIAN_DAY  );
    std::pair< double, double > timeOfFlightBounds = std::make_pair( 500.0 * physical_constants::JULIAN_DAY, 2000.0 * physical_constants::JULIAN_DAY ); 

The time-of-flight step of the grid search is 5 days, while the departure dates are parsed every TO BE COMPLETED days. For each combination of departure date and time-of-flight, different numbers of revolutions are also parsed (N=0-5), and only the one leading to the best trajectory (lowest deltaV) is saved. Because the recommended base functions depend on the value of the time-of-flight and on the number of revolutions, the base functions settings need to be define inside the for loops. More precisely, the base functions of the radial and normal velocity components depend on the time-of-flight only, while the base functions for the axial velocity component also depend on the numbe of revolutions. The global code structure for the grid search is the following one:

.. code-block:: cpp

    int numberCases = 0;

    for ( int i = 0 ; i <= ( timeOfFlightBounds.second - timeOfFlightBounds.first ) / ( 5.0 * physical_constants::JULIAN_DAY ) ; i++  )
    {
        double currentTOF = timeOfFlightBounds.first + i * 5.0 * physical_constants::JULIAN_DAY;


        double frequency = 2.0 * mathematical_constants::PI / currentTOF;
        double scaleFactor = 1.0 / currentTOF;

        // Create base function settings for the components of the radial velocity composite function.
        ...

        // Create components of the radial velocity composite function.
        ...

        // Create base function settings for the components of the normal velocity composite function.
        ...

        // Create components of the normal velocity composite function.
        ...


        for ( int j = 0 ; j <= 400; j++ )
        {
            double currentDepartureDate = departureTimeBounds.first + j * ( departureTimeBounds.second - departureTimeBounds.first ) / 400.0;

            cartesianStateAtDeparture = pointerToDepartureBodyEphemeris->getCartesianState( currentDepartureDate );
            cartesianStateAtArrival = pointerToArrivalBodyEphemeris->getCartesianState( currentDepartureDate + currentTOF );


            int bestNumberOfRevolutions = 0;
            double currentBestDeltaV = hodographicShaping.computeDeltaV( );
            for ( int k = 0 ; k <= 5 ; k++ )
            {
                // Create base function settings for the components of the axial velocity composite function.
                ...

                // Set components for the axial velocity function.
                ...

                tudat::shape_based_methods::HodographicShaping hodographicShaping = shape_based_methods::HodographicShaping(
                            cartesianStateAtDeparture, cartesianStateAtArrival, currentTOF, k, bodyMap, "Vehicle", "Sun",
                            radialVelocityFunctionComponents, normalVelocityFunctionComponents, axialVelocityFunctionComponents,
                            freeCoefficientsRadialVelocityFunction, freeCoefficientsNormalVelocityFunction, freeCoefficientsAxialVelocityFunction );

                if ( hodographicShaping.computeDeltaV( ) < currentBestDeltaV )
                {
                    currentBestDeltaV = hodographicShaping.computeDeltaV( );
                    bestNumberOfRevolutions = k;
                }
	    }

	    // Save the results
	    ...

    	    numberCases++;
    	    hodographicShapingResults[ numberCases ] = TOFdepartureTimeDeltaV;

        }
    }



This grid search based on the lowest-order hodographic shaping solution provides the following Porkchop plot (in agreement with the results presented in ADD REFERENCE):

.. figure:: images/porkchopHodographicShapingLowOrder.png


Optimisation of the shaped trajectories (high-order solution)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Here, we introduce two degrees of freedom to the trajectory design problem, by adding to extra base functions to the composite function mapping the radial velocity component. We then want to optimise their values to minimise the deltaV required by the trajectory. Because of the computational load of the optimisation process, a reduced search space is considered here. The departure dates range from TO BE COMPLETED to TO BE COMPLETED, while the time-of-flight search space is reduced to the TO BE COMPLETED interval.

.. note:: 

	Only the definition of the base functions used to shape the radial velocity component is modified compared to the lowest-order grid search presented above.

Still, the global structure is similar to the one presented for the lowest-order solution presented in the first part of this tutorial. One major difference is that we no longer iterate on the number of revolutions to find the one leading the best trajectory, but the number of revolutions is now set to 1 (again, to reduce the computational load). The code is written as follows

.. code-block:: cpp

    for ( int i = 0 ; i <= 20 ; i++  )
    {
        double currentTOF = timeOfFlightBounds.first + i * 20.0 * physical_constants::JULIAN_DAY; 

        double frequency = 2.0 * mathematical_constants::PI / currentTOF;
        double scaleFactor = 1.0 / currentTOF;

        // Create base function settings for the components of the radial velocity composite function.
        ...

        // Create components of the radial velocity composite function.
        ...

        // Create base function settings for the components of the normal velocity composite function.
        ...

        // Create components of the normal velocity composite function.
        ...


        for ( int j = 0 ; j <= 5; j++ )
        {
            double currentDepartureDate = departureTimeBounds.first +
                    j * ( departureTimeBounds.second - departureTimeBounds.first ) / 200.0; 

            cartesianStateAtDeparture = pointerToDepartureBodyEphemeris->getCartesianState( currentDepartureDate );
            cartesianStateAtArrival = pointerToArrivalBodyEphemeris->getCartesianState( currentDepartureDate + currentTOF );


            // Create base function settings for the components of the axial velocity composite function.
	    ...

            // Set components for the axial velocity function.
            ... 

            std::vector< std::vector< double > > bounds( 2, std::vector< double >( 2, 0.0 ) );

            // Define search bounds: first parameter is start date, following parameters are leg durations
            bounds[ 0 ][ 0 ] = - 600.0;
            bounds[ 1 ][ 0 ] = 800.0;
            bounds[ 0 ][ 1 ] = 0.0;
            bounds[ 1 ][ 1 ] = 1500.0;

	    // Define hodographic shaping optimisation problem.
            ...

            // Perform optimisation
            ...
          
            // Save the results.
            ...
            

            // Set the free coefficients to zero (equivalent to lowest-order solution).
            Eigen::VectorXd freeCoefficientsRadialVelocityFunction = Eigen::VectorXd::Zero( 2 );
            Eigen::VectorXd freeCoefficientsNormalVelocityFunction = Eigen::VectorXd::Zero( 0 );
            Eigen::VectorXd freeCoefficientsAxialVelocityFunction = Eigen::VectorXd::Zero( 0 );

            // Compute the lowest-order solution for 1 revolution and save the results for 
            // comparison purposes.
            ...

        }
    }


Focusing first on the definition of the base functions for the radial velocity components, five of them are defined in this example. This adds two degrees of freedom, since three base functions are required to satisfy the boundary conditions in the radial direction. The following piece of code is used to define those five radial base functions:

.. code-block:: cpp

        // Create base function settings for the components of the radial velocity composite function.
        std::shared_ptr< shape_based_methods::BaseFunctionHodographicShapingSettings > firstRadialVelocityBaseFunctionSettings =
                std::make_shared< shape_based_methods::BaseFunctionHodographicShapingSettings >( );
        std::shared_ptr< shape_based_methods::BaseFunctionHodographicShapingSettings > secondRadialVelocityBaseFunctionSettings =
                std::make_shared< shape_based_methods::PowerFunctionHodographicShapingSettings >( 1.0, scaleFactor );
        std::shared_ptr< shape_based_methods::BaseFunctionHodographicShapingSettings > thirdRadialVelocityBaseFunctionSettings =
                std::make_shared< shape_based_methods::PowerFunctionHodographicShapingSettings >( 2.0, scaleFactor );
        std::shared_ptr< shape_based_methods::BaseFunctionHodographicShapingSettings > fourthRadialVelocityBaseFunctionSettings =
                std::make_shared< shape_based_methods::PowerTimesTrigonometricFunctionHodographicShapingSettings >( 1.0, 0.5 * frequency, scaleFactor );
        std::shared_ptr< shape_based_methods::BaseFunctionHodographicShapingSettings > fifthRadialVelocityBaseFunctionSettings =
                std::make_shared< shape_based_methods::PowerTimesTrigonometricFunctionHodographicShapingSettings >( 1.0, 0.5 * frequency, scaleFactor );

        // Create components of the radial velocity composite function.
        std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > > radialVelocityFunctionComponents;
        radialVelocityFunctionComponents.push_back(
                    createBaseFunctionHodographicShaping( shape_based_methods::constant, firstRadialVelocityBaseFunctionSettings ) );
        radialVelocityFunctionComponents.push_back(
                    createBaseFunctionHodographicShaping( shape_based_methods::scaledPower, secondRadialVelocityBaseFunctionSettings ) );
        radialVelocityFunctionComponents.push_back(
                    createBaseFunctionHodographicShaping( shape_based_methods::scaledPower, thirdRadialVelocityBaseFunctionSettings ) );
        radialVelocityFunctionComponents.push_back(
                    createBaseFunctionHodographicShaping( shape_based_methods::scaledPowerSine, fourthRadialVelocityBaseFunctionSettings ) );
        radialVelocityFunctionComponents.push_back(
                    createBaseFunctionHodographicShaping( shape_based_methods::scaledPowerCosine, fifthRadialVelocityBaseFunctionSettings ) );


The hodographic shaping optimisation problem has been implemented in the class :literal:`HodographicShapingOptimisationProblem` (see ADD LINK for more details). Creating an object of this class automatically creates a PAGMO compatible optimisation problem whose design parameters are the free coefficients of the hodographic shaping method, and which aims at minimising the deltaV of the trajectory.

.. code-block:: cpp

	problem prob{ HodographicShapingOptimisationProblem( cartesianStateAtDeparture, cartesianStateAtArrival, currentTOF, 1, bodyMap, "Vehicle",
                                                                 "Sun", radialVelocityFunctionComponents, normalVelocityFunctionComponents,
                                                                 axialVelocityFunctionComponents, bounds ) };


Once the optimisation problem has been defined, the selection of the algorithm, the creation of the Island and the solving of the optimisation problem itself are done in a very similar manner to what is presented in the previous optimisation tutorials:

.. code-block:: cpp

	algorithm algo{ pagmo::sga( ) };

            // Create an island with 1024 individuals
            island isl{ algo, prob, 1024 /*24*/ };

            // Evolve for 100 generations
            for( int i = 0 ; i < 10; i++ )
            {
                isl.evolve( );
                while( isl.status( ) != pagmo::evolve_status::idle &&
                       isl.status( ) != pagmo::evolve_status::idle_error )
                {
                    isl.wait( );
                }
                isl.wait_check( ); // Raises errors
            }


The results obtained after optimising the shaped trajectory over the reduced search space are the following ones:

.. figure:: images/porkchopHodographicShapingLowVsHighOrder.png

The use of a global evolutionary algorithm (genetic algorithm) here is not ideal to tackle this kind of optimisation problem and thus does not guarantee convergence. Local optimizers are known to perform better in that case (ADD REFERENCE), but the NLOPT library used in PAGMO for local optimisation encounters issues when run on Windows, so that global optimisation has been implemented in this tutorial for system compatibility. However, it is still sufficient to see that introducing some degrees of freedom in the trajectory design can reduce the deltaV budget and thus leads to better preliminary designs.


Results
~~~~~~~ 

The output of the application should look as follows:

.. code-block:: cpp

	




