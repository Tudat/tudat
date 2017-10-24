.. _tudatFeaturesPagmo2:

External Libraries: PaGMO 2
===========================


Download and build PaGMO 2
~~~~~~~~~~~~~~~~~~~~~~~~~~
Clone or download PaGMO 2 from::

    https://github.com/esa/pagmo2

Make sure to name the source folder :literal:`pagmo2/`.

.. warning:: You need Boost version 1.55 or higher in order to build PaGMO 2.

The PaGMO project contains its own CMakeLists.txt. In order to build PaGMO as a standalone release version simply run CMake in the folder.
In order to include PaGMO in your project you will need to specify whether you want to build PaGMO or PyGMO (Python interface for PaGMO).
You can either set either of the two options :literal:`PAGMO_BUILD_PAGMO` or :literal:`PAGMO_BUILD_PYGMO` according to your needs in the CMakeLists.txt file.

.. warning:: Only one of the two flags can be set to ON at anytime, so make sure to switch the other to OFF (by default PAGMO_BUILD_PAGMO is set to ON).

If you are building PaGMO inside your project, you can override the two options using the command :literal:`-option()`

    .. code-block:: cmake

        option(PAGMO_BUILD_PAGMO ON)
        option(PAGMO_BUILD_PYGMO OFF)

You can simply include PaGMO 2 by adding the source folder as a subdirectory.

    .. code-block:: cmake
        
        unset(pagmo_LIB_DEPENDS CACHE)
        add_subdirectory( "my_path/pagmo2/" )

The command :literal:`unset(pagmo_LIB_DEPENDS CACHE)` makes sure to eliminate residual dependencies in the cache in case of re-build of the package.

Create an optimization problem with PaGMO
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Let us define the Himmelblau's function as a PaGMO problem.
The Himmelblau's function is defined as: :math:`f(x)=(x^2+y-11)^2 + (x+y^2-7)^2` and has 

    * one local maximum :math:`f(− 0.270845, − 0.923039) = 181.617`;
    * four local minima:

        * :math:`f( 3.0 , 2.0 ) = 0.0`;
        * :math:`f ( − 2.805118 , 3.131312 ) = 0.0`;
        * :math:`f ( − 3.779310 , − 3.283186 ) = 0.0`;
        * :math:`f ( 3.584428 , − 1.848126 ) = 0.0`.

A PaGMO problem in C++ is defined as a ``struct`` with the following mandatory methods:

    * ``std::vector< double > fitness( const std::vector< double > &x ) const``
    * ``std::pair< std::vector< double >, std::vector< double > > get_bounds( ) const``

In a file ``himmelblaus.h`` let us combine the himmelblau's function and the problem and allow the user to define the **box-boundaries** for the optimization search.

    .. code-block:: cpp

        //File: himmelblaus.h

        struct my_problem {

            //Empty constructor
            my_problem( ){ }

            //Actual constructor allowing the user to define the boundaries
            my_problem( const double x_min, const double x_max, const double y_min, 
                    const double y_max ) :
                x_min_( x_min ), x_max_( x_max ), y_min_( y_min ), y_max( y_max_ )
            { } 
                
            // Mandatory, computes the fitness, i.e. the Himmelblau's function
            std::vector< double > fitness( const std::vector< double > &x ) const
            {
                std::vetor< double > return_value;
                     
                return_value.push_back( pow( x[0]*x[0] + x[1] - 11, 2) + pow( x[0] + x[1]*x[1] - 7, 2 ) );  
                
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

        private:

            //Storage members
            double x_min_;
            double x_max_;
            double y_min_;
            double y_max_;
        
        };

The **empty constructor** used in the example above will avoid multi-threading identity related problems associated with the ``pthread`` library.

The use of PaGMO enbedded alghorithms and the user-defined problem is demonstrated in the following code block. Let us define a file called ``main.cpp``.

    .. code-block:: cpp
        
        //File: main.cpp

        #include<iostream>
        #include"pagmo/algorithms/sade.hpp"
        #include"pagmo/island.hpp"
        #include"pagmo/problem.hpp"
        #include"himmelblaus.h"

        int main( )
        {

	    // Create a PaGMO problem from the user-defined Himmelblau's problem
	    // Search within 0<x<5 and 0<y<5
            pagmo::problem prob{ my_problem( 0, 5, 0, 5) };
            
	    // Use the self-adaptive differential evolutionary algorithm
            pagmo::algorithm algo{ pagmo::sade( ) };
                
            // Define a PaGMO island with a population of 8 individuals
            pagmo::island isl = pagmo::island{ algo, prob, 8 };
            
            // Evolve the algorithm 20 times
            for( int i = 1; i <= 20; i++ )
            {
             
                isl.evolve( );
                
                // Check for the status of the island before iterating
                while( isl.status()!=pagmo::evolve_status::idle )
                    isl.wait();
            
            }

            std::cout << "Best x: " << isl.get_population().champion_x()[0] << std::endl;
            std::cout << "Best y: " << isl.get_population().champion_x()[1] << std::endl;
            std::cout << "Minimum: " << isl.get_population().champion_f()[0] << std::endl;
        
        return 0;
        
        }


If no operation between evolutions need to be performed, the for-loop can be simply replaced with ``isl.evolve(20)``. Because of the parallel computing nature of the island, the 
alghorithm is evolved in a different thread, therefore a check on the island status is needed while using a for-loop.

Including PaGMO in your application
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
1. Copy :literal:`FindPaGMO2.cmake` and :literal:`FindBoost.cmake` into your :literal:`CMAKE_MODULE_PATH`.
2. Include the following lines in your CMakeLists.txt file to add Boost and PaGMO2 header only library:

    .. code-block:: cmake

        # Configure Boost libraries.
        if(NOT Boost_USE_STATIC_LIBS)
          set(Boost_USE_STATIC_LIBS ON)
        endif()
        if(NOT Boost_USE_MULTITHREADED)
          set(Boost_USE_MULTITHREADED ON)
        endif()
        if(NOT Boost_USE_STATIC_RUNTIME)
          set(Boost_USE_STATIC_RUNTIME ON)
        endif()

        # Find Boost libraries on local system.
        find_package(Boost 1.55.0 COMPONENTS thread date_time system unit_test_framework 
            serialization filesystem regex REQUIRED)

        # Include Boost directories.
        # Set CMake flag to suppress Boost warnings (platform-dependent solution).
        if(NOT APPLE)
          include_directories(SYSTEM AFTER "${Boost_INCLUDE_DIRS}")
        else()
          set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${Boost_INCLUDE_DIRS}\"")
        endif()

        # Find PaGMO library on local system.
        find_package(PaGMO2)

        # Include PaGMO directories.
        if(NOT APPLE)
          include_directories(SYSTEM AFTER "${PAGMO_INCLUDE_DIR}")
        else( )
          set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${PAGMO_INCLUDE_DIR}\"")
        endif( )


3. Include the following lines to set-up the bin application:

    .. code-block:: cmake

        # Add application.
        add_executable(application_PagmoHimmelblausExample "../my_path/main.cpp")


        # Set PaGMO debug and release flags
        target_compile_options(application_PagmoHimmelblausExample PRIVATE "$<$<CONFIG:DEBUG>:${PAGMO_CXX_FLAGS_DEBUG}>" "$<$<CONFIG:RELEASE>:${PAGMO_CXX_FLAGS_RELEASE}>")
        
        # Set C++ standard to C++11
        set_property(TARGET application_PagmoHimmelblausExample PROPERTY CXX_STANDARD 11)
        set_property(TARGET application_PagmoHimmelblausExample PROPERTY CXX_STANDARD_REQUIRED YES)
        set_property(TARGET application_PagmoHimmelblausExample PROPERTY CXX_EXTENSIONS NO)
        
        # Link Boost libraries and pthread
        target_link_libraries(application_PagmoHimmelblausExample ${Boost_LIBRARIES} pthread)



        

    

        

         
            

