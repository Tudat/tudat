.. _applicationCMakeLists:

Application CMakelists
======================

This page contains a tutorial for the ``CMakeLists.txt`` files in your applications folder. First the template ``CMakeLists.txt`` file is explained. Then a step by step guide is provided to link your own code to your application. 


TemplateApplication CMake
~~~~~~~~~~~~~~~~~~~~~~~~~

The ``CMakeLists.txt`` of the TemplateApplication located in ``/tudatExampleApplications/templateApplication/TemplateApplication`` is used in this tutorial. 

.. code-block:: cmake

   # Specify minimum CMake version required.
   cmake_minimum_required(VERSION 2.6)

   # Specify project name.
   project(TemplateApplication)

These lines specify the minimum verion of CMake required for this project. In case the version of your CMake installation is too old it will give a warning. The next line is used to specify the name of your project. This also sets the ``PROJECT_SOURCE_DIR`` and ``PROJECT_BINARY_DIR`` variables used by CMake. 


.. code-block:: cmake

   # Load UserSettings.txt
   if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR})
     set(BUILD_STYLE "standalone")
     include("${CMAKE_CURRENT_SOURCE_DIR}/UserSettings.txt" OPTIONAL)
   else()
     set(BUILD_STYLE "part of ${CMAKE_PROJECT_NAME}")
     include("${CMAKE_CURRENT_SOURCE_DIR}/UserSettings.txt" OPTIONAL)
     include("${CMAKE_SOURCE_DIR}/UserSettings.txt" OPTIONAL)
     STRING(REGEX REPLACE ${CMAKE_SOURCE_DIR} "" RELATIVE_PROJECT_PATH ${CMAKE_CURRENT_SOURCE_DIR})
     set(RELATIVE_PROJECT_PATH "${RELATIVE_PROJECT_PATH}" CACHE STRING "Relative path wrt to project for function")
     # message(STATUS "Relative path (wrt to project): ${RELATIVE_PROJECT_PATH}")
   endif()

This block is used to load additional CMake settings provided in ``UserSettings.txt``,  generally there is no need to specify anything in this file. These lines should remain as is.

.. code-block:: cmake

   # Set CMake build-type. If it not supplied by the user (either directly as an argument of through
   # the "UserSettings.txt" file, the default built type is "Release".
   if((NOT CMAKE_BUILD_TYPE) OR (CMAKE_BUILD_TYPE STREQUAL "Release"))
     set(CMAKE_BUILD_TYPE Release)
   elseif(CMAKE_BUILD_TYPE STREQUAL "Debug")
     set(CMAKE_BUILD_TYPE Debug)
   endif()

   message(STATUS "<< ${PROJECT_NAME} (${CMAKE_BUILD_TYPE} - ${BUILD_STYLE}) >>")

This block of code sets the build type of your application. This will either create a "Release" build or a "Debug" build of your application. These lines should remain as is.

.. code-block:: cmake

   # Add local module path
   list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/CMakeModules")
   message(STATUS "CMake Module path(s): ${CMAKE_MODULE_PATH}")
   
   # Set compiler based on preferences (e.g. USE_CLANG) and system.
   include(compiler)

These lines add the modules of your CMake installation to the list of CMakeModules and reports their location. Furthermore, it obtains the information of the compiler you want to use as specified in the configuration of the project as described in :ref:`verifyKitsAndCMake`. These lines should remain as is.


.. code-block:: cmake

   # Define the directory with the source code.
   set(SRCROOT "${CMAKE_CURRENT_SOURCE_DIR}")

   # Define the code root directory.
   set(CODEROOT "${CMAKE_CURRENT_SOURCE_DIR}/..")

   # Set testing options based on platform.
   enable_testing()

   # Set lib and bin directories where static libraries and unit tests are built.
   if(NOT LIB_ROOT)
     set(LIB_ROOT "${CODEROOT}/lib")
   endif()
   if(NOT BIN_ROOT)
     set(BIN_ROOT "${CODEROOT}/bin")
   endif()

These lines set the variables for the local paths to be used later in the file.

.. code-block:: cmake


   # Set the global macros for setting up targets.
   macro(setup_executable_target target_name CUSTOM_OUTPUT_PATH)
     set_property(TARGET ${target_name} PROPERTY RUNTIME_OUTPUT_DIRECTORY "${BIN_ROOT}/applications")
     install(TARGETS ${target_name} RUNTIME DESTINATION "${BIN_ROOT}/applications")
   endmacro(setup_executable_target)

   macro(setup_library_target target_name CUSTOM_OUTPUT_PATH)
     set_property(TARGET ${target_name} PROPERTY LIBRARY_OUTPUT_DIRECTORY "${LIB_ROOT}")
     set_property(TARGET ${target_name} PROPERTY ARCHIVE_OUTPUT_DIRECTORY "${LIB_ROOT}")
   endmacro(setup_library_target)

   macro(setup_unit_test_target target_name CUSTOM_OUTPUT_PATH)
     set_property(TARGET ${target_name} PROPERTY RUNTIME_OUTPUT_DIRECTORY "${BIN_ROOT}/unit_tests")
     get_property(CUSTOM_TEST_PROGRAM_NAME TARGET ${target_name} PROPERTY OUTPUT_NAME)
     add_test("${target_name}" "${BIN_ROOT}/unit_tests/${target_name}")
   endmacro(setup_unit_test_target)

These lines create CMake macro's to be used later in the file.

.. code-block:: cmake

   # Find Eigen3 library on local system.
   find_package(Eigen3 REQUIRED)
   
   # Include Eigen3 directories.
   # Set CMake flag to suppress Eigen warnings (platform-dependent solution).
   if(NOT APPLE OR APPLE_INCLUDE_FORCE)
     include_directories(SYSTEM AFTER "${EIGEN3_INCLUDE_DIR}")
   else()
     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${EIGEN3_INCLUDE_DIR}\"")
   endif()
   
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
   find_package(Boost 1.45.0
                COMPONENTS thread date_time system unit_test_framework filesystem regex REQUIRED)

   # Include Boost directories.
   # Set CMake flag to suppress Boost warnings (platform-dependent solution).
   if(NOT APPLE OR APPLE_INCLUDE_FORCE)
     include_directories(SYSTEM AFTER "${Boost_INCLUDE_DIRS}")
   else()
     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${Boost_INCLUDE_DIRS}\"")
   endif()

   # Find Tudat library on local system.
   find_package(Tudat 2.0 REQUIRED)

   # Include Tudat directories.
   # Set CMake flag to suppress Tudat warnings (platform-dependent solution).
   if(NOT APPLE OR APPLE_INCLUDE_FORCE)
     include_directories(SYSTEM AFTER "${TUDAT_INCLUDE_DIR}")
   else()
     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${TUDAT_INCLUDE_DIR}\"")
   endif()

   # Find CSPICE library on local system.
   find_package(Spice)

   # Include CSpice directories.
   if(NOT APPLE OR APPLE_INCLUDE_FORCE)
     include_directories(SYSTEM AFTER "${SPICE_INCLUDE_DIR}")
   else( )
     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${SPICE_INCLUDE_DIR}\"")
   endif( )

   option(USE_NRLMSISE00 "build Tudat with NRLMSISE-00 enabled" ON)
   if(NOT USE_NRLMSISE00)
     message(STATUS "NRLMSISE-00 disabled!")
     add_definitions(-DUSE_NRLMSISE00=0)
   else()
     message(STATUS "NRLMSISE-00 enabled!")
     add_definitions(-DUSE_NRLMSISE00=1)
     # Find USE_NRLMSISE00 library on local system.
     find_package(NRLMSISE00)
   
     # Include NRLMSISE00 directories.
     if(NOT APPLE OR APPLE_INCLUDE_FORCE)
       include_directories(SYSTEM AFTER "${NRLMSISE00_INCLUDE_DIR}")
     else( )
       set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${NRLMSISE00_INCLUDE_DIR}\"")
     endif( )
   endif( )
   
   if(USE_CSPICE)
     list(APPEND TUDAT_APPLICATION_EXTERNAL_LIBRARIES
         tudat_spice_interface
         cspice)
   endif()

   if(USE_NRLMSISE00)
     list(APPEND TUDAT_APPLICATION_EXTERNAL_LIBRARIES
         nrlmsise00)
   endif()

These lines are used to find external libraries and include them in the list of directories used by CMake such that they can be found when called inside your code. It reports the to the user whether or not certain libraries are used. The lines:

.. code-block:: cmake

   # Find Boost libraries on local system.
   find_package(Boost 1.45.0
                COMPONENTS thread date_time system unit_test_framework filesystem regex REQUIRED)

are important. Boost has many components, many are not necessary for your application. Therefore, only relevant components are loaded by CMake to speed up the compilation and size of your application. In this case six components are selected. But for example when using the PaGMo library (used for optimization) some extra components (serialization, chrono, atomic) may be necessary to be included here. 

.. code-block:: cmake

   list(APPEND TUDAT_APPLICATION_PROPAGATION_LIBRARIES tudat_simulation_setup tudat_propagators
       tudat_aerodynamics tudat_geometric_shapes tudat_relativity tudat_gravitation tudat_mission_segments
       tudat_electro_magnetism tudat_propulsion tudat_ephemerides tudat_numerical_integrators tudat_reference_frames
       tudat_basic_astrodynamics tudat_input_output tudat_basic_mathematics tudat_propagators tudat_basics ${TUDAT_APPLICATION_EXTERNAL_LIBRARIES})

These lines create a list with all the libraries required for propagating in Tudat such that they can be loaded simply by inluding the ``TUDAT_APPLICATION_PROPAGATION_LIBRARIES`` variable in linking libraries as is done in the last lines of the ``CMakeLists.txt`` file:

.. code-block:: cmake

   # Add helloWorld application.
   add_executable(application_HelloWorld "${SRCROOT}/helloWorld.cpp")
   setup_executable_target(application_HelloWorld "${SRCROOT}")
   target_link_libraries(application_HelloWorld tudat_gravitation tudat_basic_astrodynamics ${Boost_LIBRARIES} )

The first line indicates the name of the "to be created" application and the main source file in this case ``helloWorld.cpp``. The second line determines where the executable is located after the application is succesfully build. The last line is used to indicate all the used libraries. 

Linking own code
~~~~~~~~~~~~~~~~

In case you write additional code to be used in your application structured in source and header files the following steps need to be taken:
   
**Step 1: set source files** 
   
.. code-block:: cmake

    set(SOURCES
      "${SRCROOT}/myCode.cpp"

**Step 2: set header files**

.. code-block:: cmake

   set(HEADERS
      "${SRCROOT}/myCode.h"


**Step 3: add static library**

.. code-block:: cmake

   add_library(libraryName STATIC ${SOURCES})
   setup_library_target(libraryName "${SRCROOT}")

**Step 4: add library to target_link list**

.. code-block:: cmake

   target_link_libraries(application_HelloWorld libraryName tudat_gravitation tudat_basic_astrodynamics ${Boost_LIBRARIES} )

.. warning::

   These libraries are read right to left. Make sure to include your created libraries to the left since it could require components of the tudat libraries.






 
  

