.. _jsonInterface_tutorials_lifetimeMaximisationCPP:

.. role:: jsontype
.. role:: jsonkey

Perturbed Earth-orbiting Satellite Lifetime Maximisation Using a Custom C++ Application
=======================================================================================

In this tutorial, the creation of a custom JSON-based application is described. Read first the general introduction on :ref:`jsonInterface_customJsonApp`.

The example described on this page is identical to that described in :ref:`jsonInterface_tutorials_lifetimeMaximisation`. The only difference is that there, the :literal:`json_interface` application is used with many input files to run the different cases, while in this case only one input file, containing the shared settings, is used, and a custom C++ application is written to manually modify a few parameters for each case. The C++ code can be found in::

  tudatBundle/tudatExampleApplications/satellitePropagatorExamples/SatellitePropagatorExamples/lifetimeMaximisation.cpp

The file :class:`lifetimeMaximisation.json` used as input for all the propagations is identical to the file :class:`shared.json` created in :ref:`jsonInterface_tutorials_lifetimeMaximisation`.

In this case, we do not write a C++ application because we want to use Tudat features that cannot be provided through a JSON file, but because we want to avoid having to generate a different input file for each propagation. The disadvantage of this choice is that we cannot use parallel processing, as all the propagations will be run sequentially by a single process, while when using the :literal:`json_interface` with different input files, we are able to run many cases concurrently.

In this case there is no need to write a derived class of :class:`JsonSimulationManager` with custom implementations for virtual methods. We can simply modify the JSON object containing the settings read from the JSON file before it is actually used to set up the simulation objects. Thus, all the code we need to write is:

.. code-block:: cpp
  :linenos:

  #include <Tudat/JsonInterface/jsonInterface.h>
  #include <SatellitePropagatorExamples/applicationOutput.h>

  int main( )
  {
    const std::string cppFilePath( __FILE__ );
    const std::string cppFolder = cppFilePath.substr( 0, cppFilePath.find_last_of("/\\") + 1 );
    tudat::json_interface::JsonSimulationManager< > jsonSimulationManager( cppFolder + "lifetimeMaximisation.json" );

    const std::string outputDirectory = tudat_applications::getOutputPath( ) + "LifetimeMaximisation/";
    const unsigned int numberOfCases = 365;
    for ( unsigned int i = 0; i < numberOfCases; ++i )
    {
        // Notify on propagation start
        std::cout << "Running propagation " << i + 1 << " of " << numberOfCases << "..." << std::endl;

        // Define the initial and final epochs
        const double initialEpoch = i * tudat::physical_constants::JULIAN_DAY;
        jsonSimulationManager[ "initialEpoch" ] = initialEpoch;
        jsonSimulationManager[ "finalEpoch" ] = initialEpoch + tudat::physical_constants::JULIAN_YEAR;

        // Define the output file
        jsonSimulationManager[ "export" ][ 0 ][ "file" ] = outputDirectory + "day" + std::to_string( i + 1 ) + ".txt";

        // Create settings objects
        jsonSimulationManager.updateSettings( );

        // Propagate
        jsonSimulationManager.runPropagation( );

        // Export results
        jsonSimulationManager.exportResults( );

        // Silence unused key warnings after first propagation
        if ( i == 0 )
        {
            jsonSimulationManager[ "options" ][ "unusedKey" ] = tudat::json_interface::continueSilently;
        }
    }

    return EXIT_SUCCESS;
  }

In this example, we first create a :class:`JsonSimulationManager`, which will contain the JSON object obtained by parsing the specified input file. Then, we write a loop in which we modify this JSON object before it is actually used to set up the objects needed for the propagation (integrator, propagators, bodies, etc.) when we call the :literal:`updateSettings` method. We can access this JSON object by using the method :literal:`getJsonObject`. However, we can also modify this object by accessing the :literal:`at` method and the :literal:`[ ]` operators of the :class:`JsonSimulationManager` directly. Thus, these two lines are equivalent:

.. code-block:: cpp

  jsonSimulationManager.getJsonObject( )[ "initialEpoch" ] = 0;
  jsonSimulationManager[ "initialEpoch" ] = 0;

And so are these two lines too:

.. code-block:: cpp

  std::cout << jsonSimulationManager.getJsonObject( ).at( "initialEpoch" ) << std::endl;
  std::cout << jsonSimulationManager.at( "initialEpoch" ) << std::endl;

Inside the loop, in which we iterate for each of the propagations to be carried out, we modify the keys :jsonkey:`initialEpoch` and :jsonkey:`finalEpoch` of the JSON object. Additionally, we want each propagation to generate an output file with a unique name, so we also modify the key :jsonkey:`export[ 0 ].file`. Then, we can set up the simulation, run the propagation and export the results.

After the first propagation has been completed, we turn off warnings for unused keys. This is done to silence warnings about the key :jsonkey:`bodies.satellite.initialState` being unused. When running the first propagation, the Keplerian state defined in :jsonkey:`bodies.satellite.initialState` is converted to Cartesian and assigned to :jsonkey:`propagators[ 0 ].initialStates`. Further propagations find that the key :jsonkey:`propagators[ 0 ].initialStates` is defined, and thus they do not use the information at :jsonkey:`bodies.satellite.initialState`, resulting in an unused key warning if we do not set the key :jsonkey:`options.unusedKey` to :literal:`"continueSilently"`.
