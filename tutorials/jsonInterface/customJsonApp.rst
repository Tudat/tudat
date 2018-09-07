.. _jsonInterface_customJsonApp:

.. role:: jsontype
.. role:: jsonkey

Writing Custom JSON-based Applications
======================================

A simple JSON-based application contains just a few lines of code:

.. code-block:: cpp

  #include <Tudat/JsonInterface/jsonInterface.h>

  int main( )
  {
    tudat::json_interface::JsonSimulationManager< > jsonSimulationManager( "main.json" );
    jsonSimulationManager.updateSettings( );
    jsonSimulationManager.runPropagation( );
    jsonSimulationManager.exportResults( );
    return EXIT_SUCCESS;
  }

which will always use the file :class:`main.json` in the current working directory as input file. If one wants to be able to run the application using different input files, then the code becomes:

.. code-block:: cpp
  :linenos:

  #include <Tudat/JsonInterface/jsonInterface.h>

  int main( int argc, char* argv[ ] )
  {
    const std::string inputPath = argv[ 1 ];
    tudat::json_interface::JsonSimulationManager< > jsonSimulationManager( inputPath );
    jsonSimulationManager.updateSettings( );
    jsonSimulationManager.runPropagation( );
    jsonSimulationManager.exportResults( );
    return EXIT_SUCCESS;
  }

which is basically the implementation of the :literal:`json_interface` application in :class:`Tudat/JsonInterface/jsonInterface.cpp`.

The four basic steps in a JSON-based applications are:

  1. Create a :class:`JsonSimulationManager` object (line 6). This class takes two template arguments, the first being the type to be used for the independent variable (the epoch) and the second the scalar type for the state variable (position, velocity, mass, etc.). When these types are not specified, the default values (:class:`double`) are used, i.e. :literal:`JsonSimulationManager< >` is equivalent to :literal:`JsonSimulationManager< double, double >`. The input argument is the absolute or relative path to a JSON file.

  .. note:: When creating a :class:`JsonSimulationManager` using this constructor, the current working directory is set to the directory in which the specified file is located.
  
  2. Set up the simulation (line 7). This uses the information contained in the JSON objects parsed from the specified input file to create all the settings objects necessary for the simulation (integrator, bodies, propagator, etc.).
  
  3. Run the propagation (line 8). In addition to integrating the equations of motion, calling the method :literal:`run` this will also check whether there are unused keys, print messages and/or generate a file with all the settings that are going be used depending on the settings specified in the key :jsonkey:`options` of the JSON object.

  4. Export the results (line 9). Calling the method :literal:`exportResults` will generate the output files with the requested results specified in the key :jsonkey:`export` of the JSON object.
  
Since not all Tudat features are supported by the JSON interface, in some cases it will be necessary to perform some additional steps between step 2 and 3, i.e. the simulation will be first set up from a JSON file, then additional settings will be provided manually, and the the simulation will be run. Then, optionally, the results can be exported using the :literal:`exportResults` method. However, defining the additional settings between lines 7 and 8 does not always lead to the desired results, since step 2 is a complex process in which some variables depend on each another.

Imagine that the aerodynamic coefficients of the body to propagate are not specified in the JSON file because they will be provided manually in the C++ file. This may result in an error, since trying to set up a simulation in which one wants to save e.g. the aerodynamic drag on a body with no aerodynamic coefficients interface will result in an error, which is generated while setting up the simulation and not when actually running it. Using a placeholder aerodynamic coefficients interface in the JSON file (e.g. zero force coefficients), this will result in a successful simulation set up, but then, when one wants to reset the aerodynamic coefficients interface in the C++ application manually, it is necessary to reset also all objects that had been set up based on the placeholder coefficients (e.g. the vehicle's flight conditions). This can be complex and is prone to leading to errors (sometimes run-time errors, sometimes successful simulations in which the results are wrong), so a different approach is followed when writing custom C++ applications.

When setting up the simulation (step 2), the following virtual method is called:

.. code-block:: cpp
  :linenos:
  
  virtual void updateSettings( )
  {
      ...
      resetIntegratorSettings( );
      resetSpice( );
      resetBodies( );              // must be called after resetIntegratorSettings and resetSpice
      resetExportSettings( );
      resetPropagatorSettings( );  // must be called after resetBodies and resetExportSettings
      resetApplicationOptions( );
      resetDynamicsSimulator( );
  }

Note that each of the methods called by this method is also virtual. This means that a derived class of :class:`JsonSimulationManager` can be created if a custom implementation of any of these methods is needed. The method :literal:`updateSettingsFromJSONObject` is generally not overridden, as it is dangerous to modify the order in which each of the virtual methods is called. If one *does* want to modify this method, the following has to be taken into account:

  - The default implementation of :literal:`resetBodies` uses the integrator settings' initial time to interpolate the ephemeris of celestial bodies if Spice is enabled and the key :jsonkey:`spice.preloadEpehemeris` is set to :literal:`true`.
  
  - The default implementation of :literal:`resetPropagatorSettings` uses the integrator settings' initial time to infer the initial state of the celestial bodies from their ephemeris, as well as some of the properties of other bodies (such as initial translational state or mass). Additionally, the variable to be computed are determined from the export settings.

In practice, this means that :literal:`resetBodies` must be called after :literal:`resetIntegratorSettings` and :literal:`resetSpice`, and :literal:`resetPropagatorSettings` must be called after :literal:`resetBodies` and :literal:`resetExportSettings`.

Thus, to avoid undefined behaviour, rather than overriding the :literal:`updateSettings` method, one would override just some of the virtual methods called therein. It is recommended to call the original implementation inside the custom implementations of these methods, and then provide additional steps. For instance, before the :literal:`main` function:

.. code-block:: cpp
  :linenos:
  
  class CustomJsonSimulationManager : public tudat::json_interface::JsonSimulationManager< >
  {
  public:
      // Inherit constructor.
      using JsonSimulationManager< >::JsonSimulationManager;

  protected:
      // Override resetBodies method.
      virtual void resetBodies( )
      {
          // First, call the original resetBodies, which uses the information in the JSON file.
          JsonSimulationManager::resetBodies( );

          // Then, provide additional steps.
          ...
      }
  };

Then, in the :literal:`main` function, we only need to change line 6 to:

.. code-block:: cpp
  
  CustomJsonSimulationManager jsonSimulationManager( inputPath );

If one wants to perform some operations on the results of the integration before exporting them, or does not want to export them to an output file, the call to the :literal:`exportResults` methods can be omitted, and the results can be retrieved from:

.. code-block:: cpp
  
  std::map< double, Eigen::VectorXd > stateHistory =
      jsonSimulationManager.getDynamicsSimulator( )->getEquationsOfMotionNumericalSolution( );

  std::map< double, Eigen::VectorXd > dependentVariablesHistory =
      jsonSimulationManager.getDynamicsSimulator( )->getDependentVariableHistory( );


A tutorial on how to write a custom JSON-based application can be found in :ref:`jsonInterface_tutorials_apolloCapsuleEntry`.
