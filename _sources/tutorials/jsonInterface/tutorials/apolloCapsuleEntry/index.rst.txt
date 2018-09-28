.. _jsonInterface_tutorials_apolloCapsuleEntry:

.. role:: jsontype
.. role:: jsonkey

Apollo Capsule Entry
====================

The example described on this page is that of Apollo on a reentry trajectory towards the surface of Earth. The code for this tutorial is located in your Tudat Bundle at: ::

  tudatBundle/tudatExampleApplications/satellitePropagatorExamples/SatellitePropagatorExamples/apolloCapsuleEntryJSON.cpp

In this tutorial, the creation of a custom JSON-based application is described. Read first the general introduction on :ref:`jsonInterface_customJsonApp` and the example :ref:`walkthroughsUnguidedCapsuleEntry`, which is identical to this one but without the use of JSON files.

The Apollo capsule entry cannot be simulated using exclusively JSON input files and the :literal:`json_interface` application. Some of the Tudat features that this application uses cannot be provided using JSON files, and thus a custom JSON-based application, in which part of the settings are read from a JSON file and others specified manually using C++ code, has to be written.

The Tudat features that this application uses that are not supported by the JSON interface are:

  - Setting a body's :class:`AerodynamicCoefficientInterface` directly instead of using an :class:`AerodynamicCoefficientSettings` object.
  - Using the function :literal:`tudat_applications::getOutputPath( )`.
  - Setting an angle of attack for a body.

This functionality has to be implemented in the C++ application. However, as for any JSON-based application, the first step is to generate the input file. This file is located in your Tudat Bundle at::

  tudatBundle/tudatExampleApplications/satellitePropagatorExamples/SatellitePropagatorExamples/apolloCapsuleEntry.json

A few remarks about the contents of this file can be made:

  - The key :jsonkey:`spice.preloadEphemeris` has to be set to :literal:`false` because the key :jsonkey:`finalEpoch` is not defined, so it is not possible to interpolate the ephemeris of celestial bodies.
  - For the body :literal:`Apollo`, only the mass is specified; the aerodynamic settings (including coefficients and reference area) will be loaded manually in the C++ application.
  - Since the key :jsonkey:`file` of the export settings is mandatory, a placeholder empty string is provided for :jsonkey:`export[0].file` and :jsonkey:`export[1].file`, which will be then redefined in the C++ application.
  - The variable representing Apollo's altitude w.r.t. Earth is used twice: as termination condition and as dependent variable to be exported to an output file. Rather than repeating it twice, the key :jsonkey:`export[1].variables[1]` has been set to reference the termination's variable by using the special string :literal:`"${termination.variable}"`.

Then, in the C++ app, as explained in :ref:`jsonInterface_customJsonApp`, a derived class of :class:`JsonSimulationManager` has been created, providing custom implementations for the virtual methods :literal:`resetBodies`, :literal:`resetExportSettings` and :literal:`resetPropagatorSettings`:

.. code-block:: cpp

  class ApolloJsonSimulationManager : public tudat::json_interface::JsonSimulationManager< >
  {
  public:

      // Inherit constructor.
      using JsonSimulationManager< >::JsonSimulationManager;

  protected:

      virtual void resetBodies( )
      {
          ...
      }

      virtual void resetExportSettings( )
      {
          ...
      }

      virtual void resetPropagatorSettings( )
      {
          ...
      }
      
  };

For the :literal:`resetBodies` method, we first call the default implementation of the base class, and then we create Apollo's aerodynamic coefficients interface manually, using a function defined in :literal:`tudat::unit_tests`:

.. code-block:: cpp

  virtual void resetBodies( )
  {
      // First, call the original resetBodies, which uses the information in the JSON file
      JsonSimulationManager::resetBodies( );

      // Then, create vehicle's aerodynamic coefficients interface
      getBody( "Apollo" )->setAerodynamicCoefficientInterface( tudat::unit_tests::getApolloCoefficientInterface( ) );
  }

Note that we use the method :literal:`getBody( const std::string& bodyName )` of :class:`JsonSimulationManager`, which retrieves the body named :literal:`bodyName` from the body map. If we want to create a new body manually, we use the method :literal:`addBody( const std::string& bodyName )` (note that, if the body already exists, it will be reset to a body created with the empty constructor).

For the :literal:`resetExportSettings` method, we first call the default implementation of the base class, and then modify the output files paths using the function :literal:`tudat_applications::getOutputPath( )`:

.. code-block:: cpp

  virtual void resetExportSettings( )
  {
      // First, call the original resetExportSettings, which uses the information in the JSON file
      JsonSimulationManager::resetExportSettings( );

      // Then, replace the output file paths (empty strings placeholders had been specified in the JSON file)
      const std::string outputDirectory = tudat_applications::getOutputPath( ) + "ApolloCapsuleExampleJSON/";
      getExportSettings( 0 )->setOutputFile( outputDirectory + "apolloPropagationHistory.dat" );
      getExportSettings( 1 )->setOutputFile( outputDirectory + "apolloDependentVariableHistory.dat" );
  }

In this case, we use the method :literal:`getExportSettings( unsigned int index )`, which returns the :class:`ExportSettings` object created from the key :jsonkey:`export[index]` in the JSON file.

Finally, we also override the method :literal:`resetPropagatorSettings`:

.. code-block:: cpp
  :linenos:

  virtual void resetPropagatorSettings( )
  {
      // First, call the original resetPropagatorSettings, which uses the information in the JSON file
      JsonSimulationManager::resetPropagatorSettings( );

      // Define constant 30 degree angle of attack
      double constantAngleOfAttack = 30.0 * tudat::mathematical_constants::PI / 180.0;
      getBody( "Apollo" )->getFlightConditions( )->getAerodynamicAngleCalculator( )->
              setOrientationAngleFunctions( boost::lambda::constant( constantAngleOfAttack ) );
  }

After calling the original implementation in line 4, we define an angle of attack for Apollo. This is done by modifying Apollo's flight conditions. This can only be done once the flight conditions have been created, i.e. after :literal:`JsonSimulationManager::resetPropagatorSettings( )` has been called, since this method will create the acceleration aerodynamic models based on the body map created previously.

Now, we can write our :literal:`main` function, which uses the custom class :class:`ApolloJsonSimulationManager`:

.. code-block:: cpp

  int main( )
  {
      const std::string cppFilePath( __FILE__ );
      const std::string cppFolder = cppFilePath.substr( 0, cppFilePath.find_last_of("/\\") + 1 );

      ApolloJsonSimulationManager jsonSimulationManager( cppFolder + "apolloCapsuleEntry.json" );
      jsonSimulationManager.updateSettings( );
      jsonSimulationManager.runPropagation( );
      jsonSimulationManager.exportResults( );

      return EXIT_SUCCESS;
  }

After running the application, the results can be found in the directory::

  tudatBundle/tudatExampleApplications/satellitePropagatorExamples/SatellitePropagatorExamples/SimulationOutput/ApolloCapsuleExampleJSON
