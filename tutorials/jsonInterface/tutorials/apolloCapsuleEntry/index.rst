.. _jsonInterface_tutorials_apolloCapsuleEntry:

.. role:: jsontype
.. role:: jsonkey

Apollo Capsule Entry
====================

The example described on this page is that of Apollo on a reentry trajectory towards the surface of Earth. The code for this tutorial is located in your Tudat Bundle at:

.. code-block:: txt

  tudatBundle/tudatExampleApplications/satellitePropagatorExamples/SatellitePropagatorExamples/apolloCapsuleEntryJSON.cpp

In this tutorial, the creation of a custom JSON-based application is described. Read first the general introduction on :ref:`jsonInterface_customJsonApp` and the example :ref:`walkthroughsUnguidedCapsuleEntry`, which is identical to this one but without the use of JSON files.

The Apollo capsule entry cannot be simulated using exclusively JSON input files and the :literal:`json_interface` application. Some of the Tudat features that this application uses cannot be provided using JSON files, and thus a custom JSON-based application, in which part of the settings are read from a JSON file and others specified manually using C++ code, has to be written.

The Tudat features that this application uses that are not supported by the JSON interface are:

  - Setting a body's :class:`AerodynamicCoefficientInterface` directly instead of using an :class:`AerodynamicCoefficientSettings` object.
  - Using the function :literal:`tudat_applications::getOutputPath( )`.
  - Settings an angle of attach for a body.
  - Defining the initial state using spherical elements.

This functionality has to be implemented in the C++ application. However, as for any JSON-based application, the first step is to generate the input file. This file is located in your Tudat Bundle at:

.. code-block:: txt

  tudatBundle/tudatExampleApplications/satellitePropagatorExamples/SatellitePropagatorExamples/apolloCapsuleEntry.json

A few remarks about the contents of this file can be made:

  - The key :jsonkey:`spice.preloadEphemeris` has to be set to :literal:`false` because the key :jsonkey:`finalEpoch` is not defined, so it is not possible to interpolate the ephemeris of celestial bodies.
  - For the body :literal:`Apollo`, only the mass is specified; the aerodynamic settings (including coefficients and reference area) will be loaded manually in the C++ application.
  - If the initial state is not defined, the JSON validator will throw an :class:`UndefinedKeyError`. Thus, a placeholder initial state, which will be modified in the C++ app, is used for key :jsonkey:`propagators[0].initialStates`, which is set to :literal:`[0, 0, 0, 0, 0, 0]`.
  - In the same way, since the key :jsonkey:`file` of the export settings is mandatory, a placeholder empty string is provided for :jsonkey:`export[0].file` and :jsonkey:`export[1].file`, which will be then redefined in the C++ application.
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

      // Set spherical elements for Apollo.
      using namespace tudat::orbital_element_conversions;
      Eigen::Vector6d apolloSphericalEntryState;
      apolloSphericalEntryState( SphericalOrbitalStateElementIndices::radiusIndex ) =
              getBody( "Earth" )->getShapeModel( )->getAverageRadius( ) + 120.0E3;
      apolloSphericalEntryState( SphericalOrbitalStateElementIndices::latitudeIndex ) = 0.0;
      apolloSphericalEntryState( SphericalOrbitalStateElementIndices::longitudeIndex ) = 1.2;
      apolloSphericalEntryState( SphericalOrbitalStateElementIndices::speedIndex ) = 7.7E3;
      apolloSphericalEntryState( SphericalOrbitalStateElementIndices::flightPathIndex ) =
              -0.9 * tudat::mathematical_constants::PI / 180.0;
      apolloSphericalEntryState( SphericalOrbitalStateElementIndices::headingAngleIndex ) = 0.6;

      // Convert apollo state from spherical elements to Cartesian elements.
      Eigen::Vector6d systemInitialState = convertSphericalOrbitalToCartesianState( apolloSphericalEntryState );
      systemInitialState = tudat::ephemerides::transformStateToGlobalFrame(
                  systemInitialState, getStartEpoch( ), getBody( "Earth" )->getRotationalEphemeris( ) );

      // Reset initial states (zero placeholder vector had been specified in the JSON file)
      getPropagatorSettings( )->resetInitialStates( systemInitialState );
  }

After calling the original implementation in line 4 (this won't result in a JSON validation error because the the key :jsonkey:`propagators[0].initialStates` was defined with a placeholder value), we perform two tasks: defining an angle of attach for Apollo, in lines 9-11; and resetting the propagator's initial states in lines 13-30.

The angle of attack is defined by modifying Apollo's flight conditions. This can only be done once the flight conditions have been created, i.e. after :literal:`JsonSimulationManager::resetPropagatorSettings( )` has been called, since this method will create the acceleration aerodynamic models based on the body map created previously.

Then, we define the spherical initial state of Apollo, convert it to Cartesian components and transform it to the global frame, using Earth's rotational ephemeris object, and the simulation start epoch, which can be retrieved using the method :literal:`getStartEpoch( )`. To get the simulation end epoch, the method :literal:`getEndEpoch( )` can be used, which will return :literal:`TUDAT_NAN` if no time termination condition has been defined. Once we have created the combined vector of initial states, :literal:`systemInitialState`, which must include all the states of all the bodies and propagators, we reset the propagator's initial states in line 30. Note that :literal:`getPropagatorSettings( )` always return a :class:`MultiTypePropagatorSettings`, even when only one propagator is used (in which case it will only contain one :class:`SingleArcPropagatorSettings`). The method :literal:`resetInitialStates` of :class:`MultiTypePropagatorSettings` will reset the states of the *children* single-arc propagators, as well as the combined states of the *parent* multi-type propagator.

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

