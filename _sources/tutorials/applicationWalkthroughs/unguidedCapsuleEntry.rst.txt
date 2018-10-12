.. _walkthroughsUnguidedCapsuleEntry:

Un-guided Capsule Entry
=======================
The example described on this page is that of Apollo on a re-entry trajectory towards the surface of Earth. The code for this tutorial is given on Github, and is also located in your tudat bundle at::

   tudatBundle/tudatExampleApplications/satellitePropagatorExamples/SatellitePropagatorExamples/apolloCapsuleEntry.cpp

For this example, we have the following problem statement:

   *Given the position and velocity of the Apollo capsule at a certain point in time with respect to the Earth, what will its position and velocity be once it reaches an altitude of 25 km over the surface of Earth?*

.. warning:: The example described in this page assumes that the user has read the :ref:`walkthroughsUnperturbedEarthOrbitingSatellite`. This page only describes the differences with respect to such example, so please go back before proceeding.

Create the Vehicle
~~~~~~~~~~~~~~~~~~
First, the vehicle is created by placing an :literal:`"Apollo"` entry vehicle in the :literal:`bodyMap`, as shown below:

.. code-block:: cpp

    // Create vehicle objects.
    bodyMap[ "Apollo" ] = std::make_shared< simulation_setup::Body >( );

Once that is done, an :class:`AerodynamicCoefficientInterface` is linked to the vehicle by means of the following command:

.. code-block:: cpp

    // Create vehicle aerodynamic coefficients
    bodyMap[ "Apollo" ]->setAerodynamicCoefficientInterface(
                unit_tests::getApolloCoefficientInterface( ) );

In this example a pre-made interface is used to define the aerodynamic coefficients. The :literal:`getApolloCoefficientInterface( )` function returns a :class:`HypersonicLocalInclinationAnalysis` object, which is a derived-class from :class:`AerodynamicCoefficientInterface`. It computes the aerodynamic coefficients during propagation by means of a local-inclination method. Finally, the mass of the Apollo capsule is set and the body creation is finalized:

.. code-block:: cpp

    bodyMap[ "Apollo" ]->setConstantBodyMass( 5.0E3 );

    // Finalize body creation.
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "J2000" );

.. tip:: If you want to include a different vehicle, you will have to define a new :class:`AerodynamicCoefficientInterface` and implement your own custom aerodynamic database. Please go to :ref:`tudatFeaturesAerodynamicGuidance` for further details.

Set Up the Acceleration Models
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
A major difference with respect to the :ref:`walkthroughsUnperturbedEarthOrbitingSatellite` is the use of a spherical-harmonic gravity model and the presence of an aerodynamic force on the vehicle. The spherical-harmonic gravity model is selected by the derived-class :class:`SphericalHarmonicAccelerationSettings` with degree and order as input parameters. Both acceleration models are added to the :literal:`accelerationMap` as follows:

.. code-block:: cpp

    // Define acceleration model settings.
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfApollo;
    accelerationsOfApollo[ "Earth" ].push_back( std::make_shared< SphericalHarmonicAccelerationSettings >( 4, 0 ) );
    accelerationsOfApollo[ "Earth" ].push_back( std::make_shared< AccelerationSettings >( aerodynamic ) );
    accelerationMap[ "Apollo" ] = accelerationsOfApollo;

A crucial step in re-entry modelling is the definition of an :class:`AerodynamicGuidance` model. Controlling the orientation of the vehicle during atmospheric flight plays an important role in the shape of the trajectory as well as on the magnitude of the aerodynamic and thermal loads. In this example, a simple fixed-angle aerodynamic guidance model is used. This is implemented using a lambda expression (explained in detail :ref:`here <exampleUtilityLambdaExpression>`). In short this function always outputs the value of :literal:`constantAngleOfAttack`, which in turn sets the orientation angles of the :literal:`"Apollo"` body:

.. code-block:: cpp

    // Define constant 30 degree angle of attack
    double constantAngleOfAttack = unit_conversions::convertDegreesToRadians(30.0);
    bodyMap.at( "Apollo" )->getFlightConditions( )->getAerodynamicAngleCalculator( )->setOrientationAngleFunctions(
                [ = ]( ){ return constantAngleOfAttack; } );

.. tip:: To view the available options for aerodynamic guidance check out the :ref:`tudatFeaturesAerodynamicGuidance` section. 

Set Up the Propagation Settings
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
In most re-entry studies, it is convenient to define the entry conditions using a spherical state. The following entry state is used:

   - Altitude: 120 km
   - Latitude: 0 deg
   - Longitude: 68.75 deg 
   - Inertial speed: 7.7 km/s
   - Flight-path angle: -0.9 deg
   - Heading angle: 34.37 deg

Such state must be defined and converted to Cartesian state variables as follows:

.. code-block:: cpp

    // Set spherical elements for Apollo.
    Eigen::Vector6d apolloSphericalEntryState;
    apolloSphericalEntryState( SphericalOrbitalStateElementIndices::radiusIndex ) =
            spice_interface::getAverageRadius( "Earth" ) + 120.0E3;
    apolloSphericalEntryState( SphericalOrbitalStateElementIndices::latitudeIndex ) =
            unit_conversions::convertDegreesToRadians( 0.0 );
    apolloSphericalEntryState( SphericalOrbitalStateElementIndices::longitudeIndex ) = 
            unit_conversions::convertDegreesToRadians( 68.75 );
    apolloSphericalEntryState( SphericalOrbitalStateElementIndices::speedIndex ) = 7.7E3;
    apolloSphericalEntryState( SphericalOrbitalStateElementIndices::flightPathIndex ) =
            unit_conversions::convertDegreesToRadians( -0.9 );
    apolloSphericalEntryState( SphericalOrbitalStateElementIndices::headingAngleIndex ) = 
            unit_conversions::convertDegreesToRadians( 34.37 );

    // Convert apollo state from spherical elements to Cartesian elements.
    Eigen::Vector6d systemInitialState = convertSphericalOrbitalToCartesianState(
                apolloSphericalEntryState );

Finally, the state needs to be converted to inertial frame. As it is customary, in fact, the spherical state above (and thus also the converted Cartesian state), is defined in the rotating frame centered at the Earth center (ECEF). However, as mentioned in :ref:`walkthroughsUnperturbedEarthOrbitingSatellite`, in Tudat propagation can only occurr in inertial frames. Hence, the following lines are included:

.. code-block:: cpp
   
    // Convert the state to the global (inertial) frame.
    std::shared_ptr< ephemerides::RotationalEphemeris > earthRotationalEphemeris =
            bodyMap.at( "Earth" )->getRotationalEphemeris( );
    systemInitialState = transformStateToGlobalFrame( systemInitialState, simulationStartEpoch, earthRotationalEphemeris );

.. note:: Note that the :literal:`latitudeIndex` of :literal:`apolloSphericalEntryState` makes reference to the **geocentric** latitude.

Create a List of Dependent Variables to Save
********************************************
In this example, a number of dependent variables are saved to plot the trajectory of Apollo after re-entry. The following dependent variables are saved:

   - Mach number
   - Altitude
   - Aerodynamic acceleration norm
   - Aerodynamic force coefficients (:math:`C_D`, :math:`C_S`, :math:`C_L`)

First, a :literal:`dependentVariablesList` needs to be created, which will list all the variables to save:

.. code-block:: cpp

    // Define list of dependent variables to save.
    std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariablesList;

Next, the list is populated with the desired dependent variables. Please go to :ref:`tudatFeaturesPropagatorSettingsDependentVariables` for further details on the various dependent variables that can be stored:

.. code-block:: cpp

    dependentVariablesList.push_back( std::make_shared< SingleDependentVariableSaveSettings >(
                                          mach_number_dependent_variable, "Apollo" ) );
    dependentVariablesList.push_back( std::make_shared< SingleDependentVariableSaveSettings >(
                                          altitude_dependent_variable, "Apollo", "Earth" ) );
    dependentVariablesList.push_back( std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                                          aerodynamic, "Apollo", "Earth", 1 ) );
    dependentVariablesList.push_back( std::make_shared< SingleDependentVariableSaveSettings >(
                                          aerodynamic_force_coefficients_dependent_variable, "Apollo" ) );

    // Create object with list of dependent variables
    std::shared_ptr< DependentVariableSaveSettings > dependentVariablesToSave =
            std::make_shared< DependentVariableSaveSettings >( dependentVariablesList );

.. tip:: If you do not want the console output generated by :class:`DependentVariableSaveSettings` to show, you can simply set the second input argument as :literal:`false`.

Define the Termination Conditions
*********************************
Finally, the termination conditions are established. The termination settings are stored in the :class:`PropagationTerminationSettings` object. In this example, the re-entry trajectory is propagated until Apollo's altitude drops below 25 km. The boolean in the constructor of the derived-class :class:`PropagationDependentVariableTerminationSettings` indicates whether the simulation is terminated when :literal:`terminationDependentVariable` goes below the supplied value (:literal:`true`) or above (:literal:`false`):

.. code-block:: cpp

    // Define termination conditions
    std::shared_ptr< SingleDependentVariableSaveSettings > terminationDependentVariable =
            std::make_shared< SingleDependentVariableSaveSettings >( altitude_dependent_variable, "Apollo", "Earth" );
    std::shared_ptr< PropagationTerminationSettings > terminationSettings =
            std::make_shared< PropagationDependentVariableTerminationSettings >( terminationDependentVariable, 25.0E3, true );

.. tip:: Please go to :ref:`tudatFeaturesPropagatorSettingsTermination` for a detailed description of the available termination conditions.

Results
~~~~~~~

Below the history of some of the saved parameters is shown. One can see the capsule skipping several times before its final descent into the atmosphere, where it reaches 25 km altitude. The dependent variable history can be obtained in Tudat from the :literal:`getDependentVariableHistory` function of the :class:`DynamicsSimulator` class. The resulting :literal:`std::map` can be saved as discussed in :ref:`tudatFeaturesInputOutput`. 

.. figure:: images/apolloResults.png

.. tip:: Open the figure in a new tab for more detail.