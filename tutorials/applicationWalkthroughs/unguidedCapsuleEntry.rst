.. _walkthroughsUnguidedCapsuleEntry:

Un-guided Capsule Entry
=======================
The example described on this page is that of Apollo on a reentry trajectory towards the surface of Earth. The code for this tutorial is located in your Tudat Bundle at::

   tudatBundle/tudatExampleApplications/satellitePropagatorExamples/SatellitePropagatorExamples/apolloCapsuleEntry.cpp

For this example, we have the following problem statement:

*Given the position and velocity of the Apollo capsule at a certain point in time with respect to the Earth, what will its position and velocity be once it reaches an altitude of 25 km over the surface of Earth?*

.. warning:: The example described in this page assumes that the user has read the :ref:`walkthroughsUnperturbedEarthOrbitingSatellite`. This page only describes the differences with respect to such example, so please go back before proceding.

Create the vehicle
~~~~~~~~~~~~~~~~~~
First, the vehicle is created by placing an :literal:`"Apollo"` entry in the :literal:`bodyMap`, as shown below:

.. code-block:: cpp

    // Create vehicle objects.
    bodyMap[ "Apollo" ] = boost::make_shared< simulation_setup::Body >( );

Once that is done, an :class:`AerodynamicCoefficientInterface` is linked to the vehicle by means of the following command:

.. code-block:: cpp

    // Create vehicle aerodynamic coefficients
    bodyMap[ "Apollo" ]->setAerodynamicCoefficientInterface(
                unit_tests::getApolloCoefficientInterface( ) );

Note that in this case a pre-made interface that is part of the following Unit Tests is used::

    /tudatBundle/tudat/Tudat/Astrodynamics/Aerodynamics/UnitTests/testApolloCapsuleCoefficients.h

In this case, the :literal:`getApolloCoefficientInterface( )` function returns a :class:`HypersonicLocalInclinationAnalysis` object, which is a derived-class from :class:`AerodynamicCoefficientInterface`. Such object computes the aerodynamic coefficients during propagation by means of a local-inclination method. Finally, the mass of the Apollo capsule is set and the body creation is finalized:

.. code-block:: cpp

    bodyMap[ "Apollo" ]->setConstantBodyMass( 5.0E3 );

    // Finalize body creation.
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "J2000" );

.. tip:: If you want to include a different vehicle, you will have to define a new :class:`AerodynamicCoefficientInterface` and implement your own custom aerodynamic database. Please go to :ref:`tudatFeaturesAerodynamicGuidance` for further details.

Set up the acceleration models
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
A major difference with respect to the :ref:`walkthroughsUnperturbedEarthOrbitingSatellite` is the use of a spherical-harmonic gravity model and the presence of an aerodynamic force on the vehicle. Such acceleration models are added to the :literal:`accelerationMap` as follows:

.. code-block:: cpp

    // Define acceleration model settings.
    std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > accelerationsOfApollo;
    accelerationsOfApollo[ "Earth" ].push_back( boost::make_shared< SphericalHarmonicAccelerationSettings >( 4, 0 ) );
    accelerationsOfApollo[ "Earth" ].push_back( boost::make_shared< AccelerationSettings >( aerodynamic ) );
    accelerationMap[  "Apollo" ] = accelerationsOfApollo;

A crucial step in reentry modelling is the definition of a :class:`AerodynamicGuidance` model. Controlling the orientation of the vehicle during atmospheric flight plays an important role in the shape of the trajectory as well as on the magnitude of the aerodynamic and thermal loads. In this example, a simple fixed-angle aerodynamic guidance is used:

.. code-block:: cpp

    // Define constant 30 degree angle of attack
    double constantAngleOfAttack = 30.0 * mathematical_constants::PI / 180.0;
    bodyMap.at( "Apollo" )->getFlightConditions( )->getAerodynamicAngleCalculator( )->setOrientationAngleFunctions(
                boost::lambda::constant( constantAngleOfAttack ) );

Set up the propagation settings
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
In most reentry studies, it is convenient to define the entry conditions using a spherical state. The following entry state is used:

- Altitude: 120 km
- Latitude: 0 deg
- Longitude: 68.75 deg 
- Inertial speed: 7.7 km/s
- Flight-path angle: -0.9 deg
- Heading angle: 34.37 deg

Such state must is defined and converted to Cartesian state variables as follows:

.. code-block:: cpp

    // Set spherical elements for Apollo.
    Eigen::Vector6d apolloSphericalEntryState;
    apolloSphericalEntryState( SphericalOrbitalStateElementIndices::radiusIndex ) =
            spice_interface::getAverageRadius( "Earth" ) + 120.0E3;
    apolloSphericalEntryState( SphericalOrbitalStateElementIndices::latitudeIndex ) = 0.0;
    apolloSphericalEntryState( SphericalOrbitalStateElementIndices::longitudeIndex ) = 1.2;
    apolloSphericalEntryState( SphericalOrbitalStateElementIndices::speedIndex ) = 7.7E3;
    apolloSphericalEntryState( SphericalOrbitalStateElementIndices::flightPathIndex ) =
            -0.9 * mathematical_constants::PI / 180.0;
    apolloSphericalEntryState( SphericalOrbitalStateElementIndices::headingAngleIndex ) = 0.6;

    // Convert apollo state from spherical elements to Cartesian elements.
    Eigen::Vector6d systemInitialState = convertSphericalOrbitalToCartesianState(
                apolloSphericalEntryState );

.. note:: Note that speed defined in the :literal:`speedIndex` makes reference to the **inertial** speed of the vehicle. Furthermore, the :literal:`latitudeIndex` makes reference to the **geocentric** latitude.

Create a list of dependent variables to save
********************************************
In this example, a number of dependent variables are saved to plot the trajectory of Apollo after reentry. The following dependent variables are saved:

- Mach number
- Altitude
- Aerodynamic acceleration norm
- Aerodynamic force coefficients (CD, CS, CL)

First, a :literal:`dependentVariablesList` needs to be created, which will list all the variables to save:

.. code-block:: cpp

    // Define list of dependent variables to save.
    std::vector< boost::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariablesList;

Next, the list is populated with the desired dependent variables. Please go to :ref:`tudatFeaturesPropagatorSettingsDependentVariables` for further details on the various dependent variables that an be stored:

.. code-block:: cpp

    dependentVariablesList.push_back(
                boost::make_shared< SingleDependentVariableSaveSettings >( mach_number_dependent_variable, "Apollo" ) );
    dependentVariablesList.push_back(
                boost::make_shared< SingleDependentVariableSaveSettings >(
                    altitude_dependent_variable, "Apollo", "Earth" ) );
    dependentVariablesList.push_back(
                boost::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                    aerodynamic, "Apollo", "Earth", 1 ) );
    dependentVariablesList.push_back(
                boost::make_shared< SingleDependentVariableSaveSettings >(
                    aerodynamic_force_coefficients_dependent_variable, "Apollo" ) );

    // Create object with list of dependent variables
    boost::shared_ptr< DependentVariableSaveSettings > dependentVariablesToSave =
            boost::make_shared< DependentVariableSaveSettings >( dependentVariablesList );

Define the termination conditions
*********************************
Finally, the termination conditions are established. In this example, the reentry trajectory is propagated until Apollo's altitude drops below 25 km:

.. code-block:: cpp

    // Define termination conditions
    boost::shared_ptr< SingleDependentVariableSaveSettings > terminationDependentVariable =
            boost::make_shared< SingleDependentVariableSaveSettings >(
                altitude_dependent_variable, "Apollo", "Earth" );
    boost::shared_ptr< PropagationTerminationSettings > terminationSettings =
            boost::make_shared< PropagationDependentVariableTerminationSettings >(
                terminationDependentVariable, 25.0E3, true );

Please go to :ref:`tudatFeaturesPropagatorSettingsTermination` for a detailed description of the available termination conditions.

