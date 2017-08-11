.. _tudatFeaturesPropagatorSettingsDependentVariables:

Propagator Settings: Dependent Variables
========================================
By default, the :class:`DynamicsSimulator` propagates the state of a body which uniquely defines its position and velocity. However, it is possible to define a number of dependent variables that are derived from such state. These dependent variables can be used as termination conditions as discussed in :ref:`tudatFeaturesPropagatorSettingsTermination` or saved for further post-processing. This page describes the dependent variables currently available and how these are retrieved from Tudat.

Saving dependent variables
~~~~~~~~~~~~~~~~~~~~~~~~~~
The dependent variables are computed during the body propagation, thus the user must provide a list of dependent variables to save prior creating the :literal:`DynamicsSimulator`:

.. code-block:: cpp

      // Define list of dependent variables to save.
      std::vector< boost::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariablesList;

where :literal:`dependentVariableList` is populated as follows:

.. code-block:: cpp

      dependentVariablesList.push_back(
                boost::make_shared< SingleDependentVariableSaveSettings >( variableType , associatedBody , secondaryBody ) );

where:

- :class:`SingleDependentVariableSaveSettings`

   Defines the derived-class being used and must match with :literal:`variableType`.

- :literal:`variableType`

   Indicates the variable type and must match with one of the :literal:`enum` values in :class:`PropagationDependentVariables`. A detailed description of how is done is given in the next section.

- :literal:`associatedBody`

   Indicates to which body the saved dependent variables are associated.

- :literal:`secondaryBody`

   Optional argument that provides a secondary body that may be necessary to save the dependent variable. By default, this argument is empty.

Once the list of dependent variables to save has been populated, a :literal:`boost::shared_ptr< DependentVariableSaveSettings >` object needs to be created and passed to :literal:`propagatorSettings`:

.. code-block:: cpp

      // Create object with list of dependent variables
      boost::shared_ptr< DependentVariableSaveSettings > dependentVariablesToSave =
            boost::make_shared< DependentVariableSaveSettings >( dependentVariablesList );

      // Create propagation settings.
      boost::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
            boost::make_shared< TranslationalStatePropagatorSettings< double > >
            ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState,
              terminationSettings, propagator, dependentVariablesToSave );

.. note:: In the example above, the :class:`TranslationalStatePropagatorSettings` derived-class is used. Please note that any of the derived-classes described in :ref:`tudatFeaturesPropagatorSettings` can be used, as long as these support dependent variable saving.


Available dependent variables
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The framework discussed in the previous section explains how the :literal:`dependentVariablesList` is populated and passed to the :class:`PropagatorSettings`. The goal of this section is to list the available dependent variables and to explain how these are pushed to the :literal:`dependentVariablesList`.

.. class:: SingleDependentVariableSaveSettings

   This base-class is a generic method to retrieve a large number of dependent variables that are not classified under a particular group. Variables are saved to the :literal:`dependentVariablesList` using the following code:

   .. code-block:: cpp

      dependentVariablesList.push_back(
                   boost::make_shared< SingleDependentVariableSaveSettings >( variableType, associatedBody, secondaryBody ) );

   where:

   - :literal:`variableType`

      :class:`PropagationDependentVariables` variable that can take the following values:

      **Variables returning a** :literal:`double`

      - :literal:`mach_number_dependent_variable`
      - :literal:`altitude_dependent_variable`
      - :literal:`airspeed_dependent_variable`
      - :literal:`local_density_dependent_variable`
      - :literal:`relative_speed_dependent_variable`
      - :literal:`relative_distance_dependent_variable`
      - :literal:`radiation_pressure_dependent_variable`
      - :literal:`total_aerodynamic_g_load_variable`
      - :literal:`stagnation_point_heat_flux_dependent_variable`
      - :literal:`local_temperature_dependent_variable`
      - :literal:`geodetic_latitude_dependent_variable`
      - :literal:`control_surface_deflection_dependent_variable`
      - :literal:`total_mass_rate_dependent_variables`
      - :literal:`periapsis_altitude_dependent_variable`

      **Variables returning an** :literal:`Eigen::VectorXd`

      - :literal:`relative_position_dependent_variable`
      - :literal:`relative_velocity_dependent_variable`
      - :literal:`body_fixed_airspeed_based_velocity_variable`
      - :literal:`total_acceleration_norm_dependent_variable`
      - :literal:`total_acceleration_dependent_variable`
      - :literal:`aerodynamic_force_coefficients_dependent_variable`
      - :literal:`aerodynamic_moment_coefficients_dependent_variable`
      - :literal:`lvlh_to_inertial_frame_rotation_dependent_variable`
      - :literal:`rotation_matrix_to_body_fixed_frame_variable`

.. class:: SingleAccelerationDependentVariableSaveSettings

   This derived-class is used to retrieve acceleration-related dependent variables. A large number of acceleration models are supported and both the acceleration-norm and the acceleration-vector can be saved. Variables are saved to the :literal:`dependentVariablesList` using the following code:

   .. code-block:: cpp

      dependentVariablesList.push_back(
                boost::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                    accelerationModeType, bodyUndergoingAcceleration, bodyExertingAcceleration, useNorm ) );

   where:

   - :literal:`accelerationModeType`
  
      :class:`AvailableAcceleration` variable that defines the type of acceleration that must be retrieved. It can take the following values:
      
      - :literal:`undefined_acceleration`
      - :literal:`central_gravity`
      - :literal:`aerodynamic`
      - :literal:`cannon_ball_radiation_pressure`
      - :literal:`spherical_harmonic_gravity`
      - :literal:`mutual_spherical_harmonic_gravity`
      - :literal:`third_body_central_gravity`
      - :literal:`third_body_spherical_harmonic_gravity`
      - :literal:`third_body_mutual_spherical_harmonic_gravity`
      - :literal:`thrust_acceleration`

   - :literal:`bodyUndergoingAcceleration`

      :literal:`std::string` variable that indicates the body that experiences the acceleration that needs to be retrieved. Make sure that the body's name is listed in :class:`BodyMap`.

   - :literal:`bodyExertingAcceleration`

      :literal:`std::string` variable that indicates the body that exerts the acceleration that needs to be retrieved on :literal:`bodyUndergoingAcceleration`. Make sure that the body's name is listed in :class:`BodyMap`.

   - :literal:`useNorm`

      :literal:`bool` variable that indicates if the norm of the acceleration (true) or the acceleration vector (false) must be retrieved.

   .. warning:: Make sure that the selected :literal:`bodyExertingAcceleration` is compatible with the :literal:`accelerationModeType`.

.. class:: IntermediateAerodynamicRotationVariableSaveSettings

   This derived-class is used to retrieve the rotation matrix between two desired frames. Variables are saved to the :literal:`dependentVariablesList` using the following code:

   .. code-block:: cpp

      dependentVariablesList.push_back(
                boost::make_shared< IntermediateAerodynamicRotationVariableSaveSettings >(
                    associatedBody, baseFrame, targetFrame ) );

   where:

   - :literal:`associatedBody`

      :literal:`std::string` variable that indicates the body for which a rotation matrix is to be saved. Make sure that the body's name is listed in :class:`BodyMap`.

   - :literal:`baseFrame`

      :class:`AerodynamicsReferenceFrames` variable indicates the frame from which the rotation is to be saved. The following frames are available:

      - :literal:`inertial_frame`
      - :literal:`corotating_frame`
      - :literal:`vertical_frame`
      - :literal:`trajectory_frame`
      - :literal:`aerodynamic_frame`
      - :literal:`body_frame`

   - :literal:`targetFrame`

      :class:`AerodynamicsReferenceFrames` variable indicates the frame to which the rotation is to be saved. The available frames are listed above.

.. class:: BodyAerodynamicAngleVariableSaveSettings

   This derived-class is used to retrieve a number of rotation angles. Variables are saved to the :literal:`dependentVariablesList` using the following code:

   .. code-block:: cpp

      dependentVariablesList.push_back(
                boost::make_shared< BodyAerodynamicAngleVariableSaveSettings >(
                    associatedBody, angle ) );

   where:

   - :literal:`associatedBody`

      :literal:`std::string` variable that indicates the body for which the :literal:`angle` is to be saved. Make sure that the body's name is listed in :class:`BodyMap`.

   - :literal:`angle`

      :class:`AerodynamicsReferenceFrameAngles` variable that provides the angle to be saved. The following angles can be saved using this method:

      - :literal:`latitude_angle`
      - :literal:`longitude_angle`
      - :literal:`heading_angle`
      - :literal:`flight_path_angle`
      - :literal:`angle_of_attack`
      - :literal:`angle_of_sideslip`
      - :literal:`bank_angle`

.. warning:: At the moment, all the multi-dimensional dependent variables are pushed to the end of the save-file. This is currently being fixed, such that the save order is the same as the declaration order. Please refer to the following Pull Request for further details: https://github.com/Tudat/tudat/pull/191
