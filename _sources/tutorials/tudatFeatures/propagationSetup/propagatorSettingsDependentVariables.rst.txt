.. _tudatFeaturesPropagatorSettingsDependentVariables:

Propagator Settings: Dependent Variables
========================================
By default, the :class:`DynamicsSimulator` propagates the state of a body which uniquely defines its position and velocity. However, it is possible to define a number of dependent variables that are derived from such state. These dependent variables can be used as termination conditions as discussed in :ref:`tudatFeaturesPropagatorSettingsTermination` or saved for further post-processing. This page describes the dependent variables currently available and how these are retrieved from Tudat.

The general procedure consists of defining an object of class :class:`DependentVariableSaveSettings`, and then feeding it to a derived class of :class:`PropagatorSettings`.

.. class:: DependentVariableSaveSettings

   .. code-block:: cpp

      // Create object with list of dependent variables
      std::shared_ptr< DependentVariableSaveSettings > dependentVariablesToSave =
      std::make_shared< DependentVariableSaveSettings >( dependentVariablesList );

      // Create propagation settings.
      std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
      std::make_shared< TranslationalStatePropagatorSettings< double > >
      ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState,
        terminationSettings, propagator, dependentVariablesToSave );

   This is the general class describing the dependent variables to be saved during propagation. Its only input is a list of dependent variables, which can be of many different types. In the sections below, you will find a thorough description of what can be added to said list. 

   .. note:: In the example above, the :class:`TranslationalStatePropagatorSettings` derived class is used. Please note that any of the derived classes described in :ref:`tudatFeaturesPropagatorSettings` can be used, as long as these support dependent variable saving.

Saving Dependent Variables
~~~~~~~~~~~~~~~~~~~~~~~~~~
The dependent variables are computed during the body propagation, thus the user must provide a list of dependent variables to save prior creating the :class:`DynamicsSimulator`:

.. code-block:: cpp

      // Define list of dependent variables to save.
      std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariablesList;

where :literal:`dependentVariableList` is populated as follows:

.. code-block:: cpp

      dependentVariablesList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >( variableType , associatedBody , secondaryBody ) );

The details of the creation of the settings :class:`SingleDependentVariableSaveSettings` object are discussed below.

Available Dependent Variables
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Below, we provide a list of all dependent variables that can be saved using Tudat, with a link to the corresponding (derived class) of :class:`SingleDependentVariableSaveSettings`, and associated input. In some cases, the requirements on the environment for a variable to be saved are not necesarilly intuitive. In those cases, we mention the requirements explicitly. 

   - **Mach number** in atmosphere. Requires an aerodynamic acceleration to be acting on the vehicle. Defined by creating a :class:`SingleDependentVariableSaveSettings` object with input :literal:`mach_number_dependent_variable` as :literal:`variableType`.

   - **Altitude** above body exerting aerodynamic acceleration. Requires an aerodynamic acceleration to be acting on the vehicle. Defined by creating a :class:`SingleDependentVariableSaveSettings` object with input :literal:`altitude_dependent_variable` as :literal:`variableType`.

   - **Airspeed** in atmosphere of body exerting aerodynamic acceleration. Requires an aerodynamic acceleration to be acting on the vehicle. Defined by creating a :class:`SingleDependentVariableSaveSettings` object with input :literal:`airspeed_dependent_variable` as :literal:`variableType`.

   - **Local density** in atmosphere of body exerting aerodynamic acceleration (at position of body undergoing acceleration). Requires an aerodynamic acceleration to be acting on the vehicle. Defined by creating a :class:`SingleDependentVariableSaveSettings` object with input :literal:`local_density_dependent_variable` as :literal:`variableType`.

   - **Local temperature** in atmosphere of body exerting aerodynamic acceleration (at position of body undergoing acceleration). Requires an aerodynamic acceleration to be acting on the vehicle. Defined by creating a :class:`SingleDependentVariableSaveSettings` object with input :literal:`local_temperature_dependent_variable` as :literal:`variableType`.

   - **Local dynamic pressure** felt by the vehicle based on the current velocity and atmospheric conditions, i.e., computed with :math:`\frac{1}{2} \rho V^2`. Requires an aerodynamic acceleration to be acting on the vehicle. Defined by creating a :class:`SingleDependentVariableSaveSettings` object with input :literal:`local_dynamic_pressure_dependent_variable` as :literal:`variableType`.

   - **Local aerodynamic heat rate** felt by the vehicle based on the current velocity and atmospheric conditions, i.e., computed with :math:`\frac{1}{2} \rho V^3`. Requires an aerodynamic acceleration to be acting on the vehicle. Defined by creating a :class:`SingleDependentVariableSaveSettings` object with input :literal:`local_aerodynamic_heat_rate_dependent_variable` as :literal:`variableType`.

   - **Relative speed** (scalar velocity) of body w.r.t. a second body (between centers of mass). Defined by creating a :class:`SingleDependentVariableSaveSettings` object with input :literal:`relative_speed_dependent_variable` as :literal:`variableType`.

   - **Relative velocity** of body w.r.t. a second body (between centers of mass). Defined by creating a :class:`SingleDependentVariableSaveSettings` object with input :literal:`relative_velocity_dependent_variable` as :literal:`variableType`.

   - **Relative distance** of body from a second body (between centers of mass). Defined by creating a :class:`SingleDependentVariableSaveSettings` object with input :literal:`relative_distance_dependent_variable` as :literal:`variableType`.

   - **Relative position** of body w.r.t. a second body (between centers of mass). Defined by creating a :class:`SingleDependentVariableSaveSettings` object with input :literal:`relative_position_dependent_variable` as :literal:`variableType`.

   - **Radiation pressure coefficient** of body, due to radiation exerted by another body. Defined by creating a :class:`SingleDependentVariableSaveSettings` object with input :literal:`radiation_pressure_dependent_variable` as :literal:`variableType`.

   - **Total acceleration** acting on a body. Defined by creating a :class:`SingleDependentVariableSaveSettings` object with input :literal:`total_acceleration_dependent_variable` as :literal:`variableType`.

   - **Total torque** acting on a body. Defined by creating a :class:`SingleDependentVariableSaveSettings` object with input :literal:`total_torque_dependent_variable` as :literal:`variableType`.        
     
   - **Total mass rate** of body. Requires mass to be one of the numerically propagated variables. Defined by creating a :class:`SingleDependentVariableSaveSettings` object with input :literal:`total_mass_rate_dependent_variables` as :literal:`variableType`.

   - **Norm of total acceleration** acting on a body. Defined by creating a :class:`SingleDependentVariableSaveSettings` object with input :literal:`total_acceleration_norm_dependent_variable` as :literal:`variableType`.

   - **Norm of total torque** acting on a body. Defined by creating a :class:`SingleDependentVariableSaveSettings` object with input :literal:`total_torque_norm_dependent_variable` as :literal:`variableType`.

   - **Norm of single acceleration** acting on a body. Defined by creating a :class:`SingleAccelerationDependentVariableSaveSettings`, with :literal:`useNorm` set to true.

   - **Single acceleration** acting on a body. Defined by creating a :class:`SingleAccelerationDependentVariableSaveSettings`, with :literal:`useNorm` set to false.

   - **Spherical harmonic acceleration terms** acting on a body, thus the full list of accelerations *for each* spherical harmonics term. Defined by creating a :class:`SphericalHarmonicAccelerationTermsDependentVariableSaveSettings` object.

   - **Variation in spherical harmomic acceleration** acting on a body. Defined by creating a :class:`SingleVariationSphericalHarmonicAccelerationSaveSettings` object.

   - **Variation in spherical harmomic acceleration terms** acting on a body, thus the full list of accelerations *for each* spherical harmonics term.. Defined by creating a :class:`SingleVariationSingleTermSphericalHarmonicAccelerationSaveSettings` object.
    
   - **Norm of single torque** acting on a body. Defined by creating a :class:`SingleTorqueDependentVariableSaveSettings`, with :literal:`useNorm` set to true.

   - **Single torque** acting on a body. Defined by creating a :class:`SingleTorqueDependentVariableSaveSettings`, with :literal:`useNorm` set to false.

   - **Aerodynamic force coefficients** of a body. Defined by creating a :class:`SingleDependentVariableSaveSettings` object with input :literal:`aerodynamic_force_coefficients_dependent_variable` as :literal:`variableType`.

   - **Aerodynamic moment coefficients** of a body. Defined by creating a :class:`SingleDependentVariableSaveSettings` object with input :literal:`aerodynamic_moment_coefficients_dependent_variable` as :literal:`variableType`.

   - **Rotation matrix to body-fixed frame** of a body. Defined by creating a :class:`SingleDependentVariableSaveSettings` object with input :literal:`rotation_matrix_to_body_fixed_frame_variable` as :literal:`variableType`. 

   - **Rotation matrix between frames**  used for aerodynamics.  Defined by creating a :class:`IntermediateAerodynamicRotationVariableSaveSettings` class, with the two frames (start and end frames) provided as input. The following frames can be used (see Mooij, 1994 for details):

      - Inertial frame
      - Body-fixed (corotating) frame of central body. 
      - Vehicle-centered vertical frame
      - Vehicle-centered trajectory frame
      - Vehicle-centered aerodynamic frame
      - Vehicle-fixed body frame

   - **Rotation angle**  used for aerodynamics.  Defined by creating a :class:`BodyAerodynamicAngleVariableSaveSettings` class, with the desired angle provided as input. The following angles can be used (see Mooij, 1994 for details):

      - Latitude angle
      - Longitude angle
      - Heading angle
      - Flight-path angle
      - Angle of attack
      - Sideslip angle
      - Bank angle

   - **Airspeed-based velocity** vector (body velocity w.r.t. wind vector, assumes corotating atmosphere if no wind model is defined). Requires an aerodynamic acceleration to be acting on the vehicle.  Defined by creating a :class:`SingleDependentVariableSaveSettings` object with input :literal:`body_fixed_airspeed_based_velocity_variable` as :literal:`variableType`.

   - **Groundspeed-based velocity** vector (equal to airspeed-based velocity in absence of wind). Requires an aerodynamic acceleration to be acting on the vehicle.  Defined by creating a :class:`SingleDependentVariableSaveSettings` object with input :literal:`body_fixed_groundspeed_based_velocity_variable` as :literal:`variableType`.

   - **G-load** induced by aerodynamic acceleration. Requires an aerodynamic acceleration to be acting on the vehicle.  Defined by creating a :class:`SingleDependentVariableSaveSettings` object with input :literal:`total_aerodynamic_g_load_variable` as :literal:`variableType`.

   - **Stagnation point-heat flux** induced by atmospheric friction. Requires an aerodynamic acceleration to be acting on the vehicle, and a nose radius to be defined on the vehicle.  Defined by creating a :class:`SingleDependentVariableSaveSettings` object with input :literal:`stagnation_point_heat_flux_dependent_variable` as :literal:`variableType`.

   - **Geodetic latitude** (w.r.t. central body). Requires an aerodynamic acceleration to be acting on the vehicle.  Defined by creating a :class:`SingleDependentVariableSaveSettings` object with input :literal:`geodetic_latitude_dependent_variable` as :literal:`variableType`. 
     
   - **Control surface deflection** of a given aerodynamic control surface of body. Defined by creating a :class:`SingleDependentVariableSaveSettings` object with input :literal:`control_surface_deflection_dependent_variable` as :literal:`variableType`. 

   - **Keplerian state** of body. Defined by creating a :class:`SingleDependentVariableSaveSettings` object with input :literal:`keplerian_state_dependent_variable` as :literal:`variableType`. 

   - **Modified equinoctial state** of body. Defined by creating a :class:`SingleDependentVariableSaveSettings` object with input :literal:`modified_equinocial_state_dependent_variable` as :literal:`variableType`. The value of the parameter I is automatically chosen as +1 or -1, depending on whether the inclination is smaller or larger than 90 degrees.

   - **Relative Cartesian position**, based on position of vehicle and orbited body, expressed in the body fixed frame. Defined by creating a :class:`SingleDependentVariableSaveSettings` object with input :literal:`body_fixed_relative_cartesian_position` as :literal:`variableType`.

   - **Relative spherical position**, based on position of vehicle and orbited body, expressed in the body fixed frame. Defined by creating a :class:`SingleDependentVariableSaveSettings` object with input :literal:`body_fixed_relative_spherical_position` as :literal:`variableType`.

   - **Rotation of LVLH to inertial frame**, Rotation matrix from Local Vertical, Local Horizontal (LVLH) frame of body to inertial frame. Defined by creating a :class:`SingleDependentVariableSaveSettings` object with input :literal:`lvlh_to_inertial_frame_rotation_dependent_variable` as :literal:`variableType`. 

   - **Periapsis altitude**, based on current osculating elements. Defined by creating a :class:`SingleDependentVariableSaveSettings` object with input :literal:`periapsis_altitude_dependent_variable` as :literal:`variableType`.
    
      .. warning:: The computaton of the periapsis altitude uses the average radius of the central body, not the local radius.

Setting Up Dependent Variables
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The framework discussed in the previous section explains how the :literal:`dependentVariablesList` is populated and passed to the :class:`PropagatorSettings`. The goal of this section is to list the available dependent variables and to explain how these are pushed to the :literal:`dependentVariablesList`.

.. class:: SingleDependentVariableSaveSettings

   This base-class is a generic method to retrieve a large number of dependent variables that are not classified under a particular group. Variables are saved to the :literal:`dependentVariablesList` using the following code:

   .. code-block:: cpp

      dependentVariablesList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >( variableType , associatedBody , secondaryBody, componentIndex ) );

   where:

   - :literal:`variableType`

      :class:`PropagationDependentVariables` variable that can take the following values:

         - Variables returning **dependent variable of size 1**:

            - :literal:`mach_number_dependent_variable`
            - :literal:`altitude_dependent_variable`
            - :literal:`airspeed_dependent_variable`
            - :literal:`local_density_dependent_variable`
            - :literal:`relative_speed_dependent_variable` (secondary body defines body w.r.t. which the relative speed is computed)
            - :literal:`relative_distance_dependent_variable` (secondary body defines body w.r.t. which the relative distance is computed)
            - :literal:`radiation_pressure_dependent_variable` (secondary body defines the source of radiation for which the readiation pressure coefficient is to be provided)
            - :literal:`total_aerodynamic_g_load_variable` (secondary body defines body with atmosphere that exerts the aerodynamic acceleration that induces the g-load)
            - :literal:`stagnation_point_heat_flux_dependent_variable`
            - :literal:`local_temperature_dependent_variable`
            - :literal:`local_dynamic_pressure_dependent_variable`
            - :literal:`local_aerodynamic_heat_rate_dependent_variable`
            - :literal:`geodetic_latitude_dependent_variable`
            - :literal:`control_surface_deflection_dependent_variable` (secondary body defines name of control surface for which deflection is to be provided)
            - :literal:`total_mass_rate_dependent_variables`
            - :literal:`periapsis_altitude_dependent_variable` (secondary body defines body w.r.t. which the periapsis altitude is computed)
            - :literal:`total_torque_norm_dependent_variable`

         - Variables returning a **multi-valued dependent variable**:

            - :literal:`relative_position_dependent_variable` (secondary body defines body w.r.t. which the relative position is computed)
            - :literal:`relative_velocity_dependent_variable` (secondary body defines body w.r.t. which the relative velocity is computed)
            - :literal:`body_fixed_airspeed_based_velocity_variable`
            - :literal:`total_acceleration_norm_dependent_variable`
            - :literal:`total_acceleration_dependent_variable`
            - :literal:`aerodynamic_force_coefficients_dependent_variable`
            - :literal:`aerodynamic_moment_coefficients_dependent_variable`
            - :literal:`lvlh_to_inertial_frame_rotation_dependent_variable` (secondary body defines body w.r.t. which the state is computed when determining the matrix, taken as SSB if left empty)
            - :literal:`rotation_matrix_to_body_fixed_frame_variable`
            - :literal:`total_torque_dependent_variable`
            - :literal:`body_fixed_groundspeed_based_velocity_variable`
            - :literal:`keplerian_state_dependent_variable` (secondary body defines body w.r.t. which the Keplerian state is computed)
            - :literal:`modified_equinocial_state_dependent_variable` (secondary body defines body w.r.t. which the modified equinoctial state is computed)
            - :literal:`body_fixed_relative_cartesian_position` (secondary body defines body w.r.t. which the Cartesian state is computed)
            - :literal:`body_fixed_relative_spherical_position` (secondary body defines body w.r.t. which the spherical state is computed)

   - :literal:`associatedBody`

      Indicates to which body the saved dependent variables are associated.

   - :literal:`secondaryBody`

      Optional argument that provides a secondary body that may be necessary to save the dependent variable. By default, this argument is empty. In the list above, it is indicated which parameters require a secondaryBody to be defined, and what this parameter represents.

   - :literal:`componentIndex`

      :literal:`int` denoting which element of the acceleration vector is to be saved. Default value is -1, which means all elements are saved. Make sure that this value is compatible with the size of :literal:`variableType`.
      
.. class:: SingleAccelerationDependentVariableSaveSettings

   This derived class is used to retrieve acceleration-related dependent variables. A large number of acceleration models are supported and both the acceleration-norm and the acceleration-vector can be saved. Variables are added to the :literal:`dependentVariablesList` using the following code:

   .. code-block:: cpp

            dependentVariablesList.push_back(
                std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                accelerationModelType, bodyUndergoingAcceleration, bodyExertingAcceleration, useNorm, componentIndex );

   where:

   - :literal:`accelerationModelType`
  
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

      :literal:`std::string` variable that indicates the body that experiences the acceleration that needs to be retrieved. Make sure that the body's name is listed in :class:`NamedBodyMap`.

   - :literal:`bodyExertingAcceleration`

      :literal:`std::string` variable that indicates the body that exerts the acceleration that needs to be retrieved on :literal:`bodyUndergoingAcceleration`. Make sure that the body's name is listed in :class:`NamedBodyMap`.

   - :literal:`useNorm`

      :literal:`bool` variable that indicates if the norm of the acceleration (true) or the acceleration vector (false) must be retrieved.

   - :literal:`componentIndex`

      :literal:`int` denoting which element of the acceleration vector is to be saved. Default value is -1, which means all elements are saved.

   .. warning:: Make sure that the selected :literal:`bodyExertingAcceleration` is compatible with the :literal:`accelerationModelType`.

.. class:: SphericalHarmonicAccelerationTermsDependentVariableSaveSettings

   This derived class is used to save contributions to spherical harmonic acceleration at separate degree/order. This type of variable is added to the :literal:`dependentVariablesList` using the following code:

   .. code-block:: cpp

            dependentVariablesList.push_back(
                std::make_shared< SphericalHarmonicAccelerationTermsDependentVariableSaveSettings >(
                bodyUndergoingAcceleration, bodyExertingAcceleration, maximumDegree, maximumOrder, componentIndex );

   - :literal:`bodyUndergoingAcceleration`

      :literal:`std::string` variable that indicates the body that experiences the gravitational acceleration. Make sure that the body's name is listed in :class:`NamedBodyMap`.

   - :literal:`bodyExertingAcceleration`

      :literal:`std::string` variable that indicates the body that exerts the gravitational acceleration that needs to be retrieved on :literal:`bodyUndergoingAcceleration`. Make sure that the body's name is listed in :class:`NamedBodyMap`.

   - :literal:`maximumDegree`

      :literal:`int` denoting the maximum degree of the spherical harmonics model, for which an acceleration has to be saved.

   - :literal:`maximumOrder`

      :literal:`int` denoting the maximum order of the spherical harmonics model, for which an acceleration has to be saved. Make sure that this value is smaller or equal to :literal:`maximumDegree`.

   - :literal:`componentIndex`

      :literal:`int` denoting which element of the acceleration vector is to be saved. Default value is -1, which means all elements are saved.

.. class:: SingleVariationSphericalHarmonicAccelerationSaveSettings

   This derived class is used to save contributions to the variation of the spherical harmonic acceleration. This type of variable is added to the :literal:`dependentVariablesList` using the following code:

   .. code-block:: cpp

            dependentVariablesList.push_back(
                std::make_shared< SingleVariationSphericalHarmonicAccelerationSaveSettings >(
                bodyUndergoingAcceleration, bodyExertingAcceleration, deformationType, identifier );

   - :literal:`bodyUndergoingAcceleration`

      :literal:`std::string` variable that indicates the body that experiences the gravitational acceleration. Make sure that the body's name is listed in :class:`NamedBodyMap`.

   - :literal:`bodyExertingAcceleration`

      :literal:`std::string` variable that indicates the body that exerts the gravitational acceleration that needs to be retrieved on :literal:`bodyUndergoingAcceleration`. Make sure that the body's name is listed in :class:`NamedBodyMap`.

   - :literal:`deformationType`

      :class:`BodyDeformationTypes` variable defining the type of gravity field variation. The supported values are:

         - :literal:`basic_solid_body`
         - :literal:`tabulated_variation`

   - :literal:`identifier`

      :literal:`std::string` variable denoting the identifier for gravity field variation.

.. class:: SingleVariationSingleTermSphericalHarmonicAccelerationSaveSettings

   This derived class is used to save contributions to the variation of the spherical harmonic acceleration at separate degree/order. This type of variable is added to the :literal:`dependentVariablesList` using the following code:

   .. code-block:: cpp

            dependentVariablesList.push_back(
                std::make_shared< SingleVariationSingleTermSphericalHarmonicAccelerationSaveSettings >(
                bodyUndergoingAcceleration, bodyExertingAcceleration, maximumDegree, maximumOrder, deformationType, identifier );

   - :literal:`bodyUndergoingAcceleration`

      :literal:`std::string` variable that indicates the body that experiences the gravitational acceleration. Make sure that the body's name is listed in :class:`NamedBodyMap`.

   - :literal:`bodyExertingAcceleration`

      :literal:`std::string` variable that indicates the body that exerts the gravitational acceleration that needs to be retrieved on :literal:`bodyUndergoingAcceleration`. Make sure that the body's name is listed in :class:`NamedBodyMap`.

   - :literal:`maximumDegree`

      :literal:`int` denoting the maximum degree of the spherical harmonics model, for which an acceleration has to be saved.

   - :literal:`maximumOrder`

      :literal:`int` denoting the maximum order of the spherical harmonics model, for which an acceleration has to be saved. Make sure that this value is smaller or equal to :literal:`maximumDegree`.

   - :literal:`deformationType`

      :class:`BodyDeformationTypes` variable defining the type of gravity field variation. The supported values are:

         - :literal:`basic_solid_body`
         - :literal:`tabulated_variation`

   - :literal:`identifier`

      :literal:`std::string` variable denoting the identifier for gravity field variation.

.. class:: SingleTorqueDependentVariableSaveSettings

   This derived class is used to retrieve torque-related dependent variables. A number of torque models are supported and both the norm and the vector form of the torque can be saved. Variables are added to the :literal:`dependentVariablesList` using the following code:

   .. code-block:: cpp

            dependentVariablesList.push_back(
                std::make_shared< SingleTorqueDependentVariableSaveSettings >(
                torqueModelType, bodyUndergoingTorque, bodyExertingTorque, useNorm, componentIndex );

   where:

   - :literal:`torqueModelType`
  
      :class:`AvailableTorque` variable that defines the type of torque that must be retrieved. It can take the following values:
      
         - :literal:`underfined_torque`
         - :literal:`second_order_gravitational_torque`
         - :literal:`aerodynamic_torque`

   - :literal:`bodyUndergoingAcceleration`

      :literal:`std::string` variable that indicates the body that experiences the torque that needs to be retrieved. Make sure that the body's name is listed in :class:`NamedBodyMap`.

   - :literal:`bodyExertingAcceleration`

      :literal:`std::string` variable that indicates the body that exerts the torque that needs to be retrieved on :literal:`bodyUndergoingAcceleration`. Make sure that the body's name is listed in :class:`NamedBodyMap`.

   - :literal:`useNorm`

      :literal:`bool` variable that indicates if the norm of the torque (true) or the torque vector (false) must be retrieved.

   - :literal:`componentIndex`

      :literal:`int` denoting which element of the torque vector is to be saved. Default value is -1, which means all elements are saved.

   .. warning:: Make sure that the selected :literal:`bodyExertingAcceleration` is compatible with the :literal:`accelerationModelType`.

.. class:: IntermediateAerodynamicRotationVariableSaveSettings

   This derived class is used to retrieve the rotation matrix between two desired frames. Variables are saved to the :literal:`dependentVariablesList` using the following code:

   .. code-block:: cpp
      
            dependentVariablesList.push_back(
                std::make_shared< IntermediateAerodynamicRotationVariableSaveSettings >(
                    associatedBody, baseFrame, targetFrame )

   where:

   - :literal:`associatedBody`

      :literal:`std::string` variable that indicates the body for which a rotation matrix is to be saved. Make sure that the body's name is listed in :class:`NamedBodyMap`.

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

   This derived class is used to retrieve a number of rotation angles. Variables are saved to the :literal:`dependentVariablesList` using the following code:

   .. code-block:: cpp

      
            dependentVariablesList.push_back(
                std::make_shared< BodyAerodynamicAngleVariableSaveSettings >(
                    associatedBody, angle )

   where:

   - :literal:`associatedBody`

      :literal:`std::string` variable that indicates the body for which the :literal:`angle` is to be saved. Make sure that the body's name is listed in :class:`NamedBodyMap`.

   - :literal:`angle`

      :class:`AerodynamicsReferenceFrameAngles` variable that provides the angle to be saved. The following angles can be saved using this method:

      - :literal:`latitude_angle`
      - :literal:`longitude_angle`
      - :literal:`heading_angle`
      - :literal:`flight_path_angle`
      - :literal:`angle_of_attack`
      - :literal:`angle_of_sideslip`
      - :literal:`bank_angle`

