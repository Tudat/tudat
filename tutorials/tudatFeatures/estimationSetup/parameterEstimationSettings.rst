.. _parameterEstimationSettings:

Setting Up Estimated Parameters
===============================

Parameter Architecture
~~~~~~~~~~~~~~~~~~~~~~

Creating Estimated Parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The framework discussed in the previous section explains how the :literal:`parameterNames` is populated. The goal of this section is to list the available parameters that can be estimated, and which environment models they are linked to.


.. class:: EstimatableParameterSettings

   This base-class is a generic method to define parameters that require no more information than their type, the associated body, and (in some cases) a secondary identifier. Variables are added to the :literal:`parameterNames` using the following code:

   .. code-block:: cpp

      parameterNames.push_back(
                   boost::make_shared< EstimatableParameterSettings >( associatedBody, parameterType, secondaryIdentifier ) );

   where:
   
   - :literal:`parameterType`

      :class:`EstimatebleParametersEnum` variable that can take the following values:
      
      - :literal:`gravitational_parameter`. Gravitational parameter of a body, linked to a :class:`GravityFieldModel` object, which may be a point-mass or (time-dependent) spherical harmonic field. Parameter size: 1. Secondary identifer: None.
      - :literal:`constant_drag_coefficient`. Drag coefficient of a body that is constant, linked to a :class:`CustomAerodynamicCoefficientInterface` object, which must have 0 independent variables for the coefficients. Parameter size: 1. Secondary identifer: None.
      - :literal:`constant_rotation_rate`. Rotation rate of a body around a fixed axis, linked to a :class:`SimpleRotationalEphemeris` object. Parameter size: 1. Secondary identifer: None.
      - :literal:`radiation_pressure_coefficient`. Constant radiation pressure coefficient of a body, linked to a :class:`RadiationPressureInterface` object. Parameter size: 1. Secondary identifer: None.
      - :literal:`rotation_pole_position`. Fixed rotation axis about which a body rotates with a fixed rotation rate, linked to a :class:`` object. Parameter size: 2 (denoting pole right ascension and declination). Secondary identifer: None.
      - :literal:`ground_station_position`. Fixed body-fixed position of a ground station on a body, linked to a :class:`GroundStationState` object (requires a :class:`GroundStationState` class). Parameter size: 3 (denoting body-fixed *x*, *y* and *z* Cartesian position). Secondary identifer: Ground station name.
      - :literal:`ppn_parameter_gamma`. Parameter :math:`\gamma` used in Parametric Post-Newtonian (PPN) framework, linked to a :class:`PPNParameterSet` object (nominally the global :literal:`relativity::ppnParameterSet` variable). Parameter size: 1. Note that the name of the associated body should be :literal:`"global_metric"`. Secondary identifer: None.
      - :literal:`ppn_parameter_gamma`. Parameter :math:`\beta` used in Parametric Post-Newtonian (PPN) framework, linked to a :class:`PPNParameterSet` object (nominally the global :literal:`relativity::ppnParameterSet` variable). Parameter size: 1. Note that the name of the associated body should be :literal:`"global_metric"`. Secondary identifer: None.

    arc_wise_initial_body_state,
    initial_body_state,
    spherical_harmonics_cosine_coefficient_block,
    spherical_harmonics_sine_coefficient_block,
    constant_additive_observation_bias,
    constant_relative_observation_bias,
    equivalence_principle_lpi_violation_parameter,
    empirical_acceleration_coefficients,
    arc_wise_empirical_acceleration_coefficients,
    full_degree_tidal_love_number,
    single_degree_variable_tidal_love_number
    
Retrieving and modifying parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    