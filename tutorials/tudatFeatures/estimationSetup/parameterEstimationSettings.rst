.. _parameterEstimationSettings:

Setting Up Estimated Parameters
===============================

.. _parameterArchitecture:

Parameter Architecture
~~~~~~~~~~~~~~~~~~~~~~

The parameter estimation framework of Tudat allows an ever increasing variety of parameters to be estimated, these parameters may be:

* Properties of a body, such as a gravitational parameter :math:`\mu`
* Properties of a ground station, such as its body-fixed position :math:`\mathbf{x}_{GS}^{(B)}`
* Global properties of the simulation, such a Parameterize Post_Newtonian (PPN) parameters :math:`\gamma` and :math:`\beta`
* Acceleration model properties, such as empirical acceleration magnitudes
* Observation model properties, such as absolute and relative observation biases

In Tudat, these parameters influence the simulation in a variety of manners, and during propagation and/or observation simulation, information of this parameter is transferred in manner different ways. To provide a unified framework for estimating any type of parameter, the :class:`EstimatableParameter` class has been set up. 

.. class:: EstimatableParameter

   This class has interfaces to retrieve and reset parameters, providing a single interface for modifying/obtaining any of the parameters that Tudat supports. For each estimated parameter, there is a dedicated derived class of :class:`EstimatableParameter`. 

.. class:: EstimatableParameterSet

   The full list of estimated parameters is stored in an object of type :class:`EstimatableParameterSet`. This class is templated by the state scalar type of the estimated initial state parameters. 

.. note::
   For the remainder of this page, we will implicitly assume that the template argument of an :class:`EstimatableParameterSet` object is double, unless explicitly mentioned otherwise.

As is the case for acceleration models, integration models, environment models, *etc.*, the parameter objects are created by defining settings for them, and subsequently calling the associated factory function. The settings and passed by creating objects of type :class:`EstimatableParameterSettings` (or one of its derived classes). Some parameters settings are provided through the :class:`EstimatableParameterSettings` base class, and some through its derived classes. A full list is provided in the section on :ref:`parameterSettingCreation`. An example of the creation of the parameter objects is given below:

   .. code-block:: cpp

       // Define parameter settings
       std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;
       parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( 
          "Vehicle", radiation_pressure_coefficient ) );
       parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( 
          "Vehicle", constant_drag_coefficient ) );
       parameterNames.push_back(  std::make_shared< EstimatableParameterSettings >(
          "Earth", rotation_pole_position ) );
          
       // Define parameter objects   
       std::shared_ptr< EstimatableParameterSet< double > > parametersToEstimate =
            createParametersToEstimate( parameterNames, bodyMap );
         
Which creates parameter objects for the radiation pressure coefficient and drag coefficient of body "Vehicle", and the orientation of the rotation axis of the body "Earth".

The :class:`EstimatableParameterSet` object contains three objects that have :class:`EstimatableParameter` as base class (one for each parameter). We distinguish two types of :class:`EstimatableParameter` objects:

* Those that represent initial conditions for dynamics (denoted as :math:`\mathbf{x}_{0}` below)
* Those that represent fixed parameters for environment, acceleration or observation models (denoted as :math:`\mathbf{q}` below)

Resetting the full parameter vector :math:`\mathbf{p}(=[\mathbf{x}_{0};\mathbf{q}])` is done as follows (for :literal:`double` state scalar type):
         
   .. code-block:: cpp

       // Create parameter set  
       std::shared_ptr< EstimatableParameterSet< double > > parametersToEstimate = ...
       
       Eigen::VectorXd parameterVector =
            parametersToEstimate->getFullParameterValues< double >( );

While resetting the full parameter vector is done as:

   .. code-block:: cpp

       // Create parameter set  
       std::shared_ptr< EstimatableParameterSet< double > > parametersToEstimate = ...
       
       // Define vector of new values of estimated parameters
       Eigen::VectorXd newParameterVector = ...
       
       // Reset parameter values
       parametersToEstimate->resetParameterValues< double >( );

When resetting the parameter vector, the change in the values in :math:`\mathbf{q}` immediately take effect. For the initial state parameters to take effect, however, the dynamics must be re-propagated. This occurs automatically when estimating parameters. It can also be performed manually by calling the :literal:`resetParameterEstimate` member function of the :class:`VariationalEquationsSolver` class. 

.. _parameterSettingCreation:

Creating Estimated Parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The framework discussed in the previous section explains how the :literal:`parameterNames` is populated. The goal of this section is to list the available parameters that can be estimated, and which environment models they are linked to.


.. class:: EstimatableParameterSettings

   This base-class is a generic method to define parameters that require no more information than their type, the associated body, and (in some cases) a secondary identifier. Variables are added to the :literal:`parameterNames` using the following code:

   .. code-block:: cpp

      parameterNames.push_back(
                   std::make_shared< EstimatableParameterSettings >( associatedBody, parameterType, secondaryIdentifier ) );

   where:
   - :literal:`associatedBody`
   
      Name of body for which the parameter is estimated, as :literal:`std::string`.
   
   - :literal:`parameterType`

      :class:`EstimatebleParametersEnum` variable that can take the following values:
      
      - :literal:`gravitational_parameter`. Gravitational parameter of a body, linked to a :class:`GravityFieldModel` object, which may be a point-mass or (time-dependent) spherical harmonic field. Parameter size: 1. Secondary identifer: None.
      - :literal:`constant_drag_coefficient`. Drag coefficient of a body that is constant, linked to a :class:`CustomAerodynamicCoefficientInterface` object derived from :class:`AerodynamicCoefficientInterface`, which must have 0 independent variables for the coefficients. Parameter size: 1. Secondary identifer: None.
      - :literal:`constant_rotation_rate`. Rotation rate of a body around a fixed axis, linked to a :class:`SimpleRotationalEphemeris` object derived from :class:`RotationalEphemeris`. Parameter size: 1. Secondary identifer: None.
      - :literal:`radiation_pressure_coefficient`. Constant radiation pressure coefficient of a body, linked to a :class:`RadiationPressureInterface` object. Parameter size: 1. Secondary identifer: None.
      - :literal:`rotation_pole_position`. Fixed rotation axis about which a body rotates with a fixed rotation rate, linked to a :class:`SimpleRotationalEphemeris` object. Parameter size: 2 (denoting pole right ascension and declination). Secondary identifer: None.
      - :literal:`ground_station_position`. Fixed body-fixed position of a ground station on a body, linked to a :class:`GroundStationState` object (requires a :class:`GroundStationState` class). Parameter size: 3 (denoting body-fixed *x*, *y* and *z* Cartesian position). Secondary identifer: Ground station name.
      - :literal:`ppn_parameter_gamma`. Parameter :math:`\gamma` used in Parametric Post-Newtonian (PPN) framework, linked to a :class:`PPNParameterSet` object (nominally the global :literal:`relativity::ppnParameterSet` variable). Parameter size: 1. Note that the name of the associated body should be :literal:`"global_metric"`. Secondary identifer: None.
      - :literal:`ppn_parameter_beta`. Parameter :math:`\beta` used in Parametric Post-Newtonian (PPN) framework, linked to a :class:`PPNParameterSet` object (nominally the global :literal:`relativity::ppnParameterSet` variable). Parameter size: 1. Note that the name of the associated body should be :literal:`"global_metric"`. Secondary identifer: None.
      - :literal:`equivalence_principle_lpi_violation_parameter`. Parameter used to compute influence of a gravitational potential on proper time rate, equals 0 in general relativity, not linked to any object, but instead the :literal:`equivalencePrincipleLpiViolationParameter` global variable (in namespace :literal:`relativity`. Parameter size: 1. Note that the name of the associated body should be :literal:`"global_metric"`. Secondary identifer: None.

   - :literal:`secondaryIdentifier`
   
      Secondary identifier to define the estimated parameter (if necessary, see above), as :literal:`std::string`. Empty by default. 

.. class:: InitialTranslationalStateEstimatableParameterSettings
   
   This derived class of :class:`EstimatableParameterSettings` is used to define settings for estimating an initial translational state of a single body. It is templated by the state scalar type of the dynamics (typically :literal:`double`). Two constructors are available for this class. The first constructor explicitly provides the current value of the initial state, while the second extracts the initial state from the current ephemeris of the body. The first constructor is called as:
   
   .. code-block:: cpp
   
            parameterNames.push_back(
                   std::make_shared< InitialTranslationalStateEstimatableParameterSettings< double > >( 
                      associatedBody, initialStateValue, centralBody, frameOrientation  ) );
                      
                      
   where:
   
   - :literal:`associatedBody`
   
      Name of body for which the initial state is to be estimated, as :literal:`std::string`.
      
   - :literal:`initialStateValue`
   
      Initial state (:literal:`Eigen::Vector6d`) of :literal:`associatedBody` w.r.t.  :literal:`centralBody` in Cartesian coordinates, expressed frame :literal:`frameOrientation`.
         
   - :literal:`centralBody`
   
      Name of body w.r.t. which the dynamics is to be estimated, as :literal:`std::string`.
         
   - :literal:`frameOrientation`
   
      Orientation in the frame in which the dynamics is to be estimated, as :literal:`std::string`.
      
      
The second constructor is called as:
   
   .. code-block:: cpp
   
            parameterNames.push_back(
                   std::make_shared< InitialTranslationalStateEstimatableParameterSettings< double > >( 
                      associatedBody, initialTime, centralBody, frameOrientation  ) );
                      
                      
   where:
   
   - :literal:`associatedBody`
   
      Name of body for which the initial state is to be estimated, as :literal:`std::string`.
      
   - :literal:`initialTime`
   
      Time (:literal:`double`) at which the body's ephemeris is to be interrogated to extract the initial state. If the ephemeris origin is not equal to :literal:`centralBody`, any required frame translations are applied.
         
   - :literal:`centralBody`
   
      Name of body w.r.t. which the dynamics is to be estimated, as :literal:`std::string`.
         
   - :literal:`frameOrientation`
   
      Orientation in the frame in which the dynamics is to be estimated, as :literal:`std::string`.      

.. class:: ArcWiseInitialTranslationalStateEstimatableParameterSettings
   
   This derived class of :class:`EstimatableParameterSettings` is used to define settings for estimating an initial translational state of a single body in an arc-wise manner. It is templated by the state scalar type of the dynamics (typically :literal:`double`). As was the case for the :class:`InitialTranslationalStateEstimatableParameterSettings` class, two constructors are available for this class. The first constructor explicitly provides the current value of the arc initial states, while the second extracts the arc initial states from the current ephemeris of the body. The first constructor is called as:
   
   .. code-block:: cpp
   
            parameterNames.push_back(
                   std::make_shared< ArcWiseInitialTranslationalStateEstimatableParameterSettings< double > >( 
                      associatedBody, concatenatedInitialStateValues, arcStartTimes, centralBody, frameOrientation  ) );
                      
                      
   where:
   
   - :literal:`associatedBody`
   
      Name of body for which the initial state is to be estimated, as :literal:`std::string`.
      
   - :literal:`concatenatedInitialStateValues`
   
      Initial states (as :literal:`Eigen::VectorXd`) of :literal:`associatedBody` w.r.t.  :literal:`centralBody` in Cartesian coordinates, expressed frame :literal:`frameOrientation`. This vector consists of the concatenated arc initial states, and must have a size :literal:`6 * arcStartTimes.size( )`. With the initial state of arc :math:`j` denotes as :math:`\mathbf{x}_{j,0}`, this input vector must be given as :math:`[\mathbf{x}_{1,0};\mathbf{x}_{1,0};...;\mathbf{x}_{N,0}]`, for :math:`N` arcs.

   - :literal:`arcStartTimes`
         
      List of times (:literal:`std::vector< double >`) at which the arcs for which the dynamics is to be estimated start. The entries of this arcs must be continuously increasing: denoting the start time of arc :math:`j` as :math:`t_{j,0}`, :math:`t_{j+1,0}>t_{j,0}` must always be satisfied.
      
   - :literal:`centralBody`
   
      Name of body w.r.t. which the dynamics is to be estimated, as :literal:`std::string`.
         
   - :literal:`frameOrientation`
   
      Orientation in the frame in which the dynamics is to be estimated, as :literal:`std::string`.
      
      
The second constructor is called as:
   
   .. code-block:: cpp
   
            parameterNames.push_back(
                   std::make_shared< ArcWiseInitialTranslationalStateEstimatableParameterSettings< double > >( 
                      associatedBody, arcStartTimes, centralBody, frameOrientation  ) );
                      
                      
   where:
   
   - :literal:`associatedBody`
   
      Name of body for which the initial state is to be estimated, as :literal:`std::string`.

   - :literal:`arcStartTimes`
         
      List of times (:literal:`std::vector< double >`) at which the arcs for which the dynamics is to be estimated start. The entries of this arcs must be continuously increasing: denoting the start time of arc :math:`j` as :math:`t_{j,0}`, :math:`t_{j+1,0}>t_{j,0}` must always be satisfied. The role of this variable is two-fold: firstly to inform the propagation on the times at which to switch to a new arc, and secondly to retrieve the initial state of each arc from the body's ephemeris.
      
   - :literal:`centralBody`
   
      Name of body w.r.t. which the dynamics is to be estimated, as :literal:`std::string`.
         
   - :literal:`frameOrientation`
   
      Orientation in the frame in which the dynamics is to be estimated, as :literal:`std::string`. 
         
.. class:: SphericalHarmonicEstimatableParameterSettings
   
   This derived class of :class:`EstimatableParameterSettings` is used to define settings for estimating spherical harmonic coefficients. Two constructors are available for this class. The first constructor allows a specific set of coefficients to be estimated and is used as follows:
   
   .. code-block:: cpp
   
            parameterNames.push_back(
                   std::make_shared< SphericalHarmonicEstimatableParameterSettings >( 
                      blockIndices, associatedBody, coefficientType ) );
                      
                      
   where:
   
   - :literal:`blockIndices`

      A list of degrees/orders at which the coefficients should be estimated, of type :literal:`std::vector< std::pair< int, int > >`, with a single pair entry denoting the degree (first) and order (second) of a coefficient that is to be estimated. The vector can have arbitrary size, but may not contain repeated coefficients.
      
   - :literal:`associatedBody`

      The name of the body for which coefficients are to be estimated (as :literal:`std::string`).
  
  - :literal:`coefficientType`
  
      The type of coefficients that are estimated, must be either :literal:`spherical_harmonics_cosine_coefficient_block` or :literal:`spherical_harmonics_sine_coefficient_block`. 
      
   The second constructor defines a 'full' block of coefficients to be estimated: it sets all degrees and orders from given minimum and maximum degrees and orders as the estimated parameters. It is created from:
     
    .. code-block:: cpp
   
         parameterNames.push_back(
                std::make_shared< SphericalHarmonicEstimatableParameterSettings >( 
                   minimumDegree, minimumOrder, maximumDegree, maximumOrder, associatedBody, coefficientType ) );
                      
                      
   where:
   
   - :literal:`minimumDegree`

      Minimum degree of coefficients that are estimated.
       
   - :literal:`minimumOrder`

      Minimum order of coefficients that are estimated.
         
   - :literal:`maximumDegree`

      Maximum degree of coefficients that are estimated.
       
   - :literal:`maximumOrder`

      Maximum order of coefficients that are estimated.
      
   - :literal:`associatedBody`

      The name of the body for which coefficients are to be estimated (as :literal:`std::string`).
  
  - :literal:`coefficientType`
  
      The type of coefficients that are estimated, must be either :literal:`spherical_harmonics_cosine_coefficient_block` or :literal:`spherical_harmonics_sine_coefficient_block`. 
      
   As an example:
   
   
   .. code-block:: cpp
   
       parameterNames.push_back(
          std::make_shared< SphericalHarmonicEstimatableParameterSettings >( 
             2, 1, 4, 3, "Earth", spherical_harmonics_cosine_coefficient_block ) );
            
            
   Sets :math:`C_{21}, C_{22}, C_{31}, C_{32},C_{33}, C_{41}, C_{42}, C_{43}` as the estimated parameters.
            
   .. note::
      
      If a full set of coefficients is to be estimated, ensure that two :class:`SphericalHarmonicEstimatableParameterSettings` objects are created: one for cosine and one for sine coefficients (unless only order 0 coefficients are estimated, which have no sine coefficients)
         
.. class:: FullDegreeTidalLoveNumberEstimatableParameterSettings

   This derived class of :class:`EstimatableParameterSettings` is used to define settings for tidal Love numbers :math:`k_{n}` at degree :math:`n` that are constant over all all orders at that degree (so for :math:`n=2, k_{20}=k_{21}=k_{22}`. The estimated parameters are a property of an :class:`BasicSolidBodyTideGravityFieldVariations` object. It is created by:
      
   .. code-block:: cpp
   
         parameterNames.push_back(
            std::make_shared< FullDegreeTidalLoveNumberEstimatableParameterSettings >( 
               associatedBody, degree, deformingBodies, useComplexValue ) );
                      
   where:
   
   - :literal:`associatedBody`

      An :literal:`std::string` that gives the name of the body for which the Love number is to be estimated.
      
   - :literal:`degree`
   
      An :literal:`int` that denotes the degree :math:`n` of the Love number that is estimated.
  
   - :literal:`deformingBodies` 
  
      List of bodies that cause tidal deformation of :literal:`associatedBody`, as an :literal:`std::vector< std::string >`. If, and only if, the body only has one :literal:`BasicSolidBodyTideGravityFieldVariations`, this list may be left empty and this single tidal model is used for estimating :math:`k_{n}`. If this list of deforming bodies is not empty, it must match *exactly* the list of deforming bodies of the tidal model. In this way, multiple Love numbers at different forcing frequencies can be estimated (by creating multiple :literal:`FullDegreeTidalLoveNumberEstimatableParameterSettings` objects.
  
   - :literal:`useComplexValue`
   
      A :literal:`bool` that denotes whether to estimate the Love number as a real value (size 1) or a complex value (size 2). In the complex case, the imaginary part represents the impact of tidal dissipation.
      
.. class:: SingleDegreeVariableTidalLoveNumberEstimatableParameterSettings

   This derived class of :class:`EstimatableParameterSettings` is used to define settings for tidal Love numbers :math:`k_{nm}` at degree :math:`n` that are vary over the orders at that degree. The estimated parameters are a property of an :class:`BasicSolidBodyTideGravityFieldVariations` object. It is created by:
      
   .. code-block:: cpp
   
         parameterNames.push_back(
            std::make_shared< SingleDegreeVariableTidalLoveNumberEstimatableParameterSettings >( 
               associatedBody, degree, orders deformingBodies, useComplexValue ) );
                      
   where:
   
   - :literal:`associatedBody`

      An :literal:`std::string` that gives the name of the body for which the Love numbers are to be estimated.
      
   - :literal:`degree`
   
      An :literal:`int` that denotes the degree :math:`n` of the Love numbers that are estimated.
   
   - :literal:`orders`
   
      An :literal:`std::vector< int >` that denotes the orders :math:`m` of the Love numbers that are to be estimated. For instance, for :literal:`degree = 3` and :literal:`orders = {2, 1, 3}`, :math:`k_{32}, k_{31}` and :math:`k_{33}` are estimated
      
   - :literal:`deformingBodies` 
  
      List of bodies that cause tidal deformation of :literal:`associatedBody`, as an :literal:`std::vector< std::string >`. If, and only if, the body only has one :literal:`BasicSolidBodyTideGravityFieldVariations`, this list may be left empty and this single tidal model is used for estimating :math:`k_{nm}`. If this list of deforming bodies is not empty, it must match *exactly* the list of deforming bodies of the tidal model. In this way, multiple Love numbers at different forcing frequencies can be estimated (by creating multiple :literal:`SingleDegreeVariableTidalLoveNumberEstimatableParameterSettings` objects.
  
   - :literal:`useComplexValue`
   
      A :literal:`bool` that denotes whether to estimate the Love numbers as real value (size 1 per Love number) or a complex value (size 2 per Love number). In the complex case, the imaginary parts represent the impact of tidal dissipation.

.. class:: ConstantObservationBiasEstimatableParameterSettings

   This derived class of :class:`EstimatableParameterSettings` is used to define settings for estimating a constant absolute or relative observation biases for a given set of :literal:`LinkEnds` and :literal:`ObservableType`. Depending on the input, an object of this class defines settings for estimating an abolute *or* a relative bias. The bias model itself is created by using the :class:`ConstantObservationBiasSettings` or :class:`ConstantRelativeObservationBiasSettings`, and the estimated parameter is a property of an :class:`ConstantObservationBias` of :class:`ConstantRelativeObservationBias` object. The parameter estimation settings are creating by:
      
   .. code-block:: cpp
   
         parameterNames.push_back(
            std::make_shared< ConstantObservationBiasEstimatableParameterSettings >( 
               linkEnds, observableType, isBiasAdditive ) );
                      
   where:
   
   - :literal:`linkEnds`

      A :literal:`LinkEnds` map that defines the link ends (receiver, transmitter, etc.) for the bias.
      
   - :literal:`observableType`

      An :literal:`ObservableType` variable that denotes the type of observable for which the bias is to be estimated.

   - :literal:`isBiasAdditive`

      A :literal:`bool` that is true if the bias is absolute and false if it is relative.
      
.. class:: EmpiricalAccelerationEstimatableParameterSettings

   This derived class of :class:`EstimatableParameterSettings` is used to define settings for estimating a set of empirical accelerations that are constant throughout the propation interval. Coefficients can be estimated for each of the three directions in the RSW frame, and as a constant term, as well a sine/cosine of true anomaly (for a total of 9 possible coefficients). The parameter is a property of an object of type :literal:`EmpiricalAcceleration`, which is created by using an object of type :class:`EmpiricalAccelerationSettings` The parameter estimation settings are creating by:
      
   .. code-block:: cpp
   
         parameterNames.push_back(
            std::make_shared< EmpiricalAccelerationEstimatableParameterSettings >( 
               associatedBody, centralBody, componentsToEstimate ) );
                      
   where:
   
   - :literal:`associatedBody`

      Name of body for which accelerations are to be estimated, as :literal:`std::string`.
      
   - :literal:`centralBody`

      Name of body central body about which the accelerated body is orbiting (e.g. the body w.r.t. which the Kepler elements are calculated), as :literal:`std::string`.

   - :literal:`componentsToEstimate`

      A list that defines the list of components that are to be estimated, of type :literal:`std::map< EmpiricalAccelerationComponents, std::vector< EmpiricalAccelerationFunctionalShapes > >`. The :literal:`EmpiricalAccelerationComponents` denotes the direction of the acceleration, and :literal:`EmpiricalAccelerationFunctionalShapes` whether a constant, sine or cosine coefficient is estimated. :literal:`EmpiricalAccelerationComponents` can take the values:
   
      * radial_empirical_acceleration_component
      * along_track_empirical_acceleration_component
      * across_track_empirical_acceleration_component
    
     :literal:`EmpiricalAccelerationFunctionalShapes` can take the values:

      * constant_empirical
      * sine_empirical
      * cosine_empirical

     For instance, to estimate a constant along-track, a sine and cosine radial, and a constant and sine across-track empirical acceleration, the :literal:`componentsToEstimate` becomes:
     
      .. code-block:: cpp
   
         std::map< EmpiricalAccelerationComponents, std::vector< EmpiricalAccelerationFunctionalShapes > > componentsToEstimate;
         componentsToEstimate[ along_track_empirical_acceleration_component ].push_back( constant_empirical );
         
         componentsToEstimate[ radial_empirical_acceleration_component ].push_back( sine_empirical );
         componentsToEstimate[ radial_empirical_acceleration_component ].push_back( cosine_empirical );

         componentsToEstimate[ across_track_empirical_acceleration_component ].push_back( constant_empirical );
         componentsToEstimate[ across_track_empirical_acceleration_component ].push_back( sine_empirical );
                      
.. class:: ArcWiseEmpiricalAccelerationEstimatableParameterSettings

   This derived class of :class:`EstimatableParameterSettings` is used to define settings for estimating a set of empirical acceleration components that are arc-wise constant. Coefficients can be estimated for each of the three directions in the RSW frame, and as a constant term, as well a sine/cosine of true anomaly (for a total of 9 possible coefficients). The parameter is a property of an object of type :literal:`EmpiricalAcceleration`, which is created by using an object of type :class:`EmpiricalAccelerationSettings` The parameter estimation settings are creating by:
      
   .. code-block:: cpp
   
         parameterNames.push_back(
            std::make_shared< EmpiricalAccelerationEstimatableParameterSettings >( 
               associatedBody, centralBody, componentsToEstimate, arcStartTimeList ) );
                      
   where:
   
   - :literal:`associatedBody`

      Name of body for which accelerations are to be estimated, as :literal:`std::string`.
      
   - :literal:`centralBody`

      Name of body central body about which the accelerated body is orbiting (e.g. the body w.r.t. which the Kepler elements are calculated), as :literal:`std::string`.

   - :literal:`componentsToEstimate`

      A list that defines the list of components that are to be estimated, of type :literal:`std::map< EmpiricalAccelerationComponents, std::vector< EmpiricalAccelerationFunctionalShapes > >`. The :literal:`EmpiricalAccelerationComponents` denotes the direction of the acceleration, and :literal:`EmpiricalAccelerationFunctionalShapes` whether a constant, sine or cosine coefficient is estimated. See above in documentation for :class:`EmpiricalAccelerationEstimatableParameterSettings` for more details.
   
   - :literal:`arcStartTimeList`
   
      A list of times at which the arcs start during which the coefficients are to be estimated, of :literal:`type std::vector< double >`. For instance, when using:
      
      .. code-block:: cpp
      
         std::vector< double > arcStartTimeList;
         arcStartTimeList.push_back( 1000.0 );
         arcStartTimeList.push_back( 7200.0 );
         arcStartTimeList.push_back( 10000.0 );

      One set of empirical accelerations will be estimated, that are used for :math:`1000<t<7200`, one set for :math:`7200<t<10000` and one set for :math:`t>10000`
      
   
