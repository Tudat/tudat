.. _tudatFeaturesPropagatorSettings:

Propagator Settings: Basics
===========================
This page presents an overview of the available options within :class:`PropagatorSettings`. As the name suggests, these settings define how the orbit is propagated within the :class:`DynamicsSimulator`.

Similarly to the :class:`IntegratorSettings` discussed in :ref:`tudatFeaturesIntegratorSettings`, various derived classes are used to implement different settings:

.. class:: PropagatorSettings

   Base class from which the other settings classes described below are derived.

.. class:: TranslationalStatePropagatorSettings

    This class implements the framework required to propagate the translation state of a body. The constructor of this derived class is overloaded allowing two types of termination conditions:

    - Termination when the simulation time reaches a predefined :literal:`endTime` (Default).
    - Termination when a predefined dependent variables meets a certain criterion.

    .. method:: Default termination settings

        .. code-block:: cpp

            TranslationalStatePropagatorSettings<StateScalarType>( centralBodies,
                                                                   accelerationsMap,
                                                                   bodiesToIntegrate,
                                                                   initialBodyStates,
                                                                   endTime,
                                                                   propagator,
                                                                   dependentVariablesToSave )

        where:

        - :literal:`StateScalarType`
   
            Template argument used to set the precision of the state, in general :literal:`double` is used. For some application where a high precision is required this can be changed to e.g. :literal`long double`. 
        
        - :literal:`centralBodies`

            :literal:`std::vector< std::string >` that contains the names of the central bodies and must match with those in the :class:`BodyMap`.

        - :literal:`accelerationsMap`

            :class:`AccelerationMap` that contains the accelerations for each body as discussed in :ref:`tudatFeaturesAccelerationIndex`.

        - :literal:`bodiesToIntegrate`

            :literal:`std::vector< std::string >` that contains the names of the bodies to integrate which must match with those in the :class:`BodyMap`.

        - :literal:`initialBodyStates`

            :literal:`Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >` that stores the states of the bodies to propagate with respect to their central bodies. 

        - :literal:`endTime`

            :literal:`double` that defines the end-time of the simulation.

        - :literal:`propagator`

            :class:`TranslationalPropagatorType` which defines the type of propagator to be used. Currently, the following propagators are supported: 

               - :literal:`cowell`
               - :literal:`encke`
               - :literal:`gauss_keplerian`
               - :literal:`gauss_modified_equinoctial`
               - :literal:`unified_state_model_quaternions`
               - :literal:`unified_state_model_modified_rodrigues_parameters`
               - :literal:`unified_state_model_exponential_map`

            By default, the :literal:`cowell` propagator is used.

        - :literal:`dependentVariablesToSave`

            :literal:`boost::shared_ptr< DependentVariableSaveSettings >` that presents a list of the dependent variables to save during propagation. How this is exactly done is explained below. By default, an empty list is used and no dependent variable is saved. See the tutorial on :class:`DependentVariableSaveSettings` for more details on this class.

        .. note:: The state variables contained in :literal:`initialBodyStates` are ordered with respect to the elements of :literal:`centralBodies` and :literal:`bodiesToIntegrate`. Please take a look at the following pseudocode:

            .. code-block:: cpp

                centralBodies = { Sun , Earth , Moon }
                bodiesToIntegrate = { Earth , Moon }
                initialBodyStates = { xEarthWrtSun , yEarthWrtSun , zEarthWrtSun , uEarthWrtSun , vEarthWrtSun , wEarthWrtSun , 
                                      xMoonWrtEarth , yMoonWrtEarth , zMoonWrtEarth , uMoonWrtEarth , vMoonWrtEarth , wMoonWrtEarth }
            

    .. method:: User-defined termination settings

        .. code-block:: cpp

            TranslationalStatePropagatorSettings<StateScalarType>( centralBodies,
                                                                   accelerationsMap,
                                                                   bodiesToIntegrate,
                                                                   initialBodyStates,
                                                                   terminationSettings,
                                                                   propagator,
                                                                   dependentVariablesToSave )

        where:

        - :literal:`terminationSettings`

            :literal:`boost::shared_ptr< PropagationTerminationSettings >` that defines the termination settings of the propagation. This is the fifth argument and replaces the :literal:`endTime` in the default constructor. See the tutorial on :class:`PropagationTerminationSettings` for more details on this class.

.. class:: RotationalStatePropagatorSettings

   This class implements the framework required to propagate the rotational dynamics of a body. The settings are constructed as follows:

   .. code-block:: cpp

      RotationalStatePropagatorSettings< StateScalarType >( torqueModelMap,
                                                            bodiesToIntegrate,
                                                            initialBodyStates,
                                                            terminationSettings,
                                                            propagator,
                                                            dependentVariablesToSave )

   where:

   - ``torqueModelMap``

      :class:`TorqueModelMap` List of torque models that are to be used in propagation.

   - :literal:`bodiesToIntegrate`

      :literal:`std::vector< std::string >` that contains the names of the bodies to integrate which must match with those in the :class:`BodyMap`.

   - :literal:`initialBodyStates`

      :literal:`Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >` that stores the states of the bodies to propagate with respect to their central bodies. 

   - :literal:`terminationSettings`

      :literal:`boost::shared_ptr< PropagationTerminationSettings >` that defines the termination settings of the propagation. See the tutorial on :class:`PropagationTerminationSettings` for more details on this class.

   - :literal:`propagator`

      :class:`RotationalPropagatorType` which defines the type of propagator to be used. Currently, the following propagators are supported: 

         - :literal:`quaternions`
         - :literal:`modified_rodrigues_parameters`
         - :literal:`exponential_map`

      By default, the :literal:`quaternions` propagator is used.

   - :literal:`dependentVariablesToSave`

      :literal:`boost::shared_ptr< DependentVariableSaveSettings >` that presents a list of the dependent variables to save during propagation. How this is exactly done is explained below. By default, an empty list is used and no dependent variable is saved. See the tutorial on :class:`DependentVariableSaveSettings` for more details on this class.

.. class:: MassPropagationSettings

    This class implements the framework required to propagate the mass of a body. The constructor of this derived class is overloaded allowing either a single mass-rate per body or multiple mass-rates per body: 

    .. method:: Single mass-rate model per body

        .. code-block:: cpp

            MassPropagationSettings< StateScalarType >( bodiesWithMassToPropagate,
                                                        massRateModels,
                                                        initialBodyMasses,
                                                        terminationSettings,
                                                        dependentVariablesToSave )

        where:

        - :literal:`bodiesWithMassToPropagate`

            :literal:`std::vector< std::string >` that provides the names of the bodies with mass that must be propagated. These names must match with those in the :class:`BodyMap`.

        - :literal:`massRateModels`

            :literal:`std::map< std::string, boost::shared_ptr< MassRateModel > >` that associates a :class:`MassRateModel` to every body with mass that needs to be propagated.

        - :literal:`initialBodyMasses`

            :literal:`Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >` passed by reference that associates an initial body mass to each body with mass to be propagated.

    .. method:: Various mass-rate models per body

        .. code-block:: cpp

            MassPropagationSettings< StateScalarType >( bodiesWithMassToPropagate,
                                                        massRateModels,
                                                        initialBodyMasses,
                                                        terminationSettings,
                                                        dependentVariablesToSave )

        where:

        - :literal:`massRateModels`

            :literal:`std::map< std::string, std::vector< boost::shared_ptr< MassRateModel > > >` that associates a :class:`std::vector` of :class:`MassRateModel` to each body with mass to be propagated.

.. class:: CustomStatePropagatorSettings

    This class allows the user to define and propagate its own state derivative function. The constructor of this derived class is overloaded allowing the user to either use a scalar state or vector state:


    .. method:: Using a scalar state
    
        .. code-block:: cpp

            CustomStatePropagatorSettings< StateScalarType, TimeType >( stateDerivativeFunction,
                                                                        initialState,
                                                                        terminationSettings,
                                                                        dependentVariablesToSave )

        where:

        - :literal:`TimeType`
   
            Template argument used to set the precision of the time, in general :literal:`double` is used. For some application where a high precision is required this can be changed to e.g. :literal`long double`. 

        - :literal:`stateDerivativeFunction`

            :literal:`boost::function< StateScalarType( const TimeType , const StateScalarType ) >` that must comply with the requirements discussed in :ref:`tudatFeaturesIntegrators`.

        - :literal:`initialState`

            :literal:`StateScalarType` that stores the initial state.

    .. method:: Using a vector state
    
        .. code-block:: cpp

            CustomStatePropagatorSettings< StateScalarType, TimeType >( stateDerivativeFunction,
                                                                        initialState,
                                                                        terminationSettings,
                                                                        dependentVariablesToSave )

        where:

        - :literal:`stateDerivativeFunction`

            :literal:`boost::function< Eigen::VectorXd( const double , const Eigen::VectorXd ) >` that must comply with the requirements discussed in :ref:`tudatFeaturesIntegrators`.

        - :literal:`initialState`

            :literal:`Eigen::VectorXd` that stores the initial state.

.. class:: MultiTypePropagatorSettings

    This class is used to propagate multiple types of :class:`PropagatorSettings` concurrently. The constructor of this class is overloaded depending on how the list of propagator settings is passed:

    .. method:: Using an std::vector

        .. code-block:: cpp

            MultiTypePropagatorSettings< StateScalarType >( propagatorSettingsMap,
                                                           terminationSettings,
                                                           dependentVariablesToSave )

        where:
   
        - :literal:`propagatorSettingsMap`

            :literal:`std::vector< boost::shared_ptr< PropagatorSettings< StateScalarType > > >` where each element contains a pointer to a :class:`PropagatorSettings` class. This class is the simplest to use, since it allows to pass a set of unsorted :class:`PropagatorSettings` derived classes by means of the :literal:`push_back` method of :literal:`std::vector`.

    .. method:: Using an std::map

        .. code-block:: cpp

            MultiTypePropagatorSettings< StateScalarType >( propagatorSettingsMap,
                                                            terminationSettings,
                                                            dependentVariablesToSave )

        where:

        - :literal:`propagatorSettingsMap`

            :literal:`std::map< IntegratedStateType, std::vector< boost::shared_ptr< PropagatorSettings< StateScalarType > > > >` where each element contains a pointer to a :class:`PropagatorSettings` class. This class requires a sorted list :class:`PropagatorSettings` derived classes.

   
   .. Warning:: When using the :class:`MultiTypePropagatorSettings` derived class note that the :literal:`dependentVariablesToSave` need to be passed in this constructor and not inside the :literal:`propagatorSettingsMap` since these will be ignored. 

.. class:: MultiArcPropagatorSettings

    This class is meant to be used together with a :class:`MultiArcDynamicsSimulator`. This allows the numerical propagation to be performed in an arc-wise manner. Dynamical model settings may be defined differently per arc. 

   .. code-block:: cpp

      MultiArcPropagatorSettings< StateScalarType >( singleArcSettings,
                                                     transferInitialStateInformationPerArc)

   where:

   - ``singleArcSettings``

      ``std::vector< boost::shared_ptr< SingleArcPropagatorSettings< StateScalarType > > >`` defines the settings for the constituent arcs. The switch times for the arcs are defined by the initial times for each of the arcs. 

   - ``transferInitialStateInformationPerArc``

      ``bool`` allows only a single initial state to be defined: that for the first arc. When this variable is true, the initial state for arc 2 is taken from interpolating arc 1 at the arc 2 start time. This allows a continuous state to be set, while still using the multi-arc interface (for instance for a first estimate when doing multi-arc propagation).

.. tip:: Please beware that all the classes belonging to Tudat libraries are declared above without their namespace. To get the code working please make use of the appropriate :literal:`#include` and :literal:`using` statements.

