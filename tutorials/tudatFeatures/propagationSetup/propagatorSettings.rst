.. _tudatFeaturesPropagatorSettings:

Propagator Settings: Basics
===========================
This page presents an overview of the available options within :class:`PropagatorSettings`. As the name suggests, these settings define how the orbit is propagated within the :class:`DynamicsSimulator`.

Similarly to the :class:`IntegratorSettings` discussed in :ref:`tudatFeaturesIntegratorSettings`, various derived classes are used to implement different settings:

.. class:: TranslationalStatePropagatorSettings

    This class implements the framework required to propagate the translation state of a body. The constructor of this derived class is overloaded allowing two types of termination conditions:

    - Termination when the simulation time reaches a predefined :literal:`endTime` (Default).
    - Termination when a predefined dependent variables meets a certain criterion.

    .. method:: Default termination settings

        .. code-block:: cpp

            TranslationalStatePropagatorSettings( centralBodies,
                                                  accelerationsMap,
                                                  bodiesToIntegrate,
                                                  initialBodyStates,
                                                  endTime,
                                                  propagator,
                                                  dependentVariablesToSave)

        where:

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

            :class:`TranslationalPropagatorType` which defines the type of propagator being used. Currently, :literal:`cowell` and :literal:`encke` are available. By default, the :literal:`cowell` propagator is used.

        - :literal:`dependentVariablesToSave`

            :literal:`boost::shared_ptr< DependentVariableSaveSettings >` that presents a list of the dependent variables to save during propagation. How this is exactly done is explained below. By default, an empty list is used and no dependent variable is saved.

        .. note:: The state variables contained in :literal:`initialBodyStates` are ordered with respect to the elements of :literal:`centralBodies` and :literal:`bodiesToIntegrate`. Please take a look at the following pseudocode:

            .. code-block:: cpp

                centralBodies = { Sun , Earth , Moon }
                bodiesToIntegrate = { Earth , Moon }
                initialBodyStates = { xEarthWrtSun , yEarthWrtSun , zEarthWrtSun , uEarthWrtSun , vEarthWrtSun , wEarthWrtSun , 
                                      xMoonWrtEarth , yMoonWrtEarth , zMoonWrtEarth , uMoonWrtEarth , vMoonWrtEarth , wMoonWrtEarth }
            

    .. method:: User-defined termination settings

        .. code-block:: cpp

            TranslationalStatePropagatorSettings( centralBodies,
                                                  accelerationsMap,
                                                  bodiesToIntegrate,
                                                  initialBodyStates,
                                                  terminationSettings,
                                                  propagator,
                                                  dependentVariablesToSave )

        where:

        - :literal:`terminationSettings`

            :literal:`boost::shared_ptr< PropagationTerminationSettings >` that defines the termination settings of the propagation. This is the fifth argument and replaces the :literal:`endTime` in the default constructor.

.. class:: MassPropagatorSettings

    This class implements the framework required to propagate the mass of a body. The constructor of this derived class is overloaded allowing either a single mass-rate per body or multiple mass-rates per body: 

    .. method:: Single mass-rate model per body

        .. code-block:: cpp

            MassPropagatorSettings(
                    bodiesWithMassToPropagate,
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

            MassPropagatorSettings(
                    bodiesWithMassToPropagate,
                    massRateModels,
                    initialBodyMasses,
                    terminationSettings,
                    dependentVariablesToSave )

        where:

        - :literal:`massRateModels`

            :literal:`std::map< std::string, std::vector< boost::shared_ptr< MassRateModel > > >` that associates a :class:`std::vector` of :class:`MassRateModel` to each body with mass to be propagated.

.. class:: CustomPropagatorSettings

    This class allows the user to define and propagate its own state derivative function. The constructor of this derived class is overloaded allowing the user to either use a scalar state or vector state:


    .. method:: Using a scalar state
    
        .. code-block:: cpp

            CustomStatePropagatorSettings(
                stateDerivativeFunction,
                initialState,
                terminationSettings,
                dependentVariablesToSave )

        where:

        - :literal:`stateDerivativeFunction`

            :literal:`boost::function< double( const double , const double ) >` that must comply with the requirements discussed in :ref:`tudatFeaturesIntegrators`.

        - :literal:`initialState`

            :literal:`double` that stores the initial state.

    .. method:: Using a vector state
    
        .. code-block:: cpp

            CustomStatePropagatorSettings(
                stateDerivativeFunction,
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

            MultiTypePropagatorSettings(
                propagatorSettingsMap,
                terminationSettings,
                dependentVariablesToSave )

        where:
   
        - :literal:`propagatorSettingsMap`

            :literal:`std::vector< boost::shared_ptr< PropagatorSettings< StateScalarType > > >` where each element contains a pointer to a :class:`PropagatorSettings` class. This class is the simplest to use, since it allows to pass a set of unsorted :class:`PropagatorSettings` derived-classes by means of the :literal:`push_back` method of :literal:`std::vector`.

    .. method:: Using an std::map

        .. code-block:: cpp

            MultiTypePropagatorSettings(
                propagatorSettingsMap,
                terminationSettings,
                dependentVariablesToSave )

        where:

        - :literal:`propagatorSettingsMap`

            :literal:`std::map< IntegratedStateType, std::vector< boost::shared_ptr< PropagatorSettings< StateScalarType > > > >` where each element contains a pointer to a :class:`PropagatorSettings` class. This class requires a sorted list :class:`PropagatorSettings` derived-classes.

.. class:: MultiArcPropagatorSettings

    This class is meant to be used together with a :class:`MultiArcDynamicsSimulator`. At the moment, the multi-arc simulator elements in Tudat are undergoing testing and are thus not yet available.

.. tip:: Please beware that all the classes belonging to Tudat libraries are declared above without their namespace. To get the code working please make use of the appropriate :literal:`#include` and :literal:`using` statements.

