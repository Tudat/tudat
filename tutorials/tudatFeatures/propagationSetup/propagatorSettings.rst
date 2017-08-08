.. _tudatFeaturesPropagatorSettings:

Propagator Settings
===================
This page presents an overview of the available options within :class:`PropagatorSettings`. As the name suggests, these settings define how the orbit is propagated within the :class:`DynamicsSimulator`.

Available settings
~~~~~~~~~~~~~~~~~~
Similarly to the :class:`IntegratorSettings` discussed in :ref:`tudatFeaturesIntegratorSettings`, various derived classes are used to implement different settings:

.. class:: TranslationalStatePropagatorSettings

    This class implements the framework required to propagate the translation state of a body. The constructor of this derived class is overloaded allowing two types of termination conditions:

    - Termination when the simulation time reaches a predefined :literal:`endTime` (Default).
    - Termination when a predefined dependent variables meets a certain criterion.

    **Default termination settings**

    .. code-block:: cpp

        TranslationalStatePropagatorSettings( &centralBodies,
                                              &accelerationsMap,
                                              &bodiesToIntegrate,
                                              &initialBodyStates,
                                              endTime,
                                              propagator,
                                              dependentVariablesToSave)

    where:

    - :literal:`&centralBodies` is a :literal:`const std::vector< std::string >` passed by reference that defines the names of the central bodies which must match with those in the :class:`BodyMap`.
    - :literal:`&accelerationsMap` is a :literal:`const basic_astrodynamics::AccelerationMap` passed by reference that contains the accelerations for each body as discussed in :ref:`tudatFeaturesAccelerationIndex`.
    - :literal:`&bodiesToIntegrate` is a :literal:`const std::vector< std::string >` passed by reference that defines the names of the bodies to integrate which must match with those in the :class:`BodyMap`.
    - :literal:`&initialBodyStates` is a :literal:`const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >` passed by reference 
    - :literal:`endTime` is a :literal:`double` that defines the end-time of the simulation.
    - :literal:`propagator` is a :literal:`const TranslationalPropagatorType` which defines the type of propagator being used. Currently, :literal:`cowell` and :literal:`encke` are available. By default, the :literal:`cowell` propagator is used.
    - :literal:`dependentVariablesToSave` is a :literal:`const boost::shared_ptr< DependentVariableSaveSettings >` that presents a list of the dependent variables to save during propagation. How this is exactly done is explained below. By default, an empty list is used and no dependent variable is saved.

    **User-defined termination settings**

    .. code-block:: cpp

        TranslationalStatePropagatorSettings( &centralBodies,
                                              &accelerationsMap,
                                              &bodiesToIntegrate,
                                              &initialBodyStates,
                                              terminationSettings,
                                              propagator,
                                              dependentVariablesToSave,

    where:

    - :literal:`terminationSettings` is a :literal:`boost::shared_ptr< PropagationTerminationSettings >` that defines the termination settings of the propagation. This is the fifth argument and replaces the :literal:`endTime` in the default constructor.

.. class:: MassPropagatorSettings

    This class implements the framework required to propagate the mass of a body. The constructor of this derived class is overloaded allowing either a single mass-rate per body or multiple mass-rates per body:

    **Single mass-rate model per body**

    .. code-block:: cpp

        MassPropagatorSettings(
                bodiesWithMassToPropagate,
                massRateModels,
                &initialBodyMasses,
                terminationSettings,
                dependentVariablesToSave )

    where:

    - :literal:`bodiesWithMassToPropagate` is a :literal:`const std::vector< std::string >` that provides the names of the bodies with mass that must be propagated. These names must match with those in the :class:`BodyMap`.
    - :literal:`massRateModels` is a :literal:`const std::map< std::string, boost::shared_ptr< basic_astrodynamics::MassRateModel > >` that associates a :class:`MassRateModel` to a every body with mass that needs to be propagated.
    - :literal:`&initialBodyMasses` is a :literal:`const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >` passed by reference that associates an initial body mass to each body with mass to be propagated.

    **Various mass-rate models per body**

    .. code-block:: cpp

        MassPropagatorSettings(
                bodiesWithMassToPropagate,
                massRateModels,
                &initialBodyMasses,
                terminationSettings,
                dependentVariablesToSave )

    where:

    - :literal:`massRateModels` is a :literal:`const std::map< std::string, std::vector< boost::shared_ptr< basic_astrodynamics::MassRateModel > > >` that associates a :literal:`std::vector` of :class:`MassRateModel` to each body with mass to be propagated.

.. class:: CustomPropagatorSettings

.. class:: MultiTypePropagatorSettings

Propagation saving
~~~~~~~~~~~~~~~~~~

Propagation termination conditions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
