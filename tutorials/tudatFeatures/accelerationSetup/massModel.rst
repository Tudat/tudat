.. _tudatFeaturesFrameworkMassRateModelSetup:

Mass Rate Model Set-up
======================
Although propagating a body's translational dynamics is the backbone of Tudat's simulations, it is also possible to propagate a vehicle's mass (either concurrently or separately). The manner in which the models that govern the 'mass dynamics', i.e. mass-rate models, are handled in the code is very similar to the acceleration models: a list of settings for the models is created by the user, which are then used to create the required objects. The list to be created by the user is:

.. code-block:: cpp

    std::map< std::string, std::vector< std::shared_ptr< MassRateModelSettings > > > massRateModelSettings;

where the map key denotes the body of which the mass-rate is to be computed.

.. class:: MassRateModelSettings

   Base class for the mass rate model setup. Currently two mass rate models are available each with its own derived class described below.

.. class:: CustomMassRateModelSettings

   Using this class, the user must provide a :literal:`std::function< double( const double ) > function`, i.e. a function returning a double, representing the mass-rate, and taking another double, representing time, as an input. The internal workings of this function are completely up to the user. If any help is required in setting up such a model please contact the Tudat support team.

.. class:: FromThrustMassModelSettings

   Using this mass-rate model, the change in vehicle mass due to the expulsion of propellant is taken into account when propagating a vehicle's dynamics. It retrieves the required data from a :class:`ThrustAcceleration` object (set by :class:`ThrustAccelerationSettings`), ensuring full consistency between the two. Two option are available when creating this type of mass-rate model:

   - Use all thrust forces acting on a single body, combined into a single mass-rate model. This will in most cases be the model of choice, as there is often no need to distinguish between thrust sources when computing the mass rate: only the total amount of propellant usage is relevant. This option is toggled by setting the :literal:`useAllThrustModels` input argument of the :class:`FromThrustMassModelSettings` constructor to true.
   - Use a single thrust model, defined by a string-identifier. When creating a thrust model, a :literal:`thrustOriginId` input is provided to the :class:`ThrustMagnitudeSettings` settings constructor. Only in the :class:`FromBodyThrustMagnitudeSettings` derived class is this thrust origin id set to anything else than an empty string: it represents the engine name.



