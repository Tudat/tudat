.. _tudatFeaturesThrustModels:

Thrust Guidance
===============
This page deals with the inclusion of a thrust force into the dynamical model. Note that when using thrust, it may often be desirable to include the mass in the state vector. A good example of incorporating mass propagation is given in :ref:`walkthroughsUseOfThrustThrustForceAlongVelocityVector`. Details of combining mass and state propagation can be found on the :ref:`tudatFeaturesPropagatorSettings` page.

In Tudat, we define the thrust force by two separate types of settings (which may or may not be linked):

    - The direction of the thrust.
    - The magnitude of the thrust.

In fact, when creating settings for a thrust force, the user needs to provide settings for these two aspects of the force model

.. code-block:: cpp
    
        std::shared_ptr< ThrustDirectionGuidanceSettings > thrustDirectionSettings;
        std::shared_ptr< ThrustMagnitudeSettings > thrustMagnitudeSettings;

        SelectedAccelerationMap accelerationSettingsMap;
        accelerationSettingsMap[ "Vehicle" ][ "Vehicle" ].push_back( std::make_shared< ThrustAccelerationSettings >( thrustDirectionSettings, thrustMagnitudeSettings ) ); 

In the above code snippet, two things may stand out. First of all, we define the thrust force as one that the vehicle exerts on itself. Secondly, to define the thrust force, the user must provide two objects: one of type (derived from) :class:`ThrustDirectionGuidanceSettings` and :class:`ThrustMagnitudeSettings`. The settings are used to create a :class:`ThrustAcceleration` acceleration object. 

.. class:: ThrustAcceleration

   Class contaning the properties of the thrust acceleration. Set by the settings classes described below.
   

Thrust direction
~~~~~~~~~~~~~~~~

For the direction of the thrust, there are presently four available types of guidance. As is done for the acceleration models, some of the types of thrust direction require a specific derived class of :class:`ThrustDirectionGuidanceSettings`, while others are defined purely by their type.

.. class:: ThrustDirectionGuidanceSettings

   Base class for the setting for the direction of the thrust, in most cases not created directly by user. Important exception: thrust direction from existing orientation (see below).

        .. note:: The thrust-direction in the inertial frame is defined by this class. The direction of the thrust in the body-fixed frame (which may be required for the computation of the thrust) is defined in the :class:`ThrustMagnitudeSettings` class.

.. class:: ThrustDirectionFromStateGuidanceSettings

   In various simplified cases, the thrust direction can be assumed to be in line with either the position or velocity w.r.t. some body. As an example: 

   .. code-block:: cpp

 	ThrustDirectionFromStateGuidanceSettings( 
		centralBodyName, isThurstInVelocityDirection, directionIsOppositeToVector );

- :literal:`centralBodyName`

      :literal:`std::string` which defines the name of the central body. Note that this need not be the same body as the center of propagation. For instance, the propagation of a vehicle may be done w.r.t. the Sun, while the thrust direction is computed from the state w.r.t. the Earth.

- :literal:`isThurstInVelocityDirection`

      :literal:`bool` defining whether the thrust is colinear with velocity (if true) or position (if false) w.r.t. the central body.
   
- :literal:`directionIsOppositeToVector`

        :literal:`bool` defining whether the thrust force is in the same direction, or opposite to the direction, of the position/velocity w.r.t. the central body.

.. class:: CustomThrustDirectionSettings

   For a generalized thrust direction guidance, the thrust can be defined as an arbitrary function of time. This allows a broad range of options to be defined, at the expense of increased complexity (somehow this thrust direction needs to be manually defined):

   .. code-block:: cpp

 	CustomThrustDirectionSettings( 
		thrustDirectionFunction );


- :literal:`thrustDirectionFunction`

        A :literal:`std::function< Eigen::Vector3d( const double ) >` returning a the thrust direction in the inertial frame as an ``Eigen::Vector3d`` (which should be of unit norm!) as a function of a ``double`` (representing time). Details on how to create such an :literal:`std::function` are given on :ref:`externalPointersExamples`. 

As a possible example of how to use this function:

   .. code-block:: cpp

      Eigen::Vector3d myThrustFunction( const double time, const NamedBodyMap& bodyMap )
      {
         Eigen::Vector3d thrustDirection = 
	    //... define algorithm to compute thurst based on current environment

         return thrustDirection;
      }

      int main( )
      {
	 // Define environment and other settings
         NamedBodyMap bodyMap = ...
 
         // ...
         // ...
         // ...

         // Create custom thrust direction settings
         std::make_shared< CustomThrustDirectionSettings >( 
            std::bind( &myThrustFunction, std::placeholders::_1, bodyMap );
      }

.. class:: CustomThrustOrientationSettings

As an alternative expression for generalized thrust direction guidance, the thrust orientation can be defined as an arbitrary function of time. As with the custom thrust direction. this allows a broad range of options to be defined, at the expense of increased complexity). 

 .. code-block:: cpp

 	CustomThrustOrientationSettings( 
		thrustOrientationFunction );


- :literal:`thrustOrientationFunction`

        A :literal:`std::function< Eigen::Quaterniond( const double ) >` returning the rotation from body-fixed to inertial state, represented as an ``Eigen::Quaterniond`` (which should be of unit norm!) as a function of a ``double`` (representing time). See :class:`CustomThrustDirectionSettings` for an example on how to realize this construction.

.. class:: MeeCostateBasedThrustDirectionSettings

By using these settings for the thrust direction, the so-called co-states of the Modified Equinoctial elements are used to determine the direction of the thrust. Details of this model are given by Kluever (2010), Boudestijn (2014) and Hogervorst (2017). 

.. method:: Thrust direction from existing orientation

   In some cases, the vehicle's orientation may be predetermined, either due to aerodynamic guidance of concurrent propagation of the rotational equations of motion. In such a case, the thrust direction is computed from the body-fixed thrust direction (defined in :class:`ThrustMagnitudeSettings`) and the existing vehicle orientation. This thrust direction does not require a specific derived class, but instead only requires the input of :literal:`thrust_direction_from_existing_body_orientation` to the :class:`ThrustDirectionGuidanceSettings` constructor, so:

 .. code-block:: cpp

 	ThrustDirectionGuidanceSettings( 
		thrust_direction_from_existing_body_orientation, "" );

Thrust magnitude
~~~~~~~~~~~~~~~~
To define the thrust magnitude, there are presently three available types of settings, each with its own dedicated derived class of :class:`ThrustMagnitudeSettings`. We note that presently, the definition of the thrust direction in the body-fixed frame is also defined through these derived classes. In essence, the :class:`ThrustMagnitudeSettings` defines all local (to the vehicle systems) settings for the thrust, while :class:`ThrustDirectionGuidanceSettings` defines how the full vehicle must orient itself in space for the required thrust direction to be achieved. At present, there is no direct option for thrust-vector control (i.e. modifying the thrust direction in the body-fixed frame). If your application requires such functionality, please contact the Tudat support team. The following thrust magnitude settings are available:

.. class:: ThrustMagnitudeSettings

   Base class for the thrust magnitude settings. 

.. class:: ConstantThrustMagnitudeSettings

   This model defines a constant thrust force and specific impulse:

 .. code-block:: cpp

    ConstantThrustMagnitudeSettings(
	thrustMagnitude, specificImpulse, bodyFixedThrustDirection ):
   
- :literal:`thrustMagnitude` Constant thrust magnitude to use (in Newtons, as ``double``)

- :literal:`specificImpulse` Specific impulse to use for the thrust, as ``double``. This quantity is used when applying a mass rate model in the propagation the vehicle dynamics, relating the thrust to the mass decrease of the vehicle.

- :literal:`bodyFixedThrustDirection` Body-fixed thrust direction (positive x-direction by default) as an ``Eigen::Vector3d``. Note that this should be a unit-vector, and represents the direction opposite to the nozzle direction.

.. class:: FromFunctionThrustMagnitudeSettings

   This model defines a thrust force and specific impulse that can vary with time. It requires the following settings as input:

 .. code-block:: cpp

   FromFunctionThrustMagnitudeSettings(
            thrustMagnitudeFunction, specificImpulseFunction,
            isEngineOnFunction, bodyFixedThrustDirection );
    
- :literal:`thrustMagnitudeFunction` Thrust magnitude to use (in Newtons), as a function of time as a function of time, defined by an `std::function< double( const double ) >` 

- :literal:`specificImpulseFunction` Specific impulse to use for the thrust as a function of time, defined by an `std::function< double( const double ) >`. This quantity is used when applying a mass rate model in the propagation the vehicle dynamics, relating the thrust to the mass decrease of the vehicle.

- :literal:`isEngineOnFunction` Boolean defining whether the thrust should be engaged at all, as a function of time (e.g. thrust is 0 N if it returns false), defined by an `std::function< bool( const double ) >`. By default this function always returns true.

- :literal:`bodyFixedThrustDirection` Body-fixed thrust direction (positive x-direction by default) as an ``Eigen::Vector3d``. Note that this should be a unit-vector, and represents the direction opposite to the nozzle direction.

Note that if you wish to use a constant value (as opposed to the :literal:`std::function` expression) for any or all of the first three arguments, you can use lambda expression. For instance, for a variable thrust magnitude, but constant specific impulse of 300 s, you can use:

 .. code-block:: cpp

   FromFunctionThrustMagnitudeSettings(
            thrustMagnitudeFunction, []( const double ){ return 300; } );

Where the last two arguments (isEngineOnFunction and bodyFixedThrustDirection) are omitted and therefore are set to their default values.

.. class:: FromBodyThrustMagnitudeSettings

   A :class:`Body` object may be endowed with a :class:`VehicleSystems` property, which defines the suite of hardware that it carries. One of the systems that may be defined in the :class:`VehicleSystems` object which is a list of :class:`EngineModel` objects (stored in a list). This class creates thrust magnitude settings that use the thrust from one or all of the :class:`EngineModel` objects that a vehicle is endowed with. In such a situation, the thrust direction, force and specific impulse are taken from the :class:`EngineModel`. It requires the following settings:

 .. code-block:: cpp

   FromBodyThrustMagnitudeSettings(
            useAllEngines, thrustOrigin );

- :literal:`useAllEngines` A ``bool`` defining whether all engines (i.e. all entries in the :literal:`engineModels` member of the :class:`VehicleSystems` object in the vehicle's :class:`Body` object) are to be used, or only a single engine.

- :literal:`thrustOrigin` A ``std::string`` with the name of the engine that is to be used for the thrust (to be empty if useAllEngines is true)

.. class:: EngineModel

   Class in which the model of the engine is saved. Settings similar to those in the previous two thrust magnitude settings may be stored. However, using the interface with an engine model allows a more integrated systems/trajectory simulation to be performed, with applications in e.g. MDO. It allows multiple engine models, each with their own properties, to be defined. 

   .. note:: Currently this class in under development, note that this is still a priliminary version.
