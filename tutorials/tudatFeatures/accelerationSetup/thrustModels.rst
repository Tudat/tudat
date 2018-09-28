.. _tudatFeaturesThrustModels:

Thrust Guidance
=============
When including thrust in the orbit propagation, it may often be desirable to modify the thrust properties of the vehicle as a function of the current environment. Depending on the settings of the system, the magnitude and direction of this force can be controlled to within certain bounds. How the thrust guidance is setup is the topic of this page.

In Tudat, we define the thrust force by two separate types of settings (which may or may not be linked):

    - The direction of the thrust.
    - The magnitude of the thrust.

In fact, when creating settings for a thrust force, the user needs to provide settings for these two aspects of the force model

.. code-block:: cpp
    
        std::shared_ptr< ThrustDirectionGuidanceSettings > thrustDirectionSettings;
        std::shared_ptr< ThrustEngineSettings > thrustMagnitudeSettings;

        SelectedAccelerationMap accelerationSettingsMap;
        accelerationSettingsMap[ "Vehicle" ][ "Vehicle" ].push_back( std::make_shared< ThrustAccelerationSettings >( thrustDirectionSettings, thrustMagnitudeSettings ) ); 

In the above code snippet, two things may stand out. First of all, we define the thrust force as one that the vehicle exerts on itself. Secondly, to define the thrust force, the user must provide two objects: one of type (derived from) :class:`ThrustDirectionGuidanceSettings` and :class:`ThrustEngineSettings`. The settings are used to create a :class:`ThrustAcceleration` acceleration object. 

.. class:: ThrustAcceleration

   Class contaning the properties of the thrust acceleration. Set by the settings classes described below.
   

Thrust direction
~~~~~~~~~~~~~~~~

For the direction of the thrust, there are presently four available types of guidance. As is done for the acceleration models, some of the types of thrust direction require a specific derived class of :class:`ThrustDirectionGuidanceSettings`, while others are defined purely by their type.

.. class:: ThrustDirectionGuidanceSettings

   Base class for the setting for the direction of the thrust.

.. class:: ThrustDirectionFromStateGuidanceSettings

   In various cases, the thrust direction is to be in line with the position or velocity w.r.t. some body. The user must specify:
   
        - The central body w.r.t. which the state is to be computed. For instance, the propagation of a vehicle may be done w.r.t. the Sun, while the thrust direction is computed from the state w.r.t. the Earth.
        - Whether the thrust is colinear with velocity or position w.r.t. the central body.
        - Whether the thrust force is in the same direction, or opposite to the direction, of the state of the vehicle w.r.t. the central body.

.. class:: CustomThrustDirectionSettings

   For a generalized thrust direction guidance, the thrust can be defined as an arbitrary function of time. This allows a broad range of options to be defined, at the expense of increased complexity (somehow this thrust direction needs to be manually defined). The user must specify:

        - A function returning a ``Eigen::Vector3d`` (which should be of unit norm!) as a function of a ``double`` (representing time). If any help is required in defining this function for your specific application, don't hesitate to contact the Tudat support team.

.. class:: CustomThrustOrientationSettings

   As an alternative expression for generalized thrust direction guidance, the thrust orientation can be defined as an arbitrary function of time. As with the custom thrust direction. this allows a broad range of options to be defined, at the expense of increased complexity). The user must specify:

        - A function returning a ``Eigen::Quateriond`` as a function of a ``double`` (representing time). The quaternion defines the rotation from the body-fixed to the propagation frame. If any help is required in defining this function for your specific application, don't hesitate to contact the Tudat support team.

        .. note:: The direction of the thrust in the body-fixed frame (which is required for the computation of the thrust) is defined in the :class:`ThrustEngineSettings` class.

.. method:: Thrust direction from existing orientation

   In some cases (discussed in more detail below), the vehicle's orientation may be predefined, either due to aerodynamic guidance of concurrent propagation of the rotational equations of motion. In such a case, the thrust direction is computed from the body-fixed thrust direction and the existing vehicle orientation. This thrust direction does not require a specific derived class, but instead only requires the input of :literal:`thrust_direction_from_existing_body_orientation` to the :class:`ThrustDirectionGuidanceSettings` constructor.

Thrust magnitude
~~~~~~~~~~~~~~~~
To define the thrust magnitude, there are presently three available types of settings, each with its own dedicated derived class of :class:`ThrustEngineSettings`. We note that presently, the definition of the thrust direction in the body-fixed frame is also defined through these derived classes. In essence, the :class:`ThrustEngineSettings` defines all local (to the vehicle systems) settings for the thrust, while :class:`ThrustDirectionGuidanceSettings` defines how the full vehicle must orient itself in space for the required thrust direction to be achieved. At present, there is no option for thrust-vector control (i.e. modifying the thrust direction in the body-fixed frame). If your application requires such functionality, please contact the Tudat support team. The following thrust magnitude settings are available:

.. class:: ThrustEngineSettings

   Base class for the thrust magnitude settings. 

.. class:: ConstantThrustEngineSettings

   This model defines a constant thrust force. It requires the following settings as input:
   
      - Thrust magnitude to use (in Newtons)
      - Specific impulse to use for the thrust. This quantity is used when applying a mass rate model in the propagation the vehicle dynamics, relating the thrust to the mass decrease of the vehicle.
      - Body-fixed thrust direction (positive x-direction by default). Note that this should be a unit-vector, and represents the direction opposite to the nozzle direction.

.. class:: FromFunctionThrustEngineSettings

   This model defines a thrust force and specific impulse that can vary with time. It requires the following settings as input:
    
        - Thrust magnitude to use (in Newtons), as a function of time as a function of time.
        - Specific impulse to use for the thrust, as a function of time. This quantity is used when applying a mass rate model in the propagation the vehicle dynamics, relating the thrust to the mass decrease of the vehicle.
        - Body-fixed thrust direction (positive x-direction by default). Note that this should be a unit-vector, and represents the direction opposite to the nozzle direction.

.. class:: FromBodyThrustEngineSettings

   A :class:`Body` object may be endowed with a :class:`VehicleSystems` property, which defines the suite of hardware that it carries. One of the systems that may be defined in the :class:`VehicleSystems` object which is a list of :class:`EngineModel` objects (stored in a list). We provide thrust engine settings that use the thrust from one or all of the :class:`EngineModel` objects that a vehicle is endowed with, in this class. It requires the following settings:

      - A boolean defining whether all engines (i.e. all entries in the :literal:`engineModels` member of the :class:`VehicleSystems` object in the vehicle's :class:`Body` object.
      - The name of the engine that is to be used for the thrust (to be empty if all engines are used)

.. class:: EngineModel

   Class in which the model of the engine is saved. Settings similar to those in the previous two thrust magnitude settings may be stored. However, using the interface with an engine model allows a more integrated systems/trajectory simulation to be performed, with applications in e.g. MDO. It allows multiple engine models, each with their own properties, to be defined. 

   .. note:: Currently this class in under development, note that this is still a priliminary version.
