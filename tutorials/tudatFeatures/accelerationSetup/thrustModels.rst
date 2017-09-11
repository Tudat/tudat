.. _tudatFeaturesThrustModels:

Thrust Models
=============
When simulating the dynamics of solar system bodies, we defined the environment and types of accelerations, and let the numerical integration produce a solution to the resulting equations of motion. When considering manmade vehicles, however, we often have direct influence on the forces acting on our body. In Tudat, controlling such forces can be done in a number of manners, which we will discuss in this tutorial. At present, the forces that can be controlled are:

    - **Aerodynamic acceleration:** For aerodynamic coefficients that are defined directly as a function of vehicle orientataion (angle of attack and/or sideslip angle), these vehicle orientation angles can be modulated, thereby altering the forces acting on the vehicle.
    - **Thrust acceleration:** By using a propulsion system, a thrust force can be exerted on the vehicle. Depending on the settings of the system, the magnitude and direction of this force can be controlled to within certain bounds.

Thrust Control
~~~~~~~~~~~~~~
In Tudat, we define the thrust force by two separate types of settings (which may or may not be linked):

    - The direction of the thrust.
    - The magnitude of the thrust.

In fact, when creating settings for a thrust force, the user needs to provide settings for these two aspects of the force model

.. code-block:: cpp
    
        boost::shared_ptr< ThrustDirectionGuidanceSettings > thrustDirectionSettings;
        boost::shared_ptr< ThrustEngineSettings > thrustMagnitudeSettings;

        SelectedAccelerationMap accelerationSettingsMap;
        accelerationSettingsMap[ "Vehicle" ][ "Vehicle" ].push_back( boost::make_shared< ThrustAccelerationSettings >( thrustDirectionSettings, thrustMagnitudeSettings ) ); 

In the above code snippet, two things may stand out. First of all, we define the thrust force as one that the vehicle exerts on itself. Secondly, to define the thrust force, the user must provide two objects: one of type (derived from) :class:`ThrustDirectionGuidanceSettings` and :class:`ThrustEngineSettings`.

Thrust direction
****************
For the direction of the thrust, there are presently four available types of guidance. As is done for the acceleration models, some of the types of thrust direction require a specific derived class of :class:`ThrustDirectionGuidanceSettings`, while others are defined purely by their type.

    **Thrust colinear with state segment:**
        In various cases, the thrust direction is to be in line with the position or velocity w.r.t. some body. In Tudat, the :class:`ThrustDirectionFromStateGuidanceSettings` class is used to define such a thrust direction. The user must specify:
   
        - The central body w.r.t. which the state is to be computed. For instance, the propagation of a vehicle may be done w.r.t. the Sun, while the thrust direction is computed from the state w.r.t. the Earth.
        - Whether the thrust is colinear with velocity or position w.r.t. the central body.
        - Whether the thrust force is in the same direction, or opposite to the direction, of the state of the vehicle w.r.t. the central body.

    **Custom thrust direction:**
        For a generalized thrust direction guidance, the thrust can be defined as an arbitrary function of time. This allows a broad range of options to be defined, at the expense of increased complexity (somehow this thrust direction needs to be manually defined). In Tudat, the CustomThrustDirectionSettings class is used to define such a thrust direction. The user must specify:

        - A function returning a Vector3d (which should be of unit norm!) as a function of a double (representing time). If any help is required in defining this function for your specific application, don't hesitate to contact the Tudat support team.

    **Custom thrust orientation:**
        As an alternative expression for generalized thrust direction guidance, the thrust orientation can be defined as an arbitrary function of time. As with the custom thrust direction. this allows a broad range of options to be defined, at the expense of increased complexity). In Tudat, the CustomThrustOrientationSettings class is used to define such a thrust direction. The user must specify:

        - A function returning a Quateriond as a function of a double (representing time). The quaternion defines the rotation from the body-fixed to the propagation frame. If any help is required in defining this function for your specific application, don't hesitate to contact the Tudat support team.

        .. note:: The direction of the thrust in the body-fixed frame (which is required for the computation of the thrust) is defined in the :class:`ThrustEngineSettings` class.

    **Thrust direction from existing orientation**
        In some cases (discussed in more detail below), the vehicle's orientation may be predefined, either due to aerodynamic guidance of concurrent propagation of the rotational equations of motion (not yet implemented). In such a case, the thrust direction is computed from the body-fixed thrust direction and the existing vehicle orientation. This thrust direction does not require a specific derived class, but instead only requires the input of :literal:`thrust_direction_from_existing_body_orientation` to the :class:`ThrustDirectionGuidanceSettings` constructor.

Thrust magnitude
~~~~~~~~~~~~~~~~
To define the thrust magnitude, there are presently three available types of settings, each with its own dedicated derived class of :class:`ThrustEngineSettings`. We note that presently, the definition of the thrust direction in the body-fixed frame is also defined through these derived classes. In essence, the :class:`ThrustEngineSettings` defines all local (to the vehicle systems) settings for the thrust, while :class:`ThrustDirectionGuidanceSettings` defines how the full vehicle must orient itself in space for the required thrust direction to be achieved. At present, there is no option for thrust-vector control (i.e. modifying the thrust direction in the body-fixed frame). If your application requires such functionality, please contact the Tudat support team. The following thrust magnitude settings are available:

    **Constant thrust magnitude:**
        This model defines a constant thrust force, with its settings defined in the :class:`ConstantThrustEngineSettings` class. It requires the following settings as input:
        - Thrust magnitude to use (in Newtons)
        - Specific impulse to use for the thrust. This quantity is used when applying a mass rate model in the propagation the vehicle dynamics, relating the thrust to the mass decrease of the vehicle.
        - Body-fixed thrust direction (positive x-direction by default). Note that this should be a unit-vector, and represents the direction opposite to the nozzle direction.

    **Variable thrust magnitude/specific impulse:**
        This model defines a thrust force and specific impulse that can vary with time. Its settings are defined in the :class:`FromFunctionThrustEngineSettings` class. It requires the following settings as input:
    
        - Thrust magnitude to use (in Newtons),as a function of time as a function of time.
        - Specific impulse to use for the thrust, as a function of time. This quantity is used when applying a mass rate model in the propagation the vehicle dynamics, relating the thrust to the mass decrease of the vehicle.
        - Body-fixed thrust direction (positive x-direction by default). Note that this should be a unit-vector, and represents the direction opposite to the nozzle direction.

    **Thrust from engine settings:**
        A Body object may be endowed with a :class:`VehicleSystems` property, which defines the suite of hardware that it carries. One of the systems that may be defined in the :class:`VehicleSystems` object is a list of :class:`EngineModel` objects (stored in a list). In an object of the :class:`EngineModel` (derived) class, settings similar to those in the previous two thrust magnitude settings may be stored. However, using the interface with an engine model allows a more integrated systems/trajectory simulation to be performed, with applications in e.g. MDO. It allows multiple engine models, each with their own properties, to be defined. We provide thrust engine settings that use the thrust from one or all of the :class:`EngineModel` objects that a vehicle is endowed with, through the :class:`FromBodyThrustEngineSettings` class. It requires the following settings:

        - A boolean defining whether all engines (i.e. all entries in the :literal:`engineModels` member of the VehicleSystems object in the vehicle's Body object.
        - The name of the engine that is to be used for the thrust (to be empty if all engines are used)
