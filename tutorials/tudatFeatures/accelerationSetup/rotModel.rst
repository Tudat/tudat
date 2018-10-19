.. _tudatFeaturesFrameworkTorqueModelSetup:    

Torque Model Set-up
=======================
The setup of the rotational state is similar to the set-up of the translational acceleration framework.

The user can define the settings for the torque acting on a vehicle in the following variable: :literal:`SelectedTorqueMap`.

.. class:: SelectedTorqueMap

   This is a :literal:`typedef` for:

   .. code-block:: cpp

         std::map< std::string, std::map< std::string, std::vector< std::shared_ptr< torqueSettings > > > >
   
   Just as in the translational acceleration setup, the first key is the body undergoing the torque, the second key is the body exerting the torque, and the third entry is a vector
   of :literal:`torqueSettings`. The available torque models are given in an :literal:`enum` called :literal:`AvailableTorque`. An example can be given using the Apollo capsule 
   example:
   
   .. code-block:: cpp

         SelectedTorqueMap selectedTorqueModelMap;
         selectedTorqueModelMap[ "Apollo" ][ "Earth" ].push_back( std::make_shared< TorqueSettings >( aerodynamic_torque ) );

   Currently, there are only two options for the torque settings, both which are not derived classes of the :literal:`TorqueSettings` class, and do not need extra information. Thus both can be defined
   in a similar way.

   When all the settings of the torque acting on the body are defined, the :literal:`TorqueModelMap` can be created as follows:

   .. code-block:: cpp

          TorqueModelMap torqueModelMap = createTorqueModelsMap( bodyMap, selectedTorqueModelMap )

   This will create a map that can be used as an input to the :literal:`RotationalStatePropagatorSettings`.

Available Models
~~~~~~~~~~~~~~~~
Currently, there are two torque models that can be used for the rotational state propagator:

.. method:: Second-order Gravitational Torque

   Torque exerted by a point mass on a body with a degree two spherical harmonics mass distribution. A degree two spherical harmonics mass distribution can be represented by an inertia tensor, thus 
   for this torque model, the body undergoing the torque needs to have an inertia tensor defined. The body exerting the torque only needs to have a point-mass gravitational model defined. 
   
   The settings are created as follows, for a torque exerted by ``"Mars"`` on body ``"Phobos"``:

   .. code-block:: cpp

       SelectedTorqueMap selectedTorqueModelMap;
       selectedTorqueModelMap[ "Apollo" ][ "Earth" ].push_back( std::make_shared< TorqueSettings >( second_order_gravitational_torque ) );

   Requires the following environment models to be defined:

   - Gravity field (at least point-mass) for body exerting torque (set by :class:`GravityFieldSettings`).
   - Inertia tensor of body undergoing torque.
   - Current state of bodies undergoing and exerting torque, either from an Ephemeris model (set by :class:`EphemerisSettings`) or from the numerical propagation.

.. class:: SphericalHarmonicTorqueSettings

   Torque exerted by a point mass on a body with an arbitrary degree/order spherical harmonics mass distribution. The body exerting the torque only needs to have a gravitational model defined. 

   The settings are created with the dedicated ``SphericalHarmonicTorqueSettings`` class. As an example, for a spherical harmonic torque, expanded to degree and order 8, exerted by ``"Mars"`` on body ``"Phobos"``:

   .. code-block:: cpp

       SelectedTorqueMap selectedTorqueModelMap;
       int maximumDegree = 8;
       int maximumOrder = 8;
       selectedTorqueModelMap[ "Apollo" ][ "Earth" ].push_back( std::make_shared< SphericalHarmonicTorqueSettings >( maximumDegree, maximumOrder ) );

   Requires the following environment models to be defined:

   - Gravity field (at least point-mass) for body exerting torque (set by :class:`GravityFieldSettings`).
   - Spherical harmonic gravity field for body undergoing torque (set by :class:`SphericalHarmonicsGravityFieldSettings`)..
   - Current state of bodies undergoing and exerting torque, either from an Ephemeris model (set by :class:`EphemerisSettings`) or from the numerical propagation.

.. method:: Aerodynamic Torque

   Torque exerted by a body with an atmosphere model and shape model on another body. It is important that the body exerting the torque has an atmosphere model defined, and a shape model, if this is 
   not the case, the torque cannot be calculated. Furthermore, the body undergoing the torque needs to have aerodynamic coefficient interface defined, and needs to have its moment coefficients defined. The settings are created as follows, for an aerodynamic torque exerted by ``"Earth"`` on body ``"Apollo"``:

   .. code-block:: cpp

       SelectedTorqueMap selectedTorqueModelMap;
       selectedTorqueModelMap[ "Apollo" ][ "Earth" ].push_back( std::make_shared< TorqueSettings >( aerodynamic_torque ) );

   Requires the following environment models to be defined:

   - Atmosphere model for body exerting torque (set by :class:`AtmosphereSettings`).
   - Shape model for body exerting torque (set by :class:`BodyShapeSettings`).
   - Aerodynamic coefficient interface for body undergoing torque (set by :class:`AerodynamicCoefficientSettings`). NOTE: In the case that the aerodynamic coefficients are defined as a function of the vehicle orientation (e.g. angle of attack and sideslip angle), these angles can be manually or automatically defined. 
   - Current state of body undergoing torque and body with atmosphere.

   





