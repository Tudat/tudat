.. _tudatFeaturesFrameworkTorqueModelSetup:    

Rotational Model Set-up
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

.. method:: aerodynamic torque

   Torque exerted by a body with an atmosphere model and shape model on another body. It is important that the body exerting the torque has an atmosphere model defined, and a shape model, if this is 
   not the case, the torque cannot be calculated. Furthermore, the body undergoing the torque needs to have aerodynamic coefficient interface defined, and needs to have its moment coefficients defined.

.. method:: second-order gravitational torque

   Torque exerted by a point mass on a body with a degree two spherical harmonics mass distribution. A degree two spherical harmonics mass distribution can be represented by an inertia tensor, thus 
   for this torque model, the body undergoing the torque needs to have an inertia tensor defined. The body exerting the torque only needs to have a gravitational model defined. 

   





