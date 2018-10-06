.. _tudatFeaturesPropagatorSettingsCoordinates:

Propagator Settings: Conventional vs. Propagated Coordinates
============================================================

In the description of some of the objects in this part of the wiki about the simulation set-up, you may have noticed the use of two names to describe the states of an object. These two names are *conventional* and *propagated*, and they describe two slightly different concepts. In this part of the Wiki, you will get to know in what they differ and how that may affect your application. 

.. - input coordinates is the conventional type
.. - acceleration model and default body description states are conventional type
.. - equations of motion are propagation type
.. - integration result is propagation type
.. - flow diagram?

.. method:: Conventional Coordinates

   These coordinates are mainly used to describe the state of the object(s) you are propagating outsite of integration. This means that you will find *conventional* coordinates in these scenarios:

      - To describe the initial conditions of an object
      - To update the acceleration model of an object
      - To update the environment model of an object; this also means that the states extracted from the :class:`Body` are expressed in the conventional coordinates
      - As an output to the :literal:`getEquationsOfMotionNumericalSolution` function of the :class:`DynamicsSimulator` object

   For sake of completeness, below you will find the conventional coordinates for each propagator type which supports a multiple propagated coordinates:

      - **Translational Motion**: Cartesian coordinates
      - **Rotational Motion**: Quaternions

.. method:: Propagated Coordinates

   The *propagated* coordinates, on the other hand, are mainly used to describe the state of an object inside of the integration environment. Thus, you can expect to see these elements here:

      - To describe the equations of motion
      - To describe the state and state derivative during integration
      - As an output to the :literal:`getEquationsOfMotionNumericalSolutionRaw` function of the :class:`DynamicsSimulator` object

As a user, you will generally only interact with the conventional coordinates, but you will have the choice over which propagated coordinate to use for propagation/integration. Note that setting a type of coordinates for propagation does not necessarily also constrain which propagator is used. For instance, Cowell and Encke propagation, both use Cartesian coordinates as the propagated states. The same is also true vice-versa: Gauss' version of Lagrange's planetary equations can be used with both Keplerian and modified equinoctial elements. 

In the remaining part of this page, you will find a description of where you need to keep a closer look at the difference between these two types of coordinates and which coordinates are available in Tudat.

Important To Keep In Mind
~~~~~~~~~~~~~~~~~~~~~~~~~

Since the conventional coordinates are used to update the environment and accelerations of the bodies, but the propagated coordinates are the ones used in propagation, you can see that whenever the conventional and propagated coordinates differ, there is a need to convert between the two at every time step (or even multiple times, if the time step is divided in multiple steps for integration). Thus, this leads to a set of extra stages to be perfomed during propagation, which may in turn lead to a longer computation time. This convertion is also necessary when outputting the conventional state at the end of propagation. 

.. note:: The fact that using a different set of propagated coordinates may lead to a longer computation time is not always true. As a matter of fact, the default translational propagator (i.e., :literal:`cowell`) is considerably slower and less accurate than other propagators, in certain situations. Check out :ref:`walkthroughsPropagatorTypesComparison` to get an idea of the difference in performance among the various translational propagators offered by Tudat.

Another fact to consider, is that sometimes there may be a difference between the size of the conventional and propagates states. For instance, a Cartesian state is expressed with 6 elements, but the USM7 state with 7. This may lead to some confusion when extracting the results, so keep this in mind. In the next section, you can find the size of each propagated type used in Tudat.

.. _tudatFeaturesPropagatorSettingsPropagatedCoordinates:

Propagated Coordinates in Tudat
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For the following state types, one can choose between the listed propagated coordinates. You can also find a similar list in :ref:`tudatFeaturesPropagatorSettings`, where each propagator available is listed.

.. method:: Translational Motion

   - Cartesian coordinates (with :literal:`cowell` and :literal:`encke` propagators); state size: 6
   - Keplerian elements (with :literal:`gauss_keplerian` propagator); state size: 6
   - Modified equinoctial elements (with :literal:`gauss_modified_equinoctial` propagator); state size: 6
   - Unified state model with quaternions, or USM7 (with :literal:`unified_state_model_quaternions` propagator); state size: 7
   - Unified state model with modified Rodrigues parameters, or USM6 (with :literal:`unified_state_model_modified_rodrigues_parameters` propagator); state size: 7
   - Unified state model with exponential map, or USMEM (with :literal:`unified_state_model_exponential_map` propagator); state size: 7

.. method:: Rotational Motion

   - Quaternions (with :literal:`quaternions` propagator); state size: 7
   - Modified Rodrigues parameters, or MRPs (with :literal:`modified_rodrigues_parameters` propagator); state size: 7
   - Exponential map (with :literal:`exponential_map` propagator); state size: 7

Any other state type that has not been mentioned, is only described by one coordinate type, and thus their conventional and propagated coordinates always match.

.. - note difference in size between some coordinates
.. - note transformation between conventional and propagated coordinates