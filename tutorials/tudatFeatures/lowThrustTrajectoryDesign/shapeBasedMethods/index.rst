.. _tudatFeaturesShapeBasedMethodsIndex:

Shape-Based Methods 
===================

Shape-based methods are used to provide a good preliminary design for a given trajectory, in an efficient way. To do so, they rely on analytical shaping of the trajectory: they assume the trajectory respects a certain shape (which can be fully described analytically), and the parameters defining this shape are then computed to ensure the trajectory satisfies some boundary conditions (states at departure and arrival, time of flight, ...). The analytical formulation reduces the computational load significantly, which is a major advantages of those shaping methods over the traditional direct and indirect methods.  

All shape-based methods derive from the class :literal:`ShapeBasedMethodLeg`, described below.

.. class:: ShapeBasedMethodLeg

This is the base class used to derive the different shaping methods. It itself inherits from the base class LowThrustLeg. Each class inherited from :literal:`ShapeBasedMethodLeg` contains the following methods:

	- :literal:`getInitialValueInpendentVariable`
		Returns the initial value of the independent variable with respect to which the trajectory is calculated. Usually, time is the independent variable but several shaping methods use polar or azimuth angles.

	- :literal:`getFinalValueInpendentVariable` 
		Returns the final value of the independent variable.

Several shape-based methods have been implemented and inherit from this :literal:`ShapeBasedMethodLeg` base class.

.. toctree::
   :numbered:
   :maxdepth: 2

   hodographicShaping
   sphericalShaping
   	

