
.. role:: jsontype
.. role:: jsonkey
.. role:: arrow

- :jsontype:`string` :jsonkey:`type` (optional). Type of the translational propagator to be used. Possible values: :literal:`"cowell"`, :literal:`"encke"`, :literal:`"gaussKeplerian"`, :literal:`"gaussModifiedEquinoctial"`, :literal:`"unifiedStateModelQuaternions"`, :literal:`"unifiedStateModelModifiedRodriguesParameters"`, :literal:`"unifiedStateModelExponentialMap"`. Default value: :literal:`"cowell"`.
- :jsontype:`string[ ]` :jsonkey:`centralBodies` (mandatory). Names of the central bodies.
- :jsontype:`object[ ]` :jsonkey:`accelerations` (mandatory). Map in which each object contains a map of acceleration lists. The keys of the outer map are the names of the bodies undergoing the accelerations, while the keys of the inner maps are the names of the bodies exerting the accelerations. For instance, :jsonkey:`accelerations.satellite.Earth` is read as: accelerations on satellite caused by Earth.

	.. container:: toggle

		.. container:: header

			:arrow:`Acceleration`

		.. include:: Acceleration.rst
