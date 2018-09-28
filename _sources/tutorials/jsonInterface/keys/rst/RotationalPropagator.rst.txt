
.. role:: jsontype
.. role:: jsonkey
.. role:: arrow

- :jsontype:`string` :jsonkey:`type` (optional). Type of the translational propagator to be used. Possible values: :literal:`"quaternions"`, :literal:`"modifiedRodriguesParameters"`, :literal:`"exponentialMap"`. Default value: :literal:`"quaternions"`.
- :jsontype:`object[]` :jsonkey:`torques` (mandatory). Map in which each object contains a map of torque lists. The keys of the outer map are the names of the bodies undergoing the torque, while the keys of the inner maps are the names of the bodies exerting the torque. For instance, :jsonkey:`torques.Earth.Moon` is read as: torque on Earth caused by the Moon.

	.. container:: toggle

		.. container:: header

			:arrow:`Torque`

		.. include:: Torque.rst
