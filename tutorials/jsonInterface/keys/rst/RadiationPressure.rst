
.. role:: jsontype
.. role:: jsonkey
.. role:: arrow

- :jsontype:`string` :jsonkey:`type` (optional). Possible values: :literal:`"cannonBall"`. Default value: :literal:`"cannonBall"`.
- :jsontype:`string` :jsonkey:`sourceBody` (optional). Name of the radiating body. If not provided, it will be retrieved from the key to which the radiation pressure object is assigned.
- :jsontype:`string[ ]` :jsonkey:`occultingBodies` (optional). List with the names of the bodies causing (partial) occultation, if any.

.. container:: toggle

	.. container:: header

		:arrow:`Cannon Ball Radiation Pressure`

	.. include:: CannonBallRadiationPressure.rst
