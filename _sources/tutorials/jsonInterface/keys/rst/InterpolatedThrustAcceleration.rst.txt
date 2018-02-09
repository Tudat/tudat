
.. role:: jsontype
.. role:: jsonkey
.. role:: arrow

- :jsontype:`object` :jsonkey:`dataInterpolation` (mandatory). Settings to create the interpolator interface.

	.. container:: toggle

		.. container:: header

			:arrow:`Data Interpolation`

		.. include:: DataInterpolation.rst
- :jsontype:`number` :jsonkey:`specificImpulse` (mandatory). Constant specific impulse.
- :jsontype:`string` :jsonkey:`frame` (mandatory). Identifier of frame in which thrust returned by interpolator is expressed. Possible values: :literal:`"intertial"`, :literal:`"lvlh"`.
- :jsontype:`string` :jsonkey:`centralBody` (mandatory). Name of the central body.
