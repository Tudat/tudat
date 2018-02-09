
.. role:: jsontype
.. role:: jsonkey
.. role:: arrow

- :jsontype:`number` :jsonkey:`initialTime` (mandatory). Initial time from which interpolated data from Spice should be created.
- :jsontype:`number` :jsonkey:`finalTime` (mandatory). Final time until which interpolated data from Spice should be created.
- :jsontype:`number` :jsonkey:`timeStep` (mandatory). Time step with which interpolated data from Spice should be created.
- :jsontype:`object` :jsonkey:`interpolator` (optional). Settings to be used for the state interpolation. Default value: :literal:`{"type":"lagrange","order":6}`.

	.. container:: toggle

		.. container:: header

			:arrow:`Interpolator`

		.. include:: Interpolator.rst
- :jsontype:`boolean` :jsonkey:`useLongDoubleStates` (optional). Whether to use long double states doubles for higher precision. Default value: :literal:`false`.
