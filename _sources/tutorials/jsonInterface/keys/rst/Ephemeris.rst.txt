
.. role:: jsontype
.. role:: jsonkey
.. role:: arrow

- :jsontype:`string` :jsonkey:`type` (mandatory). Possible values: :literal:`"constant"`, :literal:`"kepler"`, :literal:`"approximatePlanetPositions"`, :literal:`"directSpice"`, :literal:`"interpolatedSpice"`, :literal:`"tabulated"`.
- :jsontype:`string` :jsonkey:`frameOrigin` (optional). Identifier of the reference frame origin. Default value: :literal:`"SSB"`.
- :jsontype:`string` :jsonkey:`frameOrigin` (optional). Identifier of the reference frame orientation. Default value: :literal:`"ECLIPJ2000"`.

.. container:: toggle

	.. container:: header

		:arrow:`Constant Ephemeris`

	.. include:: ConstantEphemeris.rst

.. container:: toggle

	.. container:: header

		:arrow:`Kepler Ephemeris`

	.. include:: KeplerEphemeris.rst

.. container:: toggle

	.. container:: header

		:arrow:`Approximate Planet Positions Ephemeris`

	.. include:: ApproximatePlanetPositionsEphemeris.rst

.. container:: toggle

	.. container:: header

		:arrow:`Direct Spice Ephemeris`

	.. include:: DirectSpiceEphemeris.rst

.. container:: toggle

	.. container:: header

		:arrow:`Interpolated Spice Ephemeris`

	.. include:: InterpolatedSpiceEphemeris.rst

.. container:: toggle

	.. container:: header

		:arrow:`Tabulated Ephemeris`

	.. include:: TabulatedEphemeris.rst
