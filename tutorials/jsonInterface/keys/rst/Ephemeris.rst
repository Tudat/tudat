
.. role:: jsontype
.. role:: jsonkey
.. role:: arrow

- :jsontype:`string` :jsonkey:`type` (mandatory). Possible values: :literal:`"constant"`, :literal:`"kepler"`, :literal:`"approximatePlanetPositions"`, :literal:`"directSpice"`, :literal:`"interpolatedSpice"`, :literal:`"tabulated"`.
- :jsontype:`string` :jsonkey:`frameOrigin` (optional). Identifier of the reference frame origin. Default value: :literal:`"SSB"`.
- :jsontype:`string` :jsonkey:`frameOrigin` (optional). Identifier of the reference frame orientation. Default value: :literal:`"ECLIPJ2000"`.

.. container:: toggle

	.. container:: header

		:arrow:`Constant Ephemeris`

	.. include:: rst/ConstantEphemeris.rst

.. container:: toggle

	.. container:: header

		:arrow:`Kepler Ephemeris`

	.. include:: rst/KeplerEphemeris.rst

.. container:: toggle

	.. container:: header

		:arrow:`Approximate Planet Positions Ephemeris`

	.. include:: rst/ApproximatePlanetPositionsEphemeris.rst

.. container:: toggle

	.. container:: header

		:arrow:`Direct Spice Ephemeris`

	.. include:: rst/DirectSpiceEphemeris.rst

.. container:: toggle

	.. container:: header

		:arrow:`Interpolated Spice Ephemeris`

	.. include:: rst/InterpolatedSpiceEphemeris.rst

.. container:: toggle

	.. container:: header

		:arrow:`Tabulated Ephemeris`

	.. include:: rst/TabulatedEphemeris.rst
