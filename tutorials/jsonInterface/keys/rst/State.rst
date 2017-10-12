
.. role:: jsontype
.. role:: jsonkey
.. role:: arrow

- :jsontype:`string` :jsonkey:`type` (mandatory). Possible values: :literal:`"cartesian"`, :literal:`"keplerian"`, :literal:`"spherical"`.

.. container:: toggle

	.. container:: header

		:arrow:`Cartesian State`

	.. include:: CartesianState.rst

.. container:: toggle

	.. container:: header

		:arrow:`Keplerian State`

	.. include:: KeplerianState.rst

.. container:: toggle

	.. container:: header

		:arrow:`Spherical State`

	.. include:: SphericalState.rst
