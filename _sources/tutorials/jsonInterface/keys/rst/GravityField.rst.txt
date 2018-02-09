
.. role:: jsontype
.. role:: jsonkey
.. role:: arrow

- :jsontype:`string` :jsonkey:`type` (mandatory). Possible values: :literal:`"pointMass"`, :literal:`"pointMassSpice"`, :literal:`"sphericalHarmonic"`.

.. container:: toggle

	.. container:: header

		:arrow:`Point Mass Gravity Field`

	.. include:: PointMassGravityField.rst

.. container:: toggle

	.. container:: header

		:arrow:`Spherical Harmonic Gravity Field from built-in model`

	.. include:: SphericalHarmonicGravityFieldfrombuiltinmodel.rst

.. container:: toggle

	.. container:: header

		:arrow:`Spherical Harmonic Gravity Field from custom file`

	.. include:: SphericalHarmonicGravityFieldfromcustomfile.rst

.. container:: toggle

	.. container:: header

		:arrow:`Spherical Harmonic Gravity Field with manual coefficients`

	.. include:: SphericalHarmonicGravityFieldwithmanualcoefficients.rst
