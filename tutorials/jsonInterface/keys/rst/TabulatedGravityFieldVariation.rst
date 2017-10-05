
.. role:: jsontype
.. role:: jsonkey
.. role:: arrow

- :jsontype:`object` :jsonkey:`cosineCoefficientCorrections` (mandatory). Map of corrections to cosine coefficients, where the keys correspond to the epoch.
- :jsontype:`object` :jsonkey:`sineCoefficientCorrections` (mandatory). Map of corrections to sine coefficients, where the keys correspond to the epoch.
- :jsontype:`number` :jsonkey:`minimumDegree` (mandatory). Minimum degree of spherical harmonic corrections.
- :jsontype:`number` :jsonkey:`minimumOrder` (mandatory). Minimum order of spherical harmonic corrections.
- :jsontype:`object` :jsonkey:`interpolator` (mandatory). Settings that are to be used to create an interpolator for the gravity field variations immediately upon creation.

	.. container:: toggle

		.. container:: header

			:arrow:`Interpolator`

		.. include:: Interpolator.rst
