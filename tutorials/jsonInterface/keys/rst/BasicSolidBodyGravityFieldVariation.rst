
.. role:: jsontype
.. role:: jsonkey
.. role:: arrow

- :jsontype:`string[ ]` :jsonkey:`deformingBodies` (mandatory). List with the names of the bodies causing tidal deformation.
- :jsontype:`string[ ][ ]` :jsonkey:`loveNumbers` (mandatory). List of Love numbers for the deformed body. Each element is a complex number, which in JSON is represented as a string. For instance, :literal:`2 - 0.5i` is written as :literal:`"(2,-0.5)"`.
- :jsontype:`number` :jsonkey:`referenceRadius` (mandatory). Reference (typically equatorial) radius of the body being deformed.
- :jsontype:`object` :jsonkey:`modelInterpolation` (optional). Settings that are to be used to create an interpolator for the gravity field variations immediately upon creation. If not defined, no interpolation is used, and the model is evaluated during propagation.

	.. container:: toggle

		.. container:: header

			:arrow:`Model Interpolation`

		.. include:: ModelInterpolation.rst
