
.. role:: jsontype
.. role:: jsonkey
.. role:: arrow

- :jsontype:`number` :jsonkey:`initialTime` (mandatory). Start time for model interpolation.
- :jsontype:`number` :jsonkey:`finalTime` (mandatory). End time for model interpolation.
- :jsontype:`number` :jsonkey:`timeStep` (mandatory). Time step with which to evaluate model, and provide input to interpolator.
- :jsontype:`object` :jsonkey:`interpolator` (mandatory). The settings for the interpolator to be used in the model interpolation.

	.. container:: toggle

		.. container:: header

			:arrow:`Interpolator`

		.. include:: Interpolator.rst
