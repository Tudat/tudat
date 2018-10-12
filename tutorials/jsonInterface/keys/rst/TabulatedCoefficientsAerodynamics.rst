
.. role:: jsontype
.. role:: jsonkey
.. role:: arrow

- :jsontype:`path[ 3 ]` :jsonkey:`forceCoefficients` (mandatory). Paths to the files containing the values for each force coefficient as a function of the independent variables.
- :jsontype:`path[ 3 ]` :jsonkey:`momentCoefficients` (optional). Paths to the files containing the values for each moment coefficient as a function of the independent variables.
- :jsontype:`number[ ]` :jsonkey:`independentVariableValues` (optional). Values of the independent variable, if not provided in the input files, when only one independent variable is used.
- :jsontype:`object` :jsonkey:`interpolator` (optional). Settings to be used for creating the one-dimensional interpolator of data, when only one independent variable is used.

	.. container:: toggle

		.. container:: header

			:arrow:`Interpolator`

		.. include:: Interpolator.rst
