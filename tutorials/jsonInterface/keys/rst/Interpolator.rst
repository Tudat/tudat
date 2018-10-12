
.. role:: jsontype
.. role:: jsonkey
.. role:: arrow

- :jsontype:`string` :jsonkey:`type` (mandatory). Possible values: :literal:`"linear"`, :literal:`"cubicSpline"`, :literal:`"lagrange"`, :literal:`"hermiteSpline"`, :literal:`"piecewiseConstant"`.
- :jsontype:`string` :jsonkey:`lookupScheme` (optional). Selected type of lookup scheme for independent variables. Possible values: :literal:`"huntingAlgorithm"`, :literal:`"binarySearch"`. Default value: :literal:`"huntingAlgorithm"`.
- :jsontype:`boolean` :jsonkey:`useLongDoubleTimeStep` (optional). Whether the time step is to be a long double, for higher precision. Default value: :literal:`false`.
- :jsontype:`string` :jsonkey:`boundaryHandling` (optional). Boundary handling method, in case the independent variable is outside the specified range. Possible values: :literal:`"throwException"`, :literal:`"boundaryValue"`, :literal:`"boundaryValueWithWarning"`, :literal:`"extrapolate"`, :literal:`"extrapolateWithWarning"`. Default value: :literal:`"extrapolate"`.

.. container:: toggle

	.. container:: header

		:arrow:`LagrangeInterpolator`

	.. include:: LagrangeInterpolator.rst
