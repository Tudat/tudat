
.. role:: jsontype
.. role:: jsonkey
.. role:: arrow

- :jsontype:`number` :jsonkey:`order` (mandatory). Order of the Lagrange interpolator.
- :jsontype:`boolean` :jsonkey:`useLongDoubleTimeStep` (optional). Whether the time step is to be a long double, for higher precision. Default value: :literal:`false`.
- :jsontype:`string` :jsonkey:`lookupScheme` (optional). Selected type of lookup scheme for independent variables. Possible values: :literal:`"huntingAlgorithm"`, :literal:`"binarySearch"`. Default value: :literal:`"huntingAlgorithm"`.
- :jsontype:`string` :jsonkey:`lagrangeBoundaryHandling` (optional). Possible values: :literal:`"cubicSplineBoundary"`, :literal:`"noBoundary"`. Default value: :literal:`"cubicSplineBoundary"`.
