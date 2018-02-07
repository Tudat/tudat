
.. role:: jsontype
.. role:: jsonkey
.. role:: arrow

- :jsontype:`number` :jsonkey:`initialStepSize` (mandatory). Initial step-size. Must be negative for backwards propagations.
- :jsontype:`number` :jsonkey:`minimumStepSize` (mandatory). An error will be thrown if the required step-size for the specified error tolerances is smaller than this value.
- :jsontype:`number` :jsonkey:`maximumStepSize` (mandatory). An error will be thrown if the required step-size for the specified error tolerances is larger than this value.
- :jsontype:`number` :jsonkey:`relativeErrorTolerance` (optional). Relative error tolerance for step-size control. Default value: :literal:`1.0E-12`.
- :jsontype:`number` :jsonkey:`absoluteErrorTolerance` (optional). Absolute error tolerance for step-size control. Default value: :literal:`1.0E-12`.
- :jsontype:`number` :jsonkey:`minimumOrder` (optional). Minimum order of numerical integrator. Default value: :literal:`6`.
- :jsontype:`number` :jsonkey:`maximumOrder` (optional). Maximum order of numerical integrator. Default value: :literal:`11`.
