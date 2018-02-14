
.. role:: jsontype
.. role:: jsonkey
.. role:: arrow

- :jsontype:`number` :jsonkey:`initialStepSize` (mandatory). Initial step-size. Must be negative for backwards propagations.
- :jsontype:`number` :jsonkey:`maximumNumberOfSteps` (mandatory). Number of entries in the sequence, e.g. number of integrations used for a single extrapolation.
- :jsontype:`number` :jsonkey:`minimumStepSize` (mandatory). An error will be thrown if the required step-size for the specified error tolerances is smaller than this value.
- :jsontype:`number` :jsonkey:`maximumStepSize` (mandatory). An error will be thrown if the required step-size for the specified error tolerances is larger than this value.
- :jsontype:`number` :jsonkey:`relativeErrorTolerance` (optional). Relative error tolerance for step-size control. Default value: :literal:`1.0E-12`.
- :jsontype:`number` :jsonkey:`absoluteErrorTolerance` (optional). Absolute error tolerance for step-size control. Default value: :literal:`1.0E-12`.
- :jsontype:`number` :jsonkey:`safetyFactorForNextStepSize` (optional). Safety factor for step-size control. Default value: :literal:`0.7`.
- :jsontype:`number` :jsonkey:`maximumFactorIncreaseForNextStepSize` (optional). Maximum increase factor in step-size in subsequent iterations. Default value: :literal:`10.0`.
- :jsontype:`number` :jsonkey:`minimumFactorDecreaseForNextStepSize` (optional). Minimum decrease factor in step-size in subsequent iterations. Default value: :literal:`0.1`.
