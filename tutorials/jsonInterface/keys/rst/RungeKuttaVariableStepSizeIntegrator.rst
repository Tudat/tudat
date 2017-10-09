
.. role:: jsontype
.. role:: jsonkey
.. role:: arrow

- :jsontype:`string` :jsonkey:`rungeKuttaCoefficientSet` (mandatory). Identifier of the Runge-Kutta method to be used that will determine the coefficients set. Possible values: :literal:`"rungeKuttaFehlberg45"`, :literal:`"rungeKuttaFehlberg56"`, :literal:`"rungeKuttaFehlberg78"`, :literal:`"rungeKutta87DormandPrince"`.
- :jsontype:`number` :jsonkey:`initialStepSize` (mandatory). Initial step-size. Must be negative for backwards propagations.
- :jsontype:`number` :jsonkey:`minimumStepSize` (mandatory). An error will be thrown if the required step-size for the specified error tolerances is smaller than this value.
- :jsontype:`number` :jsonkey:`maximumStepSize` (mandatory). An error will be thrown if the required step-size for the specified error tolerances is larger than this value.
- :jsontype:`number` :jsonkey:`relativeErrorTolerance` (optional). Relative error tolerance for step-size control. Default value: :literal:`1.0E-12`.
- :jsontype:`number` :jsonkey:`absoluteErrorTolerance` (optional). Absolute error tolerance for step-size control. Default value: :literal:`1.0E-12`.
- :jsontype:`number` :jsonkey:`safetyFactorForNextStepSize` (optional). Safety factor for step-size control. Default value: :literal:`0.8`.
- :jsontype:`number` :jsonkey:`maximumFactorIncreaseForNextStepSize` (optional). Maximum increase factor in step-size in subsequent iterations. Default value: :literal:`4.0`.
- :jsontype:`number` :jsonkey:`minimumFactorDecreaseForNextStepSize` (optional). Minimum decrease factor in step-size in subsequent iterations. Default value: :literal:`0.1`.
