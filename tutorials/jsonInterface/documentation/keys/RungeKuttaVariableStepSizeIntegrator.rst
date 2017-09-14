.. role:: arrow

- :literal:`string` :class:`rungeKuttaCoefficientSet` (mandatory) Identifier of the Runge-Kutta method to be used that will determine the coefficients set. Possible values: :literal:`"rungeKuttaFehlberg45"`, :literal:`"rungeKuttaFehlberg56"`, :literal:`"rungeKuttaFehlberg78"`, :literal:`"rungeKutta87DormandPrince"`.
- :literal:`numeric` :class:`initialStepSize` (mandatory) Initial step-size. Must be negative for backwards propagations.
- :literal:`numeric` :class:`minimumStepSize` (mandatory) An error will be thrown if the required step-size for the specified error tolerances is smaller than this value.
- :literal:`numeric` :class:`maximumStepSize` (mandatory) An error will be thrown if the required step-size for the specified error tolerances is larger than this value.
- :literal:`numeric` :class:`relativeErrorTolerance` (optional) Relative error tolerance for step-size control. Default value: :literal:`1.0E-12`.
- :literal:`numeric` :class:`absoluteErrorTolerance` (optional) Absolute error tolerance for step-size control. Default value: :literal:`1.0E-12`.
- :literal:`numeric` :class:`safetyFactorForNextStepSize` (optional) Safety factor for step-size control. Default value: :literal:`0.8`.
- :literal:`numeric` :class:`maximumFactorIncreaseForNextStepSize` (optional) Maximum increase factor in step-size in subsequent iterations. Default value: :literal:`4.0`.
- :literal:`numeric` :class:`minimumFactorDecreaseForNextStepSize` (optional) Minimum decrease factor in step-size in subsequent iterations. Default value: :literal:`0.1`.
