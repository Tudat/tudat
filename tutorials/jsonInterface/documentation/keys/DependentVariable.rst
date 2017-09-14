.. role:: arrow

- :literal:`string` :class:`dependentVariableType` (mandatory). Possible values: :literal:`"machNumber"`, :literal:`"altitude"`, :literal:`"airspeed"`, :literal:`"localDensity"`, :literal:`"relativeSpeed"`, :literal:`"relativePosition"`, :literal:`"relativeDistance"`, :literal:`"relativeVelocity"`, :literal:`"radiationPressure"`, :literal:`"totalAccelerationNorm"`, :literal:`"accelerationNorm"`, :literal:`"totalAcceleration"`, :literal:`"acceleration"`, :literal:`"aerodynamicForceCoefficients"`, :literal:`"aerodynamicMomentCoefficients"`, :literal:`"rotationMatrixToBodyFixedFrame"`, :literal:`"intermediateAerodynamicRotationMatrix"`, :literal:`"relativeBodyAerodynamicOrientationAngle"`, :literal:`"bodyFixedAirspeedBasedVelocity"`, :literal:`"totalAerodynamicGLoad"`, :literal:`"stagnationPointHeatFlux"`, :literal:`"localTemperature"`, :literal:`"geodeticLatitude"`, :literal:`"controlSurfaceDeflection"`, :literal:`"totalMassRates"`, :literal:`"lvlhToInertialFrameRotation"`, :literal:`"periapsisAltitude"`, :literal:`"totalTorqueNorm"`, :literal:`"torqueNorm"`, :literal:`"totalTorque"`, :literal:`"torque"`, :literal:`"bodyFixedGroundspeedBasedVelocity"`.
- :literal:`string` :class:`body` (mandatory) Name of the associated body. For some variables, name of the body undergoing the acceleration or torque.
- :literal:`string` :class:`relativeToBody` (optional) For some variables, name of the body with respect to which the variable value should be computed (e.g. relative position or altitude w.r.t. Earth).
- :literal:`int` :class:`componentIndex` (optional) For vectorial variables, the index (starting from 0) to be accessed. If not specified, the whole vector is used.

.. container:: toggle

	.. container:: header

		:arrow:`Acceleration Dependent Variable`

	.. include:: keys/AccelerationDependentVariable.rst

.. container:: toggle

	.. container:: header

		:arrow:`Torque Dependent Variable`

	.. include:: keys/TorqueDependentVariable.rst

.. container:: toggle

	.. container:: header

		:arrow:`Intermediate Aerodynamic Rotation Matrix Dependent Variable`

	.. include:: keys/IntermediateAerodynamicRotationMatrixDependentVariable.rst

.. container:: toggle

	.. container:: header

		:arrow:`Relative Body Aerodynamic Orientation Angle Variable`

	.. include:: keys/RelativeBodyAerodynamicOrientationAngleVariable.rst
