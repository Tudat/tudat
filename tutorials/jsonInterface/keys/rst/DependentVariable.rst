
.. role:: jsontype
.. role:: jsonkey
.. role:: arrow

- :jsontype:`string` :jsonkey:`dependentVariableType` (mandatory). Possible values: :literal:`"machNumber"`, :literal:`"altitude"`, :literal:`"airspeed"`, :literal:`"localDensity"`, :literal:`"localTemperature"`, :literal:`"stagnationPointHeatFlux"`, :literal:`"localDynamicPressure"`, :literal:`"localAerodynamicHeatRate"`, :literal:`"relativeSpeed"`, :literal:`"relativePosition"`, :literal:`"relativeDistance"`, :literal:`"relativeVelocity"`, :literal:`"radiationPressure"`, :literal:`"totalAccelerationNorm"`, :literal:`"accelerationNorm"`, :literal:`"totalAcceleration"`, :literal:`"acceleration"`, :literal:`"sphericalHarmonicsAccelerationTerms"`, :literal:`"aerodynamicForceCoefficients"`, :literal:`"aerodynamicMomentCoefficients"`, :literal:`"rotationMatrixToBodyFixedFrame"`, :literal:`"intermediateAerodynamicRotationMatrix"`, :literal:`"relativeBodyAerodynamicOrientationAngle"`, :literal:`"bodyFixedAirspeedBasedVelocity"`, :literal:`"totalAerodynamicGLoad"`, :literal:`"geodeticLatitude"`, :literal:`"controlSurfaceDeflection"`, :literal:`"totalMassRates"`, :literal:`"lvlhToInertialFrameRotation"`, :literal:`"periapsisAltitude"`, :literal:`"totalTorqueNorm"`, :literal:`"torqueNorm"`, :literal:`"totalTorque"`, :literal:`"torque"`, :literal:`"bodyFixedGroundspeedBasedVelocity"`, :literal:`"bodyFixedRelativeCartesianPosition"`, :literal:`"bodyFixedRelativeSphericalPosition"`, :literal:`"keplerElements"`, :literal:`"modifiedEquinoctialElements"`.
- :jsontype:`string` :jsonkey:`body` (mandatory). Name of the associated body. For some variables, name of the body undergoing the acceleration or torque.
- :jsontype:`string` :jsonkey:`relativeToBody` (optional). For some variables, name of the body with respect to which the variable value should be computed (e.g. relative position or altitude w.r.t. Earth).
- :jsontype:`number` :jsonkey:`componentIndex` (optional). For vectorial variables, the index (starting from 0) to be accessed. If not specified, the whole vector is used.

.. container:: toggle

	.. container:: header

		:arrow:`Acceleration Dependent Variable`

	.. include:: AccelerationDependentVariable.rst

.. container:: toggle

	.. container:: header

		:arrow:`Spherical Harmonic Acceleration Terms Dependent Variable`

	.. include:: SphericalHarmonicAccelerationTermsVariable.rst

.. container:: toggle

	.. container:: header

		:arrow:`Variation of Spherical Harmonic Acceleration Dependent Variable`

	.. include:: SingleVariationSphericalHarmonicAcceleration.rst

.. container:: toggle

	.. container:: header

		:arrow:`Variation of Spherical Harmonic Acceleration Terms Dependent Variable`

	.. include:: SingleVariationSingleTermSphericalHarmonicAcceleration.rst

.. container:: toggle

	.. container:: header

		:arrow:`Torque Dependent Variable`

	.. include:: TorqueDependentVariable.rst

.. container:: toggle

	.. container:: header

		:arrow:`Intermediate Aerodynamic Rotation Matrix Dependent Variable`

	.. include:: IntermediateAerodynamicRotationMatrixDependentVariable.rst

.. container:: toggle

	.. container:: header

		:arrow:`Relative Body Aerodynamic Orientation Angle Variable`

	.. include:: RelativeBodyAerodynamicOrientationAngleVariable.rst
