
.. role:: jsontype
.. role:: jsonkey
.. role:: arrow

- :jsontype:`boolean` :jsonkey:`useDefaultSettings` (optional). Whether the default settings should be used. If any of the settings listed below is provided, the default settings will be overridden. Default value: :literal:`false`.
- :jsontype:`number[6]` :jsonkey:`initialState` (optional). Initial Cartesian state. Used to define a segment of the propagator's initial states (if not defined).
- :jsontype:`object` :jsonkey:`initialState` (optional). Initial state. Used to define a segment of the propagator's initial states (if not defined).

	.. container:: toggle

		.. container:: header

			:arrow:`State`

		.. include:: rst/State.rst
- :jsontype:`number` :jsonkey:`mass` (optional). Initial/constant mass. Used to define the constant mass of the body or its initial mass when it is varying, in which case it will be used to define an entry of the propagator's initial states (if not defined).
- :jsontype:`number[7]` :jsonkey:`rotationalState` (optional). Initial rotational state. Used to define a segment of the propagator's initial states (if not defined).
- :jsontype:`number` :jsonkey:`referenceArea` (optional). Area to be used for aerodynamics and radiation pressure settings (if not defined for these individual objects).
- :jsontype:`object` :jsonkey:`aerodynamics` (optional). Used to provide aerodynamics settings of the body.

	.. container:: toggle

		.. container:: header

			:arrow:`Aerodynamics`

		.. include:: rst/Aerodynamics.rst
- :jsontype:`object` :jsonkey:`atmosphere` (optional). Used to provide atmosphere settings of the body.

	.. container:: toggle

		.. container:: header

			:arrow:`Atmosphere`

		.. include:: rst/Atmosphere.rst
- :jsontype:`object` :jsonkey:`ephemeris` (optional). Used to provide ephemeris settings of the body.

	.. container:: toggle

		.. container:: header

			:arrow:`Ephemeris`

		.. include:: rst/Ephemeris.rst
- :jsontype:`object` :jsonkey:`gravityField` (optional). Used to provide gravity field settings of the body.

	.. container:: toggle

		.. container:: header

			:arrow:`Gravity Field`

		.. include:: rst/GravityField.rst
- :jsontype:`object[]` :jsonkey:`gravityFieldVariation` (optional). Used to provide a list of gravity field variation settings.

	.. container:: toggle

		.. container:: header

			:arrow:`Gravity Field Variation`

		.. include:: rst/GravityFieldVariation.rst
- :jsontype:`object` :jsonkey:`radiationPressure` (optional). Used to provide a map of radiation pressure settings, in which the keys are the names of the radiating bodies causing a radiation pressure on the body.

	.. container:: toggle

		.. container:: header

			:arrow:`Radiation Pressure`

		.. include:: rst/RadiationPressure.rst
- :jsontype:`object` :jsonkey:`rotationModel` (optional). Used to provide rotation model settings of the body.

	.. container:: toggle

		.. container:: header

			:arrow:`Rotation Model`

		.. include:: rst/RotationModel.rst
- :jsontype:`object` :jsonkey:`shapeModel` (optional). Used to provide shape model settings.

	.. container:: toggle

		.. container:: header

			:arrow:`Shape Model`

		.. include:: rst/ShapeModel.rst
