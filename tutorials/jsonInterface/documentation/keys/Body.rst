.. role:: arrow

- :literal:`bool` :class:`useDefaultSettings` (optional) Whether the default settings should be used. If any of the settings listed below is provided, the default settings will be overridden. Default value: :literal:`false`.
- :literal:`numeric[6]` :class:`initialState` (optional) Initial Cartesian state. Used to define a segment of the propagator's initial states (if not defined).
- :literal:`object` :class:`initialState` (optional) Initial state. Used to define a segment of the propagator's initial states (if not defined).

	.. container:: toggle

		.. container:: header

			:arrow:`State`

		.. include:: keys/State.rst
- :literal:`numeric` :class:`mass` (optional) Initial/constant mass. Used to define the constant mass of the body or its initial mass when it is varying, in which case it will be used to define an entry of the propagator's initial states (if not defined).
- :literal:`numeric[7]` :class:`rotationalState` (optional) Initial rotational state. Used to define a segment of the propagator's initial states (if not defined).
- :literal:`numeric` :class:`referenceArea` (optional) Area to be used for aerodynamics and radiation pressure settings (if not defined for these individual objects).
- :literal:`object` :class:`aerodynamics` (optional) Used to provide aerodynamics settings of the body.
- :literal:`object` :class:`atmosphere` (optional) Used to provide atmosphere settings of the body.
- :literal:`object` :class:`ephemeris` (optional) Used to provide ephemeris settings of the body.
- :literal:`object` :class:`gravityField` (optional) Used to provide gravity field settings of the body.
- :literal:`object[]` :class:`gravityFieldVariation` (optional) Used to provide a list of gravity field variation settings.
- :literal:`object` :class:`radiationPressure` (optional) Used to provide a map of radiation pressure settings, in which the keys are the name of the radiating bodies acting on the body.
- :literal:`object` :class:`rotationModel` (optional) Used to provide rotation model settings of the body.
- :literal:`object` :class:`shapeModel` (optional) Used to provide shape model settings.
