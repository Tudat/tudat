
.. role:: jsontype
.. role:: jsonkey
.. role:: arrow

- :jsontype:`number` :jsonkey:`initialEpoch` (optional). Used for the initial time of the integrator (if not provided inside the integrator object) and for interpolation of Spice ephemeris (if Spice enabled and set to preload ephemeris).
- :jsontype:`number` :jsonkey:`finalEpoch` (optional). Used for interpolation of Spice ephemeris (if Spice enabled and set to preload ephemeris) and to create the time-termination condition. If smaller than :jsonkey:`initialEpoch`, the dynamics will be propagated backwards in time.
- :jsontype:`string` :jsonkey:`globalFrameOrigin` (optional). Used to set the global frame origin for the propagation. Default value: :literal:`"SSB"`.
- :jsontype:`string` :jsonkey:`globalFrameOrientation` (optional). Used to set the global frame orientation for the propagation. Updates the frame orientation of epehemeris and rotation model settings if not specified. Default value: :literal:`"ECLIPJ2000"`.
- :jsontype:`object` :jsonkey:`spice` (optional). Used to provide Spice settings.

	.. container:: toggle

		.. container:: header

			:arrow:`Spice`

		.. include:: rst/Spice.rst
- :jsontype:`object` :jsonkey:`bodies` (mandatory). Used to define all the bodies to be considered in the propagation. The keys of the object are the body names.

	.. container:: toggle

		.. container:: header

			:arrow:`Body`

		.. include:: rst/Body.rst
- :jsontype:`object[]` :jsonkey:`propagator` (mandatory). Used to define the propagator(s).

	.. container:: toggle

		.. container:: header

			:arrow:`Propagator`

		.. include:: rst/Propagator.rst
- :jsontype:`object` :jsonkey:`integrator` (mandatory). Used to define the integrator.

	.. container:: toggle

		.. container:: header

			:arrow:`Integrator`

		.. include:: rst/Integrator.rst
- :jsontype:`object` :jsonkey:`termination` (optional). Used to define termination condition(s). The time termination condition will be created automatically and added to the provided conditions (if any) if the key :jsonkey:`.finalEpoch` is specified. Either a single termination condition object or a multiple termination condition object can be provided.

	.. container:: toggle

		.. container:: header

			:arrow:`Single Termination Condition`

		.. include:: rst/SingleTerminationCondition.rst

	.. container:: toggle

		.. container:: header

			:arrow:`Multiple Termination Condition`

		.. include:: rst/MultipleTerminationCondition.rst
- :jsontype:`object[]` :jsonkey:`export` (optional). Used to define the export settings. Each element represents an output file to which results will be saved.

	.. container:: toggle

		.. container:: header

			:arrow:`Export`

		.. include:: rst/Export.rst
- :jsontype:`object` :jsonkey:`options` (optional). Used to configure options for the application.

	.. container:: toggle

		.. container:: header

			:arrow:`Options`

		.. include:: rst/Options.rst
