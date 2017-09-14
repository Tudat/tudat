.. role:: arrow

- :literal:`numeric` :class:`initialEpoch` (optional) Used for the initial time of the integrator (if not provided inside the integrator object) and for interpolation of Spice ephemeris (if Spice enabled and set to preload ephemeris).
- :literal:`numeric` :class:`finalEpoch` (optional) Used for interpolation of Spice ephemeris (if Spice enabled and set to preload ephemeris) and to create the time-termination condition. If smaller than :literal:`initialEpoch`, the dynamics will be propagated backwards in time.
- :literal:`string` :class:`globalFrameOrigin` (optional) Used to set the global frame origin for the propagation. Default value: :literal:`"SSB"`.
- :literal:`string` :class:`globalFrameOrientation` (optional) Used to set the global frame orientation for the propagation. Updates the frame orientation of epehemeris and rotation model settings if not specified. Default value: :literal:`"ECLIPJ2000"`.
- :literal:`object` :class:`spice` (optional) Used to provide Spice settings.

	.. container:: toggle

		.. container:: header

			:arrow:`Spice`

		.. include:: keys/Spice.rst
- :literal:`object` :class:`bodies` (mandatory) Used to define all the bodies to be considered in the propagation. The keys of the object are the body names.

	.. container:: toggle

		.. container:: header

			:arrow:`Body`

		.. include:: keys/Body.rst
- :literal:`object[]` :class:`propagator` (mandatory) Used to define the propagator(s).

	.. container:: toggle

		.. container:: header

			:arrow:`Propagator`

		.. include:: keys/Propagator.rst
- :literal:`object` :class:`integrator` (mandatory) Used to define the integrator.

	.. container:: toggle

		.. container:: header

			:arrow:`Integrator`

		.. include:: keys/Integrator.rst
- :literal:`object` :class:`termination` (optional) Used to define termination condition(s). The time termination condition will be created automatically and added to the provided conditions (if any) if the key :literal:`~/finalEpoch` is specified.

	.. container:: toggle

		.. container:: header

			:arrow:`Single Termination Condition`

		.. include:: keys/SingleTerminationCondition.rst

	.. container:: toggle

		.. container:: header

			:arrow:`Multiple Termination Condition`

		.. include:: keys/MultipleTerminationCondition.rst
- :literal:`object[]` :class:`export` (optional) Used to define the export settings. Each element represents an output file to which results will be saved.

	.. container:: toggle

		.. container:: header

			:arrow:`Export`

		.. include:: keys/Export.rst
- :literal:`object` :class:`options` (optional) Used to configure options for the application.

	.. container:: toggle

		.. container:: header

			:arrow:`Options`

		.. include:: keys/Options.rst
