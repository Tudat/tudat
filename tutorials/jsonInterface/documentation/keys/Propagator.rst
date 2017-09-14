.. role:: arrow

- :literal:`string[]` :class:`bodiesToPropagate` (mandatory) Names of the bodies to propagate.
- :literal:`numeric[]` :class:`initialStates` (optional) Vector containing the concatenated initial states of all the bodies to be propagated. If not defined, it will be determined from the properties :literal:`initialState`, :literal:`mass` or :literal:`rotationalState` of the bodies to be propagated. If all the bodies to propagate are celestial bodies and this property is not defined, the ephemeris at the initial epoch will be used.
- :literal:`string` :class:`integratedStateType` (optional). Possible values: :literal:`"translational"`, :literal:`"mass"`, :literal:`"rotational"`. Default value: :literal:`"translational"`.

	.. container:: toggle

		.. container:: header

			:arrow:`Translational Propagator`

		.. include:: keys/TranslationalPropagator.rst
