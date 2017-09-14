.. role:: arrow

- :literal:`string` :class:`integratedStateType` (optional). Possible values: :literal:`"translational"`, :literal:`"mass"`, :literal:`"rotational"`. Default value: :literal:`"translational"`.
- :literal:`string[]` :class:`bodiesToPropagate` (mandatory) Names of the bodies to propagate.
- :literal:`numeric[]` :class:`initialStates` (optional) Vector containing the concatenated initial states of all the bodies to be propagated. If not defined, it will be determined from the properties :literal:`initialState`, :literal:`mass` or :literal:`rotationalState` of the bodies to be propagated. If all the bodies to propagate are celestial bodies and this property is not defined, the ephemeris at the initial epoch will be used.

.. container:: toggle

	.. container:: header

		:arrow:`Translational Propagator`

	.. include:: keys/TranslationalPropagator.rst

.. container:: toggle

	.. container:: header

		:arrow:`Mass Propagator`

	.. include:: keys/MassPropagator.rst

.. container:: toggle

	.. container:: header

		:arrow:`Rotational Propagator`

	.. include:: keys/RotationalPropagator.rst
