
.. role:: jsontype
.. role:: jsonkey
.. role:: arrow

- :jsontype:`string` :jsonkey:`integratedStateType` (optional). Possible values: :literal:`"translational"`, :literal:`"mass"`, :literal:`"rotational"`. Default value: :literal:`"translational"`.
- :jsontype:`string[ ]` :jsonkey:`bodiesToPropagate` (mandatory). Names of the bodies to propagate.
- :jsontype:`number[ ]` :jsonkey:`initialStates` (optional). Vector containing the concatenated initial states of all the bodies to be propagated. If not defined, it will be determined from the properties :jsonkey:`initialState`, :jsonkey:`mass` or :jsonkey:`rotationalState` of the bodies to be propagated. If all the bodies to propagate are celestial bodies and this property is not defined, the ephemeris at the initial epoch will be used.

.. container:: toggle

	.. container:: header

		:arrow:`Translational Propagator`

	.. include:: TranslationalPropagator.rst

.. container:: toggle

	.. container:: header

		:arrow:`Mass Propagator`

	.. include:: MassPropagator.rst

.. container:: toggle

	.. container:: header

		:arrow:`Rotational Propagator`

	.. include:: RotationalPropagator.rst
