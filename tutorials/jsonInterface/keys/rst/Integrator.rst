
.. role:: jsontype
.. role:: jsonkey
.. role:: arrow

- :jsontype:`string` :jsonkey:`type` (optional). Method to use for the integration of the equations of motion. Possible values: :literal:`"euler"`, :literal:`"rungeKutta4"`, :literal:`"rungeKuttaVariableStepSize"`. Default value: :literal:`"rungeKutta4"`.
- :jsontype:`number` :jsonkey:`initialTime` (optional). Initial value of the independent variable (typically, seconds since J2000). If not provided, will be obtained from :jsonkey:`.initialEpoch`.
- :jsontype:`number` :jsonkey:`saveFrequency` (optional). Frequency with which the state is saved to the results map. For instance, a value of :literal:`5` means that the results are only saved for 1 in every 5 integration steps. Default value: :literal:`1`.
- :jsontype:`boolean` :jsonkey:`assessPropagationTerminationConditionDuringIntegrationSubsteps` (optional). Whether the termination conditions defined in :jsonkey:`.termination` should be assessed after each integrator substep (e.g. for an RK4 integrator, after computing :literal:`k1`, :literal:`k2`, :literal:`k3`...). Not applicable to the Euler integrator. Default value: :literal:`false`.

.. container:: toggle

	.. container:: header

		:arrow:`Euler Integrator`

	.. include:: EulerIntegrator.rst

.. container:: toggle

	.. container:: header

		:arrow:`Runge-Kutta 4 Integrator`

	.. include:: RungeKutta4Integrator.rst

.. container:: toggle

	.. container:: header

		:arrow:`Runge-Kutta Variable Step-Size Integrator`

	.. include:: RungeKuttaVariableStepSizeIntegrator.rst
