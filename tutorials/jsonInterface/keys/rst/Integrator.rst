
.. role:: jsontype
.. role:: jsonkey
.. role:: arrow

- :jsontype:`string` :jsonkey:`type` (mandatory). Method to use for the integration of the equations of motion. Possible values: :literal:`"euler"`, :literal:`"rungeKutta4"`, :literal:`"rungeKuttaVariableStepSize"`. Default value: :literal:`"rungeKutta4"`.
- :jsontype:`number` :jsonkey:`initialTime` (mandatory). Initial value of the independent variable (typically, seconds since J2000).
- :jsontype:`number` :jsonkey:`saveFrequency` (optional). Frequency with which the state is saved to the results map. For instance, a value of :literal:`5` means that the results are only saved for 1 in every 5 integration steps. Default value: :literal:`1`.
- :jsontype:`boolean` :jsonkey:`assessPropagationTerminationConditionDuringIntegrationSubsteps` (optional). Whether the termination conditions defined in :jsonkey:`.termination` should be assessed after each integrator substep (e.g. for an RK4 integrator, after computing :literal:`k1`, :literal:`k2`, :literal:`k3`...). Not applicable to the Euler integrator. Default value: :literal:`false`.

.. container:: toggle

	.. container:: header

		:arrow:`Euler Integrator`

	.. include:: rst/EulerIntegrator.rst

.. container:: toggle

	.. container:: header

		:arrow:`Runge-Kutta 4 Integrator`

	.. include:: rst/RungeKutta4Integrator.rst

.. container:: toggle

	.. container:: header

		:arrow:`Runge-Kutta Variable Step-Size Integrator`

	.. include:: rst/RungeKuttaVariableStepSizeIntegrator.rst
