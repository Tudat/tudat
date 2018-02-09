
.. role:: jsontype
.. role:: jsonkey
.. role:: arrow

- :jsontype:`object` :jsonkey:`variable` (mandatory). Variable the termination condition is associated to. Currently, a termination condition can only be defined as a function of the independent variable, the CPU-time variable or a dependent variable. For the independent and CPU-time variables, the key :jsonkey:`lowerLimit` is ignored.

	.. container:: toggle

		.. container:: header

			:arrow:`Variable`

		.. include:: Variable.rst
- :jsontype:`number` :jsonkey:`lowerLimit` (mandatory if :jsonkey:`upperLimit` undefined). The propagation will terminate if the value of :jsonkey:`variable` is smaller than this value (or larger than :jsonkey:`upperLimit` if defined).
- :jsontype:`number` :jsonkey:`upperLimit` (mandatory if :jsonkey:`lowerLimit` undefined). The propagation will terminate if the value of :jsonkey:`variable` is larger than this value (or smaller than :jsonkey:`lowerLimit` if defined).
