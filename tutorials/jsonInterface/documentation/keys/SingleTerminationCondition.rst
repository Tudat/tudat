.. role:: arrow

- :literal:`object` :class:`variable` (mandatory) Variable the termination condition is associated to.

	.. container:: toggle

		.. container:: header

			:arrow:`Variable`

		.. include:: keys/Variable.rst
- :literal:`numeric` :class:`lowerLimit` (mandatory if :literal:`! upperLimit`) The propagation will terminate if the value of :literal:`variable` is smaller than this value (or larger than :literal:`upperLimit` if defined).
- :literal:`numeric` :class:`upperLimit` (mandatory if :literal:`! lowerLimit`) The propagation will terminate if the value of :literal:`variable` is larger than this value (or smaller than :literal:`lowerLimit` if defined).
