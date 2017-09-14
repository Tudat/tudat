.. role:: arrow

- :literal:`object[]` :class:`anyOf` (mandatory if :literal:`! allOf`) Object containing a list of single-variable termination conditions. The propagation will terminate if any of the conditions is met.

	.. container:: toggle

		.. container:: header

			:arrow:`Single Termination Condition`

		.. include:: keys/SingleTerminationCondition.rst
- :literal:`object[]` :class:`allOf` (mandatory if :literal:`! anyOf`) Object containing a list of single-variable termination conditions. The propagation will terminate only if all the conditions are met.

	.. container:: toggle

		.. container:: header

			:arrow:`Single Termination Condition`

		.. include:: keys/SingleTerminationCondition.rst
