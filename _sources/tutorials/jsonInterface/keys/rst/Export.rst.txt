
.. role:: jsontype
.. role:: jsonkey
.. role:: arrow

- :jsontype:`path` :jsonkey:`file` (mandatory). Path to the file in which to save the results. The file will contain a matrix in which each row corresponds to an integration step and each column to the values of the variables in the same order as provided in :jsonkey:`variables`.
- :jsontype:`object[ ]` :jsonkey:`variables` (mandatory). List of variables that are going to be saved.

	.. container:: toggle

		.. container:: header

			:arrow:`Variable`

		.. include:: Variable.rst
- :jsontype:`string` :jsonkey:`header` (optional). Header to be included in the first line of the output file.
- :jsontype:`boolean` :jsonkey:`epochsInFirstColumn` (optional). Whether to include the epochs in the first column of the results matrix. Default value: :literal:`true`.
- :jsontype:`boolean` :jsonkey:`onlyInitialStep` (optional). Whether to save only the values corresponding to the initial integration step. Does not override :jsonkey:`onlyFinalStep`. Default value: :literal:`false`.
- :jsontype:`boolean` :jsonkey:`onlyFinalStep` (optional). Whether to save only the values corresponding to the initial integration step. Does not override :jsonkey:`onlyInitialStep`. Default value: :literal:`false`.
- :jsontype:`number` :jsonkey:`numericalPrecision` (optional). Maximum number of significant digits to be used for decimal numbers. Default value: :literal:`15`.
- :jsontype:`number` :jsonkey:`printVariableIndicesToTerminal` (optional). Boolean that defines whether indices at which in output vector, at which variables are saved, will be printed to the terminal Default value: :literal:`false`.
