.. role:: arrow

- :literal:`path` :class:`file` (mandatory) Path to the file in which to save the results. The file will contain a matrix in which each row corresponds to an integration step and each column to the values of the variables in the same order as provided in :literal:`variables`.
- :literal:`object[]` :class:`variables` (mandatory) List of variables that are going to be saved.

	.. container:: toggle

		.. container:: header

			:arrow:`Variable`

		.. include:: keys/Variable.rst
- :literal:`string` :class:`header` (optional) Header to be included in the first line of the output file.
- :literal:`bool` :class:`epochsInFirstColumn` (optional) Whether to include the epochs in the first column of the results matrix. Default value: :literal:`true`.
- :literal:`bool` :class:`onlyInitialStep` (optional) Whether to save only the values corresponding to the initial integration step. Does not override :literal:`onlyFinalStep`. Default value: :literal:`false`.
- :literal:`bool` :class:`onlyFinalStep` (optional) Whether to save only the values corresponding to the initial integration step. Does not override :literal:`onlyInitialStep`. Default value: :literal:`false`.
- :literal:`int` :class:`numericalPrecision` (optional) Maximum number of significant digits to be used for decimal numbers. Default value: :literal:`16`.
