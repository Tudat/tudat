
.. role:: jsontype
.. role:: jsonkey
.. role:: arrow

- :jsontype:`boolean` :jsonkey:`notifyOnPropagationStart` (optional). Whether a message should be printed when the propagation starts. Default value: :literal:`false`.
- :jsontype:`boolean` :jsonkey:`notifyOnPropagationTermination` (optional). Whether a message should be printed when the propagation is terminated with no errors. Default value: :literal:`false`.
- :jsontype:`number` :jsonkey:`printInterval` (optional). Period of time (if any) after which the current state should be outputted to the command-line.
- :jsontype:`boolean` :jsonkey:`defaultValueUsedForMissingKey` (optional). Behaviour of the JSON validator when a default value is used for a missing optional key. Possible values: :literal:`continueSilently`, :literal:`printWarning`, :literal:`throwError`. Default value: :literal:`continueSilently`.
- :jsontype:`boolean` :jsonkey:`unusedKey` (optional). Behaviour of the JSON validator when the propagation starts and there are keys defined in the JSON file that have not been used. Possible values: :literal:`continueSilently`, :literal:`printWarning`, :literal:`throwError`. Default value: :literal:`printWarning`.
- :jsontype:`path` :jsonkey:`fullSettingsFile` (optional). Path where the file (if any) containing the full settings actually used in the propagation, including all the used default values, should be saved.
- :jsontype:`boolean` :jsonkey:`tagOutputFilesIfPropagationFails` (optional). Whether the generated output files should contain the header :literal:`FAILURE` if the propagation terminates before reaching the termination condition. Default value: :literal:`true`.
