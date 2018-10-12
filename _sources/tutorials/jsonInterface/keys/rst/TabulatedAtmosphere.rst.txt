
.. role:: jsontype
.. role:: jsonkey
.. role:: arrow

- :jsontype:`map` :jsonkey:`file` (mandatory). Map of files containing atmospheric properties, the map should contain as many files as independent variables and each file should contain as many dimensions as dependent variables.
- :jsontype:`string[ ]` :jsonkey:`independentVariablesNames` (mandatory). Names of independent variables included in the tabulated atmosphere. Possible values: :literal:`"altitude"`, :literal:`"longitude"`, :literal:`"latitude"`, :literal:`"time"`. Default value: :literal:`"altitude"`.
- :jsontype:`string[ ]` :jsonkey:`dependentVariablesNames` (mandatory). Names of dependent variables included in the tabulated atmosphere. Possible values: :literal:`"density"`, :literal:`"pressure"`, :literal:`"temperature"`, :literal:`"gasConstant"`, :literal:`"specificHeatRatio"`, :literal:`"molarMass"`. Default values: :literal:`"density"`, :literal:`"pressure"`, :literal:`"temperature"`.
- :jsontype:`number` :jsonkey:`specificGasConstant` (optional). Specific gas constant for (constant) atmospheric chemical composition. Default value: 287.
- :jsontype:`number` :jsonkey:`ratioOfSpecificHeats` (optional). Ratio of specific heats for (constant) atmospheric chemical composition. Default value: 1.4.
- :jsontype:`string[ ]` :jsonkey:`boundaryHandling` (optional). Boundary handling methods, in case an independent variable is outside the specified range. Possible values: :literal:`"throwException"`, :literal:`"boundaryValue"`, :literal:`"boundaryValueWithWarning"`, :literal:`"extrapolate"`, :literal:`"extrapolateWithWarning"`. Default value: :literal:`"boundaryValue"`.