
.. role:: jsontype
.. role:: jsonkey
.. role:: arrow

- :jsontype:`object` :jsonkey:`map` (mandatory if :jsonkey:`independentVariableValues`, :jsonkey:`dependentVariableValues` or :jsonkey:`file` undefined). Object containing the values of the dependent variable (scalar or vectorial) to be interpolated. The keys are the string representation of the values of the independent variable (scalar).
- :jsontype:`number[ ]` :jsonkey:`independentVariableValues` (optional). List of values of the dependent variable (scalar).
- :jsontype:`number[ ][ ]` :jsonkey:`dependentVariableValues` (optional). List of values of the independent variable (scalar or vectorial).
- :jsontype:`number[ ][ ]` :jsonkey:`dependentVariableFirstDerivativeValues` (mandatory if the interpolator type is :literal:`"hermiteSpline"`). List of values of the first derivatives of the independent variable (scalar or vectorial).
- :jsontype:`path` :jsonkey:`file` (optional). Path to the file containing the data to be interpolated. Each line corresponds to a pair of dependent-independent values. The first column contains the values of the independent variables, while the second column (and subsequent columns for vectorial dependent variables) contains the values of the dependent variable. Empty characters (such as space or tabulator) can be used to separate values.
