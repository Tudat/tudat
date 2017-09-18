
.. role:: jsontype
.. role:: jsonkey
.. role:: arrow

- :jsontype:`boolean` :jsonkey:`useAllThrustModels` (optional). Whether all engines of the associated body are to be combined into a single thrust model. Default value: :literal:`true`.
- :jsontype:`string` :jsonkey:`associatedThrustSource` (mandatory if :jsonkey:`useAllThrustModels` set to :literal:`false`). Name of engine model from which thrust is to be derived.
