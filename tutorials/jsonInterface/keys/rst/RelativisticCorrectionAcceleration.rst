
.. role:: jsontype
.. role:: jsonkey
.. role:: arrow

- :jsontype:`boolean` :jsonkey:`calculateSchwarzschildCorrection` (optional). Whether the Schwarzschild term is used. Default value: :literal:`true`.
- :jsontype:`boolean` :jsonkey:`calculateLenseThirringCorrection` (optional). Whether the Lense-Thirring term is used. Default value: :literal:`false`.
- :jsontype:`boolean` :jsonkey:`calculateDeSitterCorrection` (optional). Whether the de Sitter term is used. Default value: :literal:`false`.
- :jsontype:`string` :jsonkey:`primaryBody` (mandatory). Name of primary body (e.g. Sun for acceleration acting on an Earth-orbiting satellite).
- :jsontype:`number[ 3 ]` :jsonkey:`centralBodyAngularMomentum` (optional). Constant angular momentum of central body. Default value: :literal:`[0, 0, 0]`.
