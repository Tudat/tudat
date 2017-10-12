
.. role:: jsontype
.. role:: jsonkey
.. role:: arrow

- :jsontype:`number` :jsonkey:`centralBodyAverageRadius` (optional). If not specified, will be retrieved from the corresponding central body specified in the propagator.
- :jsontype:`number` :jsonkey:`epoch` (optional). If not specified, will be retrieved from the initial epoch of the integrator.
- :jsontype:`number` :jsonkey:`radius` (mandatory if :jsonkey:`altitude` undefined).
- :jsontype:`number` :jsonkey:`altitude` (mandatory if :jsonkey:`radius` undefined).
- :jsontype:`number` :jsonkey:`latitude` (optional). Default value: :literal:`0`.
- :jsontype:`number` :jsonkey:`longitude` (optional). Default value: :literal:`0`.
- :jsontype:`number` :jsonkey:`speed` (mandatory).
- :jsontype:`number` :jsonkey:`flightPathAngle` (optional). Default value: :literal:`0`.
- :jsontype:`number` :jsonkey:`headingAngle` (optional). Default value: :literal:`0`.
