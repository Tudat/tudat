
.. role:: jsontype
.. role:: jsonkey
.. role:: arrow

- :jsontype:`number` :jsonkey:`initialTime` (mandatory). Time at which :jsonkey:`initialOrientation` represents the instantaneous rotation.
- :jsontype:`number` :jsonkey:`rotationRate` (mandatory). Rotation rate of body about its local z-axis.
- :jsontype:`number[ 3 ][ 3 ]` :jsonkey:`initialOrientation` (optional). Rotation from base to target frame at :jsonkey:`initialTime`. If not specified, will be computed using Spice.
