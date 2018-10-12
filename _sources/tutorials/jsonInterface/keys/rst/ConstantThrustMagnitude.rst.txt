
.. role:: jsontype
.. role:: jsonkey
.. role:: arrow

- :jsontype:`number` :jsonkey:`constantMagnitude` (mandatory). Constant thrust magnitude that is to be used.
- :jsontype:`number` :jsonkey:`specificImpulse` (mandatory). Constant specific impulse that is to be used.
- :jsontype:`number[ 3 ]` :jsonkey:`bodyFixedDirection` (optional). Direction of thrust force in body-fixed frame. Default value: :literal:`[1, 0, 0]`.
