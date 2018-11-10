
.. role:: jsontype
.. role:: jsonkey
.. role:: arrow

- :jsontype:`string` :jsonkey:`bodyExertingAcceleration` (mandatory). Name of body exerting acceleration.
- :jsontype:`number[ ][ ]` :jsonkey:`componentIndices` (mandatory). List of indices for degree/order that is to be saved. For instance, when saving the degree/order terms 0/0, 2/0 and 2/2, use ``"componentIndices": [ [ 0, 0 ], [ 2, 0 ], [ 2, 2 ] ]``

