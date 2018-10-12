
.. role:: jsontype
.. role:: jsonkey
.. role:: arrow

- :jsontype:`object[ ]` :jsonkey:`anyOf` (mandatory if :jsonkey:`allOf` undefined). Object containing a list of single and/or multiple termination conditions. The propagation will terminate if any of the conditions is met.
- :jsontype:`object[ ]` :jsonkey:`allOf` (mandatory if :jsonkey:`anyOf` undefined). Object containing a list of single and/or multiple termination conditions. The propagation will terminate only if all the conditions are met.
