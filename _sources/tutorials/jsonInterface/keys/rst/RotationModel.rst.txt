
.. role:: jsontype
.. role:: jsonkey
.. role:: arrow

- :jsontype:`string` :jsonkey:`type` (mandatory). Possible values: :literal:`"simple"`, :literal:`"spice"`.
- :jsontype:`string` :jsonkey:`originalFrame` (mandatory). Identifier of the base frame of rotation model.
- :jsontype:`string` :jsonkey:`targetFrame` (mandatory). Identifier of the target frame of rotation model.

.. container:: toggle

	.. container:: header

		:arrow:`Simple Rotation Model`

	.. include:: SimpleRotationModel.rst
