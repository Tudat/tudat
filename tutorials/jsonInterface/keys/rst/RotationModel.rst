
.. role:: jsontype
.. role:: jsonkey
.. role:: arrow

- :jsontype:`string` :jsonkey:`type` (mandatory). Possible values: :literal:`"simple"`, :literal:`"spice"`.
- :jsontype:`string` :jsonkey:`originalFrame` (mandatory). Identifier of the base frame of rotation model.
- :jsontype:`string` :jsonkey:`targetFrame` (mandatory; except for GCRS to ITRS, where "GCRS" is hard-coded). Identifier of the target frame of rotation model.

.. container:: toggle

	.. container:: header

		:arrow:`Simple Rotation Model`

	.. include:: SimpleRotationModel.rst

.. container:: toggle

	.. container:: header

		:arrow:`GCRS to ITRS (high-accuracy Earth model)`

	.. include:: ITRS.rst
