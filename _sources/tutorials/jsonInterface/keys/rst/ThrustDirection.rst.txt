
.. role:: jsontype
.. role:: jsonkey
.. role:: arrow

- :jsontype:`string` :jsonkey:`type` (mandatory). Possible values: :literal:`"colinearWithStateSegment"`, :literal:`"fromExistingBodyOrientation"`.
- :jsontype:`string` :jsonkey:`relativeBody` (mandatory). Name of the body relative to which thrust guidance algorithm is defined.

.. container:: toggle

	.. container:: header

		:arrow:`Colinear With State Segment Thrust Direction`

	.. include:: ColinearWithStateSegmentThrustDirection.rst
