
.. role:: jsontype
.. role:: jsonkey
.. role:: arrow

- :jsontype:`string` :jsonkey:`bodyIdentifier` (mandatory). Parameter identifying for which body an ephemeris is to be created. Possible values: :literal:`"mercury"`, :literal:`"venus"`, :literal:`"earthMoonBarycenter"`, :literal:`"mars"`, :literal:`"jupiter"`, :literal:`"saturn"`, :literal:`"uranus"`, :literal:`"neptune"`, :literal:`"pluto"`.
- :jsontype:`boolean` :jsonkey:`useCircularCoplanarApproximation` (mandatory). Boolean defining whether a circular, coplanar orbit of the body is to be assumed (creating an :class:`ApproximatePlanetPositionsCircularCoplanar` object), or whether a non-zero inclination and long-period changes in the orbit are to be included (creating an :class:`ApproximatePlanetPositions` object).
