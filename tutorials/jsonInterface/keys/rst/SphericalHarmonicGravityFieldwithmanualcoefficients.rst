
.. role:: jsontype
.. role:: jsonkey
.. role:: arrow

- :jsontype:`number` :jsonkey:`gravitationalParameter` (mandatory). Value of the associated gravitational parameter.
- :jsontype:`number` :jsonkey:`referenceRadius` (mandatory). Value of the associated reference radius.
- :jsontype:`number[ ][ ]` :jsonkey:`cosineCoefficients` (mandatory). Matrix containing the cosine spherical harmonic coefficients (geodesy-normalised). Rows correspond to degrees and columns to orders.
- :jsontype:`number[ ][ ]` :jsonkey:`sineCoefficients` (mandatory). Matrix containing the sine spherical harmonic coefficients (geodesy-normalised). Rows correspond to degrees and columns to orders.
- :jsontype:`string` :jsonkey:`associatedReferenceFrame` (mandatory). Identifier of the associated reference frame, e.g. :literal:`IAU_Earth`.
