
.. role:: jsontype
.. role:: jsonkey
.. role:: arrow

- :jsontype:`path` :jsonkey:`file` (mandatory). Path to a PDS file containing the cosine and sine coefficients. Each row contains a degree, order, cosine coefficient and sine coefficient.
- :jsontype:`string` :jsonkey:`associatedReferenceFrame` (mandatory). Identifier of the associated reference frame, e.g. :literal:`IAU_Earth`.
- :jsontype:`number` :jsonkey:`maximumDegree` (mandatory). Maximum degree of the coefficients to be loaded from the file.
- :jsontype:`number` :jsonkey:`maximumOrder` (mandatory). Maximum order of the coefficients to be loaded from the file.
- :jsontype:`number` :jsonkey:`gravitationalParameterIndex` (mandatory). Index (starting from zero) at which the associated gravitational parameter of the model can be found in a vector resulting from splitting the first line of the input file (aka header). Must be set to :literal:`-1` if the file has no header. Default value: :literal:`0`.
- :jsontype:`number` :jsonkey:`referenceRadiusIndex` (mandatory). Index (starting from zero) at which the associated reference radius of the model can be found in a vector resulting from splitting the header. Must be set to :literal:`-1` if the file has no header. Default value: :literal:`1`.
- :jsontype:`number` :jsonkey:`gravitationalParameter` (mandatory if :jsonkey:`gravitationalParameterIndex` set to :literal:`-1`). Value of the gravitational parameter, to be provided if the file has no header.
- :jsontype:`number` :jsonkey:`referenceRadius` (mandatory if :jsonkey:`referenceRadiusIndex` set to :literal:`-1`). Value of the reference radius, to be provided if the file has no header.
