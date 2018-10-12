
.. role:: jsontype
.. role:: jsonkey
.. role:: arrow

- :jsontype:`string` :jsonkey:`deformationType` (mandatory). Variable defining the type of gravity field variation. Possible values: :literal:`"basicSolidBody"`, :literal:`"tabulatedVariation"`.
- :jsontype:`string` :jsonkey:`identifier` (optional). Variable denoting the identifier for gravity field variation.
- :jsontype:`number` :jsonkey:`maximumDegree` (mandatory). Value of the maximum degree of the spherical harmonics model, for which an acceleration has to be saved.
- :jsontype:`number` :jsonkey:`maximumOrder` (mandatory). Value of the maximum order of the spherical harmonics model, for which an acceleration has to be saved.
