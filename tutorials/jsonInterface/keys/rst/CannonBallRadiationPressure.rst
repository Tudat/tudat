
.. role:: jsontype
.. role:: jsonkey
.. role:: arrow

- :jsontype:`number` :jsonkey:`radiationPressureCoefficient` (mandatory). Radiation pressure coefficient.
- :jsontype:`number` :jsonkey:`referenceArea` (optional). Surface area that undergoes radiation pressure. If not defined, it will be retrieved from the body's :jsonkey:`referenceArea` key to which the map of radiation pressure objects is assigned. Note that for the cannon ball radiation pressure to be calculated, the reference area has to be defined at least once, either at the body level or in the cannon ball radiation pressure settings. The definition of the reference area is optional only if it has already been defined in the body's keys.
