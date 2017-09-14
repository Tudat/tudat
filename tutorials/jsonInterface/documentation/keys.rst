.. _jsonInterface_documentation_keys:

List of keys
============

.. role:: showall

.. role:: hideall

In this page, ...

The type :literal:`path` denotes that a :literal:`string` such as :literal:`"@path(relPath)"` has to be provided when using relative paths inside modular files. When providing absolute paths or not using modular files, both :literal:`"@path(file)"` and :literal:`"file"` will work.

:literal:`T[]` denotes that an array of :literal:`T` is expected. If unidimensional array inference is supported, a single :literal:`T` is also valid. If no number is indicated, there is no constraint on the size of the array. For instance, :literal:`numeric[6]` is used for state vectors.

Some keys appear twice because they accept different value types (and a different description for each case is provided). For instance, the body's initial state can be of value type :literal:`numeric[6]` (where a Cartesian state is assumed), or an :literal:`object` containing keys such as :class:`x`, :class:`vy`, :class:`semiMajorAxis`, :class:`inclination`, etc.

:showall:`Expand All` :hideall:`Collapse All`

.. include:: keys/Simulation.rst

