.. _jsonInterface_keys:

.. role:: jsontype
.. role:: jsonkey

.. role:: showall
.. role:: hideall
.. role:: arrow

List of keys
============

In this page, an exhaustive list with the possible keys that can be defined in input JSON files for the :literal:`json_interface` application is provided.

The type :jsontype:`path` before a key name denotes that a :jsontype:`string` such as :literal:`"@path(relPath)"` has to be provided when using relative paths inside modular files. When providing absolute paths or not using modular files, both :literal:`"@path(file)"` and :literal:`"file"` will work.

:jsontype:`T[]` denotes that an array of :jsontype:`T` is expected. If unidimensional array inference is supported, a single :jsontype:`T` is also valid. If no number is indicated, there is no constraint on the size of the array. For instance, :jsontype:`number[6]` is used for state vectors. A second :jsontype:`[]` is used to represent matrices, e.g. :jsontype:`number[][]` will be converted to :class:`Eigen::MatrixXd`, and :jsontype:`number[3][3]` will be converted to :class:`Eigen::Matrix3d`.

Some keys appear twice because they accept different value types (and a different description for each case is provided). For instance, the body's initial state can be of value type :jsontype:`number[6]` (where a Cartesian state is assumed), or an :jsontype:`object` containing keys such as :jsonkey:`x`, :jsonkey:`vy`, :jsonkey:`semiMajorAxis`, :jsonkey:`inclination`, etc.

In some parts, key paths are used to refer to the key of a given object in the JSON tree. The different keys are separated by dots. For instance :jsonkey:`atmosphere.type` refers to the key :jsonkey:`type` of an :jsontype:`object` assigned to the key :jsonkey:`atmosphere`. When the key path is absolute (or in other words relative to the JSON defined in the main input file), the key path starts with a dot (for instance, :jsonkey:`.initialEpoch` refers to the key :jsonkey:`initialEpoch` defined directly in the root object of the main input file). When the key path does not start with a dot, it is relative to the context in which it is mentioned. For arrays, the index (starting from zero) is used, so :jsonkey:`centralBodies[1]` refers to the second body in the :jsontype:`array` assigned to the key :jsonkey:`centralBodies`.

:showall:`Unfold All` :hideall:`Fold All`

.. include:: rst/root.rst

