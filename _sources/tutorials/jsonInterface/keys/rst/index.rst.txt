.. _jsonInterface_keys:

.. role:: jsontype
.. role:: jsonkey

.. role:: showall
.. role:: hideall
.. role:: arrow

List of Keys
============

In this page, an exhaustive list with the possible keys that can be defined in input JSON files for the :literal:`json_interface` application is provided. These are used to create the settings as discussed in :ref:`tudatFeaturesIndex`. 

The type :jsontype:`path` before a key name denotes that a :jsontype:`string` such as :literal:`"@path(relPath)"` has to be provided when using relative paths inside :ref:`jsonInterface_modularFiles`. When providing absolute paths or not using modular files, both :literal:`"@path(file)"` and :literal:`"file"` will work.

:jsontype:`T[ ]` denotes that an array of :jsontype:`T` is expected. If no number is indicated, there is no constraint on the size of the array. For instance, :jsontype:`number[ 6 ]` is used for state vectors. A second :jsontype:`[ ]` is used to represent matrices, e.g. :jsontype:`number[ ][ ]` corresponds to :class:`Eigen::MatrixXd`, and :jsontype:`number[ 3 ][ 3 ]` corresponds to :class:`Eigen::Matrix3d`. When a vector is expected (i.e. an :class:`Eigen::Matrix` with only one row and/or one column), it is possible to provide it as a an array of numbers (e.g. :literal:`[ 1, 2, 3 ]`) only if at least one of the dimensions of the expected :class:`Eigen::Matrix` is of fixed size, and this size matches the length of the provided array. If the expected matrix has dynamic size (both for the number of rows and columns), it is mandatory to provide an array of arrays (e.g. :literal:`[ [ 1, 2, 3 ] ]` will be converted to a matrix of dynamic size 1x3, while :literal:`[ [ 1 ], [ 2 ], [ 3 ] ]` will be converted to a matrix of dynamic size 3x1).

Some keys appear twice because they accept different value types (and a different description for each case is provided). For instance, the body's initial state can be of value type :jsontype:`number[ 6 ]` (where a Cartesian state is assumed), or an :jsontype:`object` containing keys such as :jsonkey:`x`, :jsonkey:`vy`, :jsonkey:`semiMajorAxis`, :jsonkey:`inclination`, etc.

In some parts, key paths are used to refer to the key of a given object in the JSON tree. The different keys are separated by dots. For instance :jsonkey:`atmosphere.type` refers to the key :jsonkey:`type` of an :jsontype:`object` assigned to the key :jsonkey:`atmosphere`. When the key path is absolute (or in other words relative to the JSON defined in the main input file), the key path starts with a dot (for instance, :jsonkey:`.initialEpoch` refers to the key :jsonkey:`initialEpoch` defined directly in the root object of the main input file). When the key path does not start with a dot, it is relative to the context in which it is mentioned. For arrays, the index (starting from zero) is used, so :jsonkey:`centralBodies[ 1 ]` refers to the second body in the :jsontype:`array` assigned to the key :jsonkey:`centralBodies`.

:showall:`Unfold All` :hideall:`Fold All`

.. include:: root.rst
