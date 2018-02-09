.. _extendingJSON_otherTypes:

.. role:: jsontype
.. role:: jsonkey

Other types
===========

Some types that are not used to store settings require different :literal:`to_json` and :literal:`from_json` implementations. If these types are prone to be used in different parts of the :literal:`json_interface`, their :literal:`to_json` and :literal:`from_json` should be declared in :class:`Tudat/JsonInterface/Support/valueConversions.h`. As discussed previously, there are already some functions in that file, such as custom implementations for :class:`std::map` and :class:`std::unordered_map` that ignore the special keys, for :class:`std::vector` containing non-primitive types, and for :class:`Eigen::Matrix`.

Additionally, there are functions for :class:`std::complex`, so that the following is possible:

.. code-block:: cpp

  nlohmann::json j = "(1,-0.5)";
  std::complex< double > complexNumber = j;         // 1 - 0.5i

For :class:`std::pair`:

.. code-block:: cpp

  nlohmann::json j = R"(
    [
      6,
      "keplerian"
    ]
  )"_json;
  std::pair< int, std::string > pair = j;         // { 6, "keplerian" }

And for :class:`Eigen::Quaterniond`, so that it can be created directly from an :class:`Eigen::Matrix3d` when the provided :class:`nlohmann::json` object is of value type :jsontype:`array`, or using the function :literal:`spice_interface::computeRotationQuaternionBetweenFrames` when it is of value type :jsontype:`object`. Thus, it is possible to create an :class:`Eigen::Quaterniond` from any of these two JSON files:

.. code-block:: json

  [
    [ 1,  0,  0 ],
    [ 0,  1, -1 ],
    [ 0, -1,  1 ]
  ]


.. code-block:: json

  {
    "originalFrame": "ECLIPJ2000",
    "targetFrame": "IAU_Earth",
    "initialTime": 0
  }

