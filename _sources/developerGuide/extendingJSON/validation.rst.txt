.. _extendingJSON_validation:

.. role:: jsontype
.. role:: jsonkey

Validation
==========

The validation of the JSON input files provided by the user is made possible by the :literal:`getValue< ExpectedType >` function and special keys and :class:`KeyPath` s, as explained in :ref:`extendingJSON_enhancedValueAccess`. This validator is customisable, i.e. the user can define a set of keys in their input file that will determine the behaviour of the validator under certain circumstances.


Validator options
~~~~~~~~~~~~~~~~~

The following keys of the :literal:`mainJson` object can be defined to customise the behaviour of the validator:

  - :jsonkey:`options.defaultValueUsedForMissingKey`: determines the behaviour of the validator when the :literal:`getValue( const nlohmann::json& j, const KeyPath& keyPath, const ValueType& defaultValue )` function is called and :literal:`keyPath` is undefined. By default, :literal:`defaultValue` is used without informing the user.

  - :jsonkey:`options.unusedKey`: determines the behaviour of the validator when the :literal:`void checkUnusedKeys( const nlohmann::json& mainJson, ... )` function is called and some of the keys defined in :literal:`mainJson` have not been accessed. This function is called by the :literal:`json_interface` application right before integrating the equations of motion, when all the settings objects have been created. By default, a warning is printed for each unused key.

These three options are converted to a value of the enum :class:`ExceptionResponseType` defined in :class:`Tudat/JsonInterface/Support/errorHandling.h`. The three possible values are:

  - :literal:`continueSilently`: allow this validator feature and do not inform the user.
  - :literal:`printWarning`: allow this validator feature but print a warning when using it.
  - :literal:`throwError`: do not allow this validator feature and terminate throwing an error when trying to use this feature. For :jsonkey:`defaultValueUsedForMissingKey`, an :class:`UndefinedKeyError` will be thrown (i.e. the behaviour of the :literal:`getValue` function with and without third default value argument will be equivalent). For :jsonkey:`unusedKey`, the error message :literal:`Validation failed because there are unused keys.` will be printed after providing the list of unused keys.


Access history
~~~~~~~~~~~~~~

The validator feature to check the keys that have not been used is useful to inform the user about the existence of redundant information in their input file, the use of old key identifiers that are not used anymore in the current version of Tudat, or the detection of typos in optional keys when manually writing the input file. Consider the following:

.. code-block:: cpp

    nlohmann::json j = { { "tyep", "euler" }, { "stepSize": 20 } };  // not typo in key "tyep"
    std::string type = getValue( j, "type", "rungeKutta4" );

This would result in the integrator being used to be of the default type :literal:`rungeKutta4` even though the intention of the user was to use of type :literal:`euler`. However, the user will not know about this unless they have configured the validator to warn about the use of default values (which is off by default). Thus, before integrating the equations of motion, when all the required information has been retrieved from the :literal:`mainJson` object, the call to :literal:`checkUnusedKeys` will print the following warning (by default):

.. code-block:: txt

    Unused key: integrator.tyep

allowing the user to quickly identify and fix the problem. The validator can also be configured not to print a warning (this is highly discouraged) or to terminate in case that there are unused keys.

In order to make this feature possible, the :literal:`json_interface` has to keep track of all the key paths of :literal:`mainJson` that have been accessed, and then iterate over all its keys to check whether any of them has not been used. Keeping track of the accessed key paths is done by calling the :literal:`getValue` function (i.e. using the default value access methods will not keep track of the used keys). The original implementation proposal consisted in defining a special key, e.g. :literal:`#accessedKeyPaths`, for the :literal:`mainJson` object, in which an array of :class:`KeyPath` s would be stored. However, the :literal:`mainJson` object is not always accessible in the :literal:`from_json` functions. For instance, in the following code it *is* accessible:

.. code-block:: cpp
  :caption: :class:`simulation.h`
  :name: simulation-h
  
  #include "integrator.h"
  
  void from_json( const nlohmann::json& mainJson, Simulation& simulation )
  {
      simulation.integrator =
        getValue< Integrator >( mainJson, "integrator" );  // integrators::from_json called
  }

so the :literal:`getValue` function could potentially modify :literal:`mainJson` to add :literal:`integrator` as an accessed key, but here the :literal:`mainJson` object is not accessible:

.. code-block:: cpp
  :caption: :class:`integrator.h`
  :name: integrator-h
  
  namespace integrators
  {
      void from_json( const nlohmann::json& jsonIntegrator, Integrator& integrator )
      {
          integrator.type = getValue< std::string >( jsonIntegrator, "type" );
          integrator.stepSize = getValue< double >( jsonIntegrator, "stepSize" );
      }
  }

and thus is is not possible for the :literal:`getValue` to add the :literal:`integrator.type` and :literal:`integrator.stepSize` key paths to the :literal:`mainObject`'s :literal:`#accessedKeyPaths` special key.

A possible walk-around this issue could consist in defining the :literal:`#accessedKeyPaths` for each :class:`nlohmann::json` object, and not only for the :literal:`mainJson`. However, the updated sub-:class:`nlohmann::json` objects would still have to be passed back to the :literal:`mainJson` somehow. Even if this was possible, the definition of the :literal:`from_json` function, with the first argument as a :literal:`const nlohmann::json&` makes this approach unfeasible, as the :literal:`getValue` function would not be capable of updating it (it will not be possible to update its :literal:`#accessedKeyPaths` because it is a constant object). Thus, a completely different approach had to be followed, making use of a global variable. This global variable is declared in :class:`Tudat/JsonInterface/Support/valueAccess.h`:

.. code-block:: cpp

  extern std::set< KeyPath > accessedKeyPaths;

This variable is automatically updated when calling the :literal:`getValue` function. In order to clear the contents of this variable, the function :literal:`clearAccessHistory` must be called. In the :literal:`json_interface` application, this is done after reading the input JSON file and just before starting to use the read :class:`nlohmann::json` object to update the settings objects.

.. warning:: The current implementation has one limitation: it is not possible to keep track of the accessed keys of multiple :class:`nlohmann::json` objects simultaneously. Currently, this is not done anywhere in the JSON Interface library, but if in the feature this is required, it will be necessary to create a derived class of :class:`nlohmann::json` (e.g. :class:`EnhancedJSON`) with the :literal:`accessedKeyPaths` as property and the functions :literal:`getValue` and :literal:`checkUnusedKeys` as methods (probably other global functions declared in :class:`Tudat/JsonInterface/Support/valueAccess.h` would also have to be moved to this class). The :literal:`to_json` and :literal:`from_json` methods would have to be updated to take objects of this class as first argument instead of the basic :class:`nlohmann::json` objects. An attempt to implement this was done at one point during the development of the :literal:`json_interface`, but it was unsuccessful due to the existence of an `inheritance bug <https://github.com/nlohmann/json/issues/608>`_ in the JSON library.
