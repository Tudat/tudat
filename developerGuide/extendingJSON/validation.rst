.. _extendingJSON_validation:

Validation
==========

The validation of the JSON input files provided by the user is made possible by using the :literal:`getValue< ExpectedType >` function and special keys and :class:`KeyPath` s, as explained in [REF]. This validator is customizable, i.e. the user can define a set of keys in theit input file that will determine the behaviour of the validator under certain events.


Validator options
~~~~~~~~~~~~~~~~~

The following keys of the :literal:`mainJson` object can be defined to customize the behaviour of the validator:

  - :literal:`options.defaultValueUsedForMissingKey`: determines the behaviour of the validator when the :literal:`getValue( const json& j, const KeyPath& keyPath, const ValueType& defaultValue )` function is called and :literal:`keyPath` is undefined. By default, :literal:`defaultValue` is used without informing the user.

  - :literal:`options.unusedKey`: determines the behaviour of the validator when the :literal:`void heckUnusedKeys( const json& mainJson, ... )` function is called and some of the keys defined in :literal:`mainJson` have not been accessed. This function is called right before integrating the equations of motion, when all the settings objects have been created. By default, a warning is printed for each unused key.

  - :literal:`options.unidimensionalArrayInference`: determines the behaviour of the validator when the expected type is :literal:`std::vector< T >` and the object provided by the user is not convertible to :literal:`std::vector< T >` but *is* convertible to :literal:`T`. By default, an :literal:`std::vector< T >` containing only one element (the provided element of type :literal:`T`) is returned without informing the user.

These three options are converted to a value of the enum :class:`ExceptionResponseType` defined in :class:`Tudat/External/JsonInterface/Support/errorHandling.h`. The three possible values are:

  - :literal:`continueSilently`: allow this validator feature and do not inform the user.
  - :literal:`printWarning`: allow this validator feature but print a warning when using it.
  - :literal:`throwError`: do not allow this validator feature and terminate throwing an error when trying to use this feature. For :literal:`defaultValueUsedForMissingKey`, an :class:`UndefinedKeyError` will be thrown (i.e. the behaviour of the :literal:`getValue` function with and without third default value argument will be equivalent). For :literal:`unidimensionalArrayInference`, an :class:`IllegalValueError` will be thrown. For :literal:`unusedKey`, the error :literal:`std::runtime_error( "Validation failed because there are unsued keys." )` will be thrown.


Access history
~~~~~~~~~~~~~~

The validator feature to check the keys that have not been used is usefult to inform the user about the existence of redundant information in their input file, the use of old key identifeirs that are not used anymore in the current version of Tudat, or the detection of typos when manually writing the input file. This feature would not be necessary if defaultable properties did not exist. However, consider the following:

.. code-block:: cpp

    json j = { { "tyep", "euler" }, { "stepSize": 20 } };    // not the typo in the key "tyep"
    std::string type = getValue( j, "type", "rungeKutta4" );

This would result in the integrator being used to by of the default type :literal:`rungeKutta4` even though the intention of the user was to use of type :literal:`euler`. However, the user will not know about this unless they have configured the validator to warn about the use of default values (wich is off by default). Thus, before integratin the equations of motion, when all the required information has been obtained from the :literal:`mainJson` object, the call to :literal:`checkUnusedKeys` will print the following warning (by default):

.. code-block:: txt

    Unused key: integrator.tyep

allowing the user to quickly identify and fix the problem. The validator can also be configured not to print a warning (this is highly discouraged) or to terminate in case that there are unused keys.

In order to make this feature possible, the :literal:`json_interface` has to keep track of all the key paths of :class:`mainJson` that have been accessed, and then iterate over all its keys to check whether any of them has not been used. Keeping track of the accessed key paths is done by calling the :literal:`getValue` function (i.e. using the default value access methods will not keep track of the used keys). The original implementation proposal consisted in defining a special key, e.g. :literal:`#accessedKeyPaths`, for the :literal:`mainJson` object, in which an array of :class:`KeyPath` s would be stored. However, the :literal:`mainJson` object is not always accessible in the :literal:`from_json` functions. For instance, in the following code it *is* accessible:

.. code-block:: cpp
  :caption: :class:`simulation.h`
  :name: simulation-h
  
  #include "integrator.h"
  
  void from_json( const json& mainJson, Simulation& simulation )
  {
      simulation.integrator =
        getValue< Integrator >( mainJson, "integrator" );  // integrators::from_json called
  }

so the getValue function could potentially modify mainJson to add :literal:`integrator` as an accessed key, but here the :literal:`mainJson` object is not accessible:

.. code-block:: cpp
  :caption: :class:`integrator.h`
  :name: integrator-h
  
  namespace integrators
  {
      void from_json( const json& jsonIntegrator, Integrator& integrator )
      {
          integrator.type = getValue< std::string >( jsonIntegrator, "type" );
          integrator.stepSize = getValue< double >( jsonIntegrator, "stepSize" );
      }
  }

and thus is is not possible for the :literal:`getValue` to add the :literal:`integrator.type` and :literal:`integrator.stepSize` key paths to the :literal:`mainObject`'s :literal:`#accessedKeyPaths` special key.

A possible walk-around this issue could consist in defining the :literal:`#accessedKeyPaths` for each :class:`json` object, and not only for the :literal:`mainJson`. However, the updated sub-:class:`json` objects would still have to be passed back to the :literal:`mainJson` somehow. Even if this was possible, the definition of the :literal:`from_json` function, with the first argument as a :literal:`const json&` renders this approach unfeasible, as the :literal:`getValue` function would not be capable of updating it (it will not be possible to update its :literal:`#accessedKeyPaths` because it is a constant object). Thus, a completely different approach had to be followed, making use of a global variable. This global variable is declared in :class:`Tudat/External/JsonInterface/Support/valueAccess.h`:

.. code-block:: cpp

  extern std::set< KeyPath > accessedKeyPaths;

This variable is automatically updated when calling the :literal:`getValue` function. In order to clear the contents of this variable, the function :literal:`clearAccessHistory` must be called. This is done right after reading a JSON file.

.. warning:: The current implementation has one limitation: it is not possible to keep track of the accessed keys of multiple :literal:`mainJson` objects simultaneously. Currently, this is not done anywhere in the :literal:`json_interface`, but if in the feature this is required, it will be neccessary to create a derived class of :class:`json` (e.g. :class:`EnhancedJSON`) with the :literal:`accessedKeyPaths` as property and the functions :literal:`getValue` and :literal:`checkUnusedKeys` as methods (probably other global functions declared in :class:`Tudat/External/JsonInterface/Support/valueAccess.h` would also have to be moved to this class). The :literal:`to_json` and :literal:`from_json` methods would have to be updated to take objects of this class as first argument instead of the basic :class:`json` objects. An attempt to implement this was done at one point during the development of the :literal:`json_interface`, but it was unsuccessful due to the existance of an `inheritance bug <https://github.com/nlohmann/json/issues/608>`_ in the JSON library.


Unidimensional array inference
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

As explained above, unidimensional array inference is a capability of the validator to generate a :literal:`std::vector< T >` with one single element when the user provides an object convertible to type :literal:`T` but an object convertible to type :literal:`std::vector< T >` is expected. This feature is not implemented in the :literal:`getValue` function, but in the custom implementation of the :class:`std::vector`'s :literal:`from_json` function. This means that the following will work:

.. code-block:: cpp

  json person = { { "name", "Peter" }, { "children", "Marc" } };
  std::vector< std::string > childrenEnhancedAccess =
    getValue< std::vector< std::string > >( person, "children" );            // { "Marc" }
  std::vector< std::string > childrenBaiscAccess = person[ "children" ];     // { "Marc" }

However, the :literal:`getValue` function *is* responsible for checking whether unidimensional array inference has been applied, and for printing a warning or throwing an error depending on the user setting for the key :literal:`options.unidimensionalArrayInference` of the main object. This is done by performing the following steps:

  - Get the :literal:`json` object at the requested key path -> :literal:`originalSubjson`.
  - Convert :literal:`originalSubjson` to the requested :literal:`ValueType` -> :literal:`returnObject`.
  - Convert :literal:`returnObject` back to :literal:`json` (referred to as :literal:`deconvertedSubjson`).
  - Before returning :literal:`returnObject`, check if :literal:`! originalSubjson.is_array( ) && deconvertedSubjson.is_array( )`. If this evals to :literal:`true`, then unidimensional array inference has taken place, and depending on the value of the :class:`json` object at the key path :literal:`#root.options.unidimensionalArrayInference`, a warning may be printed, an error may be thrown or execution may continue silently.

This means that, in the previous example, is unidimensional array inference is disabled, when creating the variable :literal:`childrenEnhancedAccess` using the :literal:`getValue` function, an :class:`IllegalValueError` will be thrown, but when creating :literal:`childrenBaiscAccess` using the :literal:`[]` operator or the :literal:`at` method, no error will be thrown and the user will not be informed about the fact that unidimensional array inference took place.

.. warning:: Unidimensional array inference is currently only implemented for :class:`std::vector`, or types that use the :class:`std::vector`'s :literal:`from_json` function in their :literal:`from_json` function, such as :class:`Eigen::Matrix`. In the future, if this feature is also wanted for other container types, such as :class:`std::set`, an overriden :literal:`from_json` function should be provided.

.. warning:: In order to check whether unidimensional array inference has taken place during the call to a :literal:`from_json` function, the :literal:`getValue` function converts the converted object of :literal:`ValueType` back to :class:`json` implicitly using the :literal:`to_json` function. This means that trying to use the :literal:`getValue` function for a type that does not have a :literal:`to_json` function will result in a compile error. Consequently, in the :literal:`json_interface`, all classes for which a :literal:`from_json` function is declared should also have a :literal:`to_json`. If the object is never going to be converted to :literal:`json`, this function could be left empty:

  .. code-block:: cpp

    void to_json( json& j, NonJsonableClass& nonJsonableObject ) { }

  leading to the generation of a :class:`json` object of value type :literal:`null`, for which the check on whether unidimensional array inference has been applied will always evaluate to :literal:`false`.


.. note:: Unidimensional array inference is widely used when working with :class:`Eigen::Vector`. An :class:`Eigen::Vector` is an :class:`Eigen::Matrix` with just one column. The JSON representation of a matrix is an array of arrays (with each array corresponding to a matrix row). Thus, the JSON representation of a row vector is an array of unidimensional arrays. For instance:

  .. code-block:: cpp

    Eigen::Vector3d zeroVector = Eigen::Vector3d::Zero( );
    std::cout << json( zeroVector ) << std::endl;

  yields:
    
  .. code-block:: json
  
    [
      [ 0 ],
      [ 0 ],
      [ 0 ]
    ]
  
  Thus, when the user provides e.g. the JSON array :literal:`[0, 0, 0]` and this is converted to an :literal:`Eigen::Vector3d`, unidimensional array inference is applied for each element, as an array of numbers is expected for each row but a number is found instead. If :literal:`options.unidimensionalArrayInference` is set to print warnings, this results messages of the type:
  
  .. code-block:: txt

    Unidimensional array inferred for key: keyWhereVectorIsStored.0
    Unidimensional array inferred for key: keyWhereVectorIsStored.1
    Unidimensional array inferred for key: keyWhereVectorIsStored.2

  This could be prevented by providing directly a row vector (in MATLAB, :literal:`rowVector = [0; 0; 0]`) instead of a column vector (in MATLAB :literal:`colVector = [0 0 0]`). However, the built-in MATLAB function :literal:`jsonencode` returns the same encoded JSON object for both :literal:`rowVector` and :literal:`colVector` (i.e. :literal:`[0, 0, 0]`). Thus, when using the JSON interface in combination with the MATLAB interface, unidimensional array inference will be applied frequently, since the vectors encoded by MATLAB are always column-vectors and Tudat expects row-vectors almost everywhere when using :class:`Eigen`.

