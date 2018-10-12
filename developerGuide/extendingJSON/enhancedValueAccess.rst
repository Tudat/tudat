.. _extendingJSON_enhancedValueAccess:

.. role:: jsontype
.. role:: jsonkey

Enhanced value access
=====================

Once that the :ref:`main-json` file has been deserialized into the :literal:`mainJson` object using the :literal:`getDeserializedJSON` function, its keys can be accessed using the default value-access operators and functions, e.g.:

.. code-block:: cpp

    double initialEpoch = mainJson[ "initialEpoch" ];                       // 0.0
    std::string integratorType = mainJson.at( "integrator" ).at( "type" );  // "rungeKutta4"

The main problem of this approach is that, in case that the user didn't define the key :jsonkey:`type` for the :literal:`integrator` object, the thrown value-access error will not provide sufficient information to uniquely identify the source of the problem:

.. code-block:: txt

  libc++abi.dylib: terminating with uncaught exception of type nlohmann::detail::out_of_range:
  [json.exception.out_of_range.403] key 'type' not found

As it can be seen, for a simple input file such as :ref:`main-json` the key :jsonkey:`type` is defined for six objects. If we forget to define it for any of them, the message :literal:`key 'type' not found` will be useless as we cannot know for which object it is missing. The only way to find out is to determine which objects require a mandatory :literal:`type` key and to manually check the input file for any of those objects with this key missing. For more complex simulations, in which the JSON input file is hundreds of lines long, or splitted into several modular files, identification of the source of the problem can be difficult. A message that *is* informative would be something like :literal:`key 'integrator.type' not found`. In the JSON Interface library, :literal:`integrator.type` is known as a key path.

Key paths
~~~~~~~~~

A key path is a list of keys of a :class:`nlohmann::json` that are accessed sequentialy one after the other. The class :class:`KeyPath` is declared in the file :literal:`Tudat/JsonInterface/Support/keys.h`. This class derives from :literal:`std::vector< std::string >` and has some additional features, such as the capability to be initialised directly from a single :literal:`std::string` or being outputted as text:

.. code-block:: cpp

  KeyPath simpleKeyPath = "integrator";                   // { "integrator" }
  std::cout << simpleKeyPath << std::endl;                // integrator

  KeyPath compoundKeyPath = { "integrator", "type" };     // { "integrator", "type" }
  std::cout << compoundKeyPath << std::endl;              // integrator.type

  KeyPath compoundKeyPath2 = "integrator.type";           // { "integrator", "type" }
  std::cout << compoundKeyPath2 << std::endl;             // integrator.type

Additionally, the operator :literal:`/` is overloaded for :class:`KeyPath` and :class:`std::string`:

.. code-block:: cpp

  KeyPath simpleKeyPath = "integrator";                  // { "integrator" }
  KeyPath compoundKeyPath = simpleKeyPath / "type";      // { "integrator", "type" }


Now, image that the we want to access the values of the files to which to export the results. In our :ref:`main-json` example, using basic value access methods, this would be done by writing:

.. code-block:: cpp

  std::string file0 = mainJson.at( "export" ).at( 0 ).at( "file" );   // "epochs.txt"
  std::string file1 = mainJson.at( "export" ).at( 1 ).at( "file" );   // "states.txt"

However, a key path has been defined as a list of strings. If we want to define the key paths for those objects, we need to convert :literal:`0` and :literal:`1` to an unambiguous string representation. The special keys :literal:`@0`, :literal:`@1`, etc. are used, which means that the character :literal:`@` should not be used for the keys of the objects in the input files to avoid conflicts. Thus, we can write:

.. code-block:: cpp

  KeyPath keyPathFile0 = "export" / 0 / "file";          // { "export", "@0", "file" }
  std::cout << keyPathFile0 << std::endl;                // export[0].file
  KeyPath keyPathFile1 = "export[1].file";               // { "export", "@1", "file" }
  std::cout << keyPathFile1 << std::endl;                // export[1].file

Note that the :literal:`/` operator is also overloaded for combinations of :class:`unsigned int` with :class:`std::string` or :class:`KeyPath` (but not for two pairs of :class:`unsigned int`), so :literal:`0` and :literal:`1` are implicitly converted to strings.


Error handling
~~~~~~~~~~~~~~

The templated function :literal:`ValueType getValue( const nlohmann::json& jsonObject, const KeyPath& keyPath )` declared in :class:`Tudat/JsonInterface/Support/valueAccess.h` returns the value of :literal:`jsonObject` defined at :literal:`keyPath` as a :class:`ValueType`. For instance:

.. code-block:: cpp

  getValue< std::string >( mainJson, "integrator" / "type" );   // "rungeKutta4"
  getValue< std::string >( mainJson, "export[1].file" );        // "states.txt"

In addition to recursively accessing the keys contained in :literal:`keyPath` and eventually transforming the last retrieved :class:`nlohmann::json` object to :class:`ValueType`, this function adds support for comprehensive value-access and value-conversion errors. For instance:

.. code-block:: cpp

  getValue< double >( mainJson, "integrator" / "errorTolerance" );

throws an :class:`UndefinedKeyError`:

.. code-block:: txt

  libc++abi.dylib: terminating with uncaught exception of type tudat::json_interface::UndefinedKeyError:
  Undefined key: integrator.errorTolerance

And the following:

.. code-block:: cpp

  getValue< double >( mainJson, "export" / 1 / "file" );

throws an :class:`IllegalValueError`:

.. code-block:: txt

  libc++abi.dylib: terminating with uncaught exception of type tudat::json_interface::IllegalValueError:
  Could not convert value to expected type double
  Illegal value for key: export[1].file

since we are requesting to convert a :class:`nlohmann::json` object of value type :jsontype:`string` to :class:`double`.

Now, image that we have an :class:`Integrator` class and we define its :literal:`from_json` function so that it can be created from :class:`nlohmann::json` objects:

.. code-block:: cpp
  :caption: :class:`integrator.h`
  :name: integrator-h
  
  class Integrator
  {
  public:
      std::string type;
      double stepSize;
    
      Integrator( const std::string& type = "", const double stepSize = 0.0 ) :
        type( type ), stepSize( stepSize ) { }
  };
  
  inline void from_json( const nlohmann::json& jsonIntegrator, Integrator& integrator )
  {
      integrator.type = getValue< std::string >( jsonIntegrator, "type" );
      integrator.stepSize = getValue< double >( jsonIntegrator, "stepSize" );
  }

If, somewhere else, we write:

.. code-block:: cpp

  Integrator integrator = mainJson.at( "integrator" );

the default value access functions will be used, leading to error messages in case of missing keys or illegal values that are difficult to debug. Thus, the following should be used instead:

.. code-block:: cpp

  Integrator integrator = getValue< Integrator >( mainJson, "integrator" );

Note that, in both cases, a :class:`nlohmann::json` object containing only the integrator object is passed to the :literal:`from_json` function. However, this object is not the same in both cases. When using the default basic value access, the following object is passed:

.. code-block:: json

  {
    "type": "rungeKutta4",
    "stepSize": 0
  }
  
When using the :literal:`getValue` function, the following :class:`nlohmann::json` object is passed:

.. code-block:: json

  {
    "type": "rungeKutta4",
    "stepSize": 10,
    "#keypath": ["~", "integrator"],
    "#root": {
      "initialEpoch": 0,
      "finalEpoch": 3600,
      "spice": {
        "useStandardKernels": true
      },
      "bodies": {
        "Earth": {
          "useDefaultSettings": true,
          "ephemeris": {
            "type": "constant",
            "constantState": [0, 0, 0, 0, 0, 0]
          }
        },
        "asterix": {
          "initialState": {
            "type": "keplerian",
            "semiMajorAxis": 7.5E+6,
            "eccentricity": 0.1,
            "inclination": 1.4888
          }
        }
      },
      "propagators": [
        {
          "integratedStateType": "translational",
          "bodiesToPropagate": ["asterix"],
          "centralBodies": ["Earth"],
          "accelerations": {
            "asterix": {
              "Earth": [
                {
                  "type": "pointMassGravity"
                }
              ]
            }
          }
        }
      ],
      "integrator": {
        "type": "rungeKutta4",
        "stepSize": 10
      },
      "export": [
        {
          "file": "@path(epochs.txt)",
          "variables": [
            {
              "type": "independent"
            }
          ],
          "epochsInFirstColumn": false
        },
        {
          "file": "@path(states.txt)",
          "variables": [
            {
              "type": "state"
            }
          ],
          "epochsInFirstColumn": false
        }
      ],
      "options": {
        "defaultValueUsedForMissingKey": "continueSilently",
        "unusedKey": "printWarning"
      }
    }
  }

I.e. the original :class:`mainObject` and the key path from which the integrator can be retrieved are also passed. Although the first two keys are not necessary (they contain redundant information that could be retrieved from :literal:`#root` at :literal:`#keypath`), they are also included in the returned object so that it can be used with the default basic value access (i.e. the :literal:`[]` operator and the :literal:`at` method).


Special keys
~~~~~~~~~~~~

In order to make possible this advanced error handling in which the full key path is printed, a set of special keys are defined in the JSON Interface library. These special keys are subdivided into two categories:

  - Object-related. These special keys are assigned to :class:`nlohmann::json` objects by the :literal:`getValue` function. These keys must never be used in a JSON input file.

    - :literal:`#root`: stored the contents of the root :class:`nlohmann::json` object.
    - :literal:`#keypath`: stored the (absolute) key path from which a :class:`nlohmann::json` object is retrieved.

  - Path-related. These special keys are used only in the elements of :class:`KeyPath` objects.
    
    - :literal:`~`: known as root key. Used to denote that a key path is absolute (i.e. relative to the root :class:`nlohmann::json` object). Relative paths start with a key other than :literal:`~`.
    - :literal:`<-`: known as up key. Used to navigate up one level in the key tree.
      
For instance, imagine that our :class:`Integrator` has an :literal:`initialTime` property. If we want this property to be retrieved from the :literal:`initialEpoch` key of the :literal:`mainJson` object in case it is not defined for the integrator, we can update its :literal:`from_json` to look like this:

.. code-block:: cpp
  
  inline void from_json( const nlohmann::json& jsonIntegrator, Integrator& integrator )
  {
      ...
      try
      {
          integrator.initialTime = getValue< double >( jsonIntegrator, "initialTime" );
      }
      catch ( UndefinedKeyError& e )
      {
          integrator.initialTime = getValue< double >( jsonIntegrator, "<-" / "initialEpoch" );
      }
  }

In this case, we navigate one level up in the key by using the special key :jsonkey:`<-`, so we end up in the :literal:`mainJson` object, and then we *can* access the key :jsonkey:`initialEpoch`. However, in case it is possible for integrator objects to be defined at different levels in the key tree (i.e. not always defined in one of the keys immediately under the :literal:`mainJson`), it is better to use an absolute key path:

.. code-block:: cpp
  
  integrator.initialTime = getValue< double >( jsonIntegrator, "~" / "initialEpoch" );

.. note:: When printing a :class:`KeyPath`, either during the generation of an error or for debugging, its canonical representation is used. Canonical key paths are always absolute (i.e. relative to the :literal:`mainJson`), so there is no need to print the initial root key (:literal:`~`) as there is no possible ambiguity. Additionally, the up keys (:literal:`<-`) are removed, popping back the previous key. For instance:

  .. code-block:: cpp
    
    std::cout << "integrator" / "type" << std::endl;                   // integrator.type
    std::cout << "~" / "integrator" / "type" << std::endl;             // integrator.type
    std::cout << "integrator" / "<-" / "initialEpoch" << std::endl;    // initialEpoch

.. note:: The :literal:`from_json` function of :literal:`std::map` and :literal:`std::unordered_map` have been overridden in :class:`Tudat/JsonInterface/Support/valueConversions.h`, so that the special keys :jsonkey:`#root` and :jsonkey:`#keypath` are not assigned to the converted map.


Multi-source properties
~~~~~~~~~~~~~~~~~~~~~~~

As illustrated in the previous section, some properties (referred to as multi-source properties), such as the integrator's initial time, can be retrieved from different key paths (:literal:`initialEpoch` or :literal:`integrator.initialTime`) This is so because the value at :literal:`initialEpoch` can also be used by other parts of the simulation (such as by Spice, to load the ephemeris from an initial time until a final time). In order to avoid repeating information, this value can be omitted for the individual Spice and integrator objects and retrieved from the :literal:`mainJson` object instead. However, if Spice is not being used, it makes more sense to define it inside the integrator object. Catching the :literal:`UndefinedKeyError`, as shown above, allows to use any of these values. Since this is done frequently in many parts of the :literal:`json_interface`, the :literal:`getValue` function has been overloaded to allow an :literal:`std::vector< KeyPath >` as second argument. The value will be retrieved from the first defined key path in that list. If none of the key paths in the list are defined, an :literal:`UndefinedKeyError` will be thrown printing the last key path in the vector. Thus, the :literal:`try-catch` block shown above can be replaced by:

.. code-block:: cpp
  
  integrator.initialTime = getValue< double >(
    jsonIntegrator, { "initialTime", "~" / "initialEpoch" } );



Defaultable properties
~~~~~~~~~~~~~~~~~~~~~~

Imagine that the :literal:`type` property of our :class:`Integrator` class can take several values, e.g. :literal:`"euler"` and :literal:`"rungeKutta4"`. If we want the Runge-Kutta 4 integrator to be used when the key :jsonkey:`type` is not defined by the user, we could modify the :literal:`from_json` method:

.. code-block:: cpp
  
  inline void from_json( const nlohmann::json& jsonIntegrator, Integrator& integrator )
  {
      try
      {
          integrator.type = getValue< std::string >( jsonIntegrator, "type" );
      }
      catch ( UndefinedKeyError& e )
      {
          integrator.type = "rungeKutta4";
      }
      ...
  }

This is also done frequently throughout the JSON Interface library, so an overload for the :literal:`getValue` function with a third argument (the default value) is provided. Since the template type is inferred by the compiler from the third argument's type, in most cases we can safely remove the template arguments (e.g. :literal:`< std::string >`) from the function call. Thus, the previous :literal:`try-catch` block becomes:

.. code-block:: cpp
  
  integrator.type = getValue( jsonIntegrator, "type", "rungeKutta4" );

Note that, if the user *does* provide a value for the integrator's type, but it is not of the expected type (i.e. not a string), the default value will not be used, and an :literal:`IllegalValueError` will be thrown.


Other value-access functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In this subsection, a few functions widely used in the :literal:`json_interface`, all defined in :class:`Tudat/JsonInterface/Support/valueAccess.h`, are mentioned together with an example. For more information, see the Doxygen documentation.

- :literal:`bool isDefined( const nlohmann::json& jsonObject, const KeyPath& keyPath )`

  .. code-block:: cpp
    
    isDefined( mainJson, "integrator" / "initialTime" );   // false

- :literal:`ValueType getAs( const nlohmann::json& jsonObject )`

  .. code-block:: cpp
    
    Integrator integrator = getAs< Integrator >( jsonIntegrator );

- :literal:`nlohmann::json getRootObject( const nlohmann::json& jsonObject )`

  .. code-block:: cpp
    
    getRootObject( jsonIntegrator );   // mainJson

- :literal:`KeyPath getKeyPath( const nlohmann::json& jsonObject )`

  .. code-block:: cpp
    
    getKeyPath( jsonIntegrator );   // { "~", "integrator" }

- :literal:`std::string getParentKey( const nlohmann::json& jsonObject )`

  .. code-block:: cpp
    
    getParentKey( jsonIntegrator );   // "integrator"


Enhanced value access for arrays
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Imagine that our propagation requires the use of three different integrators for different periods of time. The user provides settings for the different integrators by defining the :literal:`integrator` key of :literal:`mainJson` to be an array with three elements:

.. code-block:: json

  [
    {
      "type": "rungeKutta4",
      "stepSize": 10,
    },
    {
      "type": "rungeKutta4",
      "stepSize": 20,
    },
    {
      "type": "rungeKutta4"
    }
  ]

Thus, when we try to access:

.. code-block:: cpp
  
  nlohmann::json integrators = mainJson.at( "integrator" );
  double thirdIntegratorStepSize = integrators.at( 2 ).at( "stepSize" );

we get the following error:

.. code-block:: txt
  
  libc++abi.dylib: terminating with uncaught exception of type nlohmann::detail::out_of_range:
  [json.exception.out_of_range.403] key 'stepSize' not found

which is not very informative when trying to identify the source of the problem.

If we use enhanced value access:

.. code-block:: cpp
  
  double thirdIntegratorStepSize = getValue< nlohmann::json >( mainJson, "integrator[2].stepSize" );

the used does get a comprehensible error message:

.. code-block:: txt
  
  libc++abi.dylib: terminating with uncaught exception of type tudat::json_interface::UndefinedKeyError:
  Undefined key: integrator[2].stepSize

Although the functionality is identical for :class:`nlohmann::json` objects of value type :jsontype:`object` and :jsontype:`array`, the internal implementation is different and can have consequences for a developer extending the JSON interface. As explained before, the comprehensible error messages are generated by defining the special key :jsonkey:`#keypath` of the :class:`nlohmann::json` objects retrieved by using the :literal:`getValue` function. If the retrieved object is of value type :jsontype:`object`, defining this key is trivial. However, if the retrieved object is of value type :jsontype:`array`, the key cannot be defined directly, as string keys cannot be defined for :class:`nlohmann::json` arrays.

To overcome this problem, :class:`nlohmann::json` objects of value type :jsontype:`array` containing structured objects are converted first to :class:`nlohmann::json` objects of value type :jsontype:`object`. Only then, it *is* possible to set the special keys :jsonkey:`#keypath` and :jsonkey:`#root`. This means that, when calling :literal:`getValue< nlohmann::json >( ... )`, the returned :class:`nlohmann::json` may be of value type :jsontype:`object` even though the original objects was of value type :jsontype:`array`. For :jsontype:`primitive` types and arrays containing only (arrays of) primitive types, there is no need to convert them to :jsontype:`object` and define the special keys.

The process of converting :class:`nlohmann::json` objects from value type :jsontype:`array` to value type :jsontype:`object` is done inside the :literal:`getValue` function automatically. For instance:

.. code-block:: cpp
  
  nlohmann::json integrators = getValue< nlohmann::json >( mainJson, "integrator" );

generate the following :class:`nlohmann::json` object:

.. code-block:: json

  {
    "@0": {
      "type": "rungeKutta4",
      "stepSize": 10
    },
    "@1": {
      "type": "rungeKutta4",
      "stepSize": 20
    },
    "@2": {
      "type": "rungeKutta4"
    },
    "#keypath": ["~", "integrator"],
    "#root": {}
  }

where the root object, i.e. :literal:`mainJson`, has been omitted in this code-block. If one wants to check whether the returned object actually represents an array, instead of using the built-in :literal:`is_array` method, one has to use the function :literal:`bool isConvertibleToArray( const nlohmann::json& j )` defined in :class:`Tudat/JsonInterface/Support/valueAccess.h`. This function returns :literal:`true` if :literal:`j` is of value type :jsontype:`object` and all its non-special keys match the expression :literal:`@i`, with :literal:`i` convertible to :class:`int`, or if :literal:`j` is already of value type :jsontype:`array`. After this check, it is safe to call the function :literal:`nlohmann::json getAsArray( const nlohmann::json& jsonObject )` to convert the object back to array type. During this process, the information stored in the special keys is lost, so this is rarely done. Instead, the :literal:`from_json` function of :class:`std::vector` has been overridden so that it is possible to write:

.. code-block:: cpp
  
  nlohmann::json jsonObject = getValue< nlohmann::json >( mainJson, "integrator" );
  jsonObject.is_array( );                    // false
  isConvertibleToArray( jsonObject );        // true
  std::vector< Integrator > integrators = getAs< std::vector< Integrator > >( jsonObject );

.. warning:: No custom implementation of the :literal:`from_json` function for :literal:`std::set` is provided by :literal:`json_interface`, since this type is not used by Tudat (as of now). In the future, if one wants to use the :literal:`getValue` function with :literal:`std::set` as template argument, the default :literal:`from_json` function for :literal:`std::set` will have to be overridden to allow conversion of :class:`nlohmann::json` objects of value type :jsontype:`object` to :literal:`std::set`, in a similar way as has been done for :literal:`std::vector` in :class:`Tudat/JsonInterface/Support/valueConversions.h`.
