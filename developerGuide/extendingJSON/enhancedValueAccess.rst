.. _extendingJSON_enhancedValueAccess:

Enhanced value access
=====================

Once that the :ref:`main-json` file has been deserialized into the :literal:`mainJson` object using the :literal:`getDeserializedJSON` function, its keys can be accessed using the default value-access operators and functions, e.g.:

.. code-block:: cpp

    double initialEpoch = mainJson[ "initialEpoch" ];                       // 0.0
    std::string integratorType = mainJson.at( "integrator" ).at( "type" );  // "rungeKutta4"

The main problem of this approach is that, in case that the user didn't define the key :literal:`type` for the :literal:`integrator` object, the thrown value-access error will not provide sufficient information to uniquely identify the source of the problem:

.. code-block:: cpp

  libc++abi.dylib: terminating with uncaught exception of type nlohmann::detail::out_of_range:
  [json.exception.out_of_range.403] key 'type' not found

As it can be seen, for a simple input file such as :ref:`main-json` the key :literal:`type` is defined for six objects. If we forget to define it for any of them, the message :literal:`key 'type' not found` will be useless as we cannot know for which object it is missing. For more complex simulations, in which the JSON input file is hundreds of lines long, identification of the source of the problem becomes even more difficult. A message that *is* informative would be something like :literal:`key 'integrator.type' not found`. In the :literal:`json_interface`, :literal:`integrator.type` is known as a key path.


Key paths
~~~~~~~~~

A key path is a list of keys of a :class:`json` that are accessed sequentialy one after the other. The class :class:`KeyPath` is declared in the file :literal:`Tudat/External/JsonInterface/Support/keys.h`. This class derives from :class:`std::vector< std::string >` and has some additional features, such as the possibility of being initialized directly from a single :literal:`std::string` or being outputted as text:

.. code-block:: cpp

  KeyPath simpleKeyPath = "integrator";                   // { "integrator" }
  std::cout << simpleKeyPath << std::endl;                // integrator

  KeyPath compoundKeyPath = { "integrator", "type" };     // { "integrator", "type" }
  std::cout << compoundKeyPath << std::endl;              // integrator.type

Additionally, the operator :literal:`/` is overloaded for :class:`KeyPath` and :class:`std::string`:

.. code-block:: cpp

  KeyPath simpleKeyPath = "integrator";                     // { "integrator" }
  KeyPath compoundKeyPath = simpleKeyPath / "integrator";   // { "integrator", "type" }

Even the following is possible:

.. code-block:: cpp

  KeyPath compoundKeyPath = "integrator" / "type";          // { "integrator", "type" }


Now, image that the we want to access the values of the files to which to export the results. In our :ref:`main-json` example, using basic value access methods, this would be done by:

.. code-block:: cpp

  std::string file0 = mainJson.at( "export" ).at( 0 ).at( "file" );   // "epochs.txt"
  std::string file1 = mainJson.at( "export" ).at( 1 ).at( "file" );   // "states.txt"

However, a key path has been defined as a list of strings. If we want to define the key paths for those objects, we need to convert :literal:`0` and :literal:`1` to :literal:`"0"` and :literal:`"1"`. Thus:

.. code-block:: cpp

  KeyPath keyPathFile0 = "export" / 0 / "file";          // { "export", "0", "file" }
  KeyPath keyPathFile1 = "export" / 1 / "file";          // { "export", "1", "file" }
  std::cout << keyPathFile1 << std::endl;                // export.1.file

Note that the :literal:`/` operator is also overloaded for combinations of :class:`unsigned int` with :class:`string` or :class:`KeyPath` (but not for two pairs of :class:`unsigned int`), so :literal:`0` and :literal:`1` are implicitely converted to strings.

Although there is ambiguity when printing :literal:`keyPathFile1`, since :literal:`export.1` can be read as "the second element of an array defined at key :literal:`export`" or "the element at key :literal:`"1"` of an object defined at key :literal:`export`", there is no ambiguity when actually using the key path to access the values. If the :class:`json` object obtained when accessing the key :literal:`export` is of value type :literal:`object`, the key :literal:`"1"` will be used, while the integer value :literal:`1` will be used when the :class:`json` object at the key :literal:`export` is of value type :literal:`array`.


Error handling
~~~~~~~~~~~~~~

The templated function :literal:`ValueType getValue( const json& jsonObject, const KeyPath& keyPath )` declared in :class:`Tudat/External/JsonInterface/Support/valueAccess.h` returns the value of :literal:`jsonObject` defined at :literal:`keyPath` as a :class:`ValueType`. For instance:

.. code-block:: cpp

  getValue< std::string >( mainJson, "integrator" / "type" );   // "rungeKutta4"
  getValue< std::string >( mainJson, "export" / 1 / "file" );   // "states.txt"

In addition to recursively accessing the keys contained in :literal:`keyPath` and eventually transforming the last retrieved :class:`json` object to :class:`ValueType`, this function adds support for comprehensive value-access and value-conversion errors. For instance:

.. code-block:: cpp

  getValue< double >( mainJson, "integrator" / "errorTolerance" );

throws an :class:`UndefinedKeyError`:

.. code-block:: cpp

  libc++abi.dylib: terminating with uncaught exception of type tudat::json_interface::UndefinedKeyError:
  Undefined key: integrator.errorTolerance

And the following:

.. code-block:: cpp

  getValue< double >( mainJson, "export" / 1 / "file" );

throws an :class:`IllegalValueError`:

.. code-block:: cpp

  libc++abi.dylib: terminating with uncaught exception of type tudat::json_interface::IllegalValueError:
  Could not convert value to expected type double
  Illegal value for key: export.1.file


Now, image that we have an :class:`Integrator` class and we define its :literal:`from_json` function so that it can be created from :class:`json` objects:

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
  
  inline void from_json( const json& jsonIntegrator, Integrator& integrator )
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

Note that, in both cases, a :class:`json` object containing only the integrator object is passed to the :literal:`from_json` function. However, this object is not the same in both cases. When using the default basic value access, the following object is passed:

.. code-block:: json

  {
    "type": "rungeKutta4",
    "stepSize": 0
  }
  
When using the :literal:`getValue` function, the following :class:`json` object is passed:

.. code-block:: json

  {
    "type": "rungeKutta4",
    "stepSize": 10,
    "#keypath": [ "~", "integrator" ],
    "#root": {
      "initialEpoch": 0,
      "finalEpoch": 3600,
      "spice": {
        "kernels": [
          "${SPICE_KERNELS_PATH}/pck00009.tpc",
          "${SPICE_KERNELS_PATH}/de-403-masses.tpc",
          "${SPICE_KERNELS_PATH}/de421.bsp"
        ]
      },
      "bodies": {
        "Earth": {
          "useDefaultSettings": true,
          "ephemeris": {
            "constantState": [
              0,
              0,
              0,
              0,
              0,
              0
            ],
            "type": "constant"
          }
        },
        "asterix": {
          "initialState": {
            "semiMajorAxis": 7.5E+6,
            "eccentricity": 0.1,
            "inclination": 1.4888,
            "type": "keplerian"
          }
        }
      },
      "propagator": {
        "centralBodies": "Earth",
        "accelerations": {
          "asterix": {
            "Earth": {
              "type": "pointMassGravity"
            }
          }
        },
        "integratedStateType": "translational",
        "bodiesToPropagate": "asterix"
      },
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
        "unusedKey": "printWarning",
        "unidimensionalArrayInference": "continueSilently"
      }
    }
  }

I.e. the original :class:`mainObject` and the key path from which the integrator can be retrieved are also passed. Although the first two keys are not necessary (they contain redundant information), they are also included in the passed object in order to make it possible to use the default basic value access (i.e. the :literal:`[]` operator and the :literal:`at` method) on the passed :literal:`jsonIntegrator` object.


Special keys
~~~~~~~~~~~~

In order to make possible this advanced error handling in which the full key path is printed, a set of special keys are defined in :literal:`json_interface`. These special keys are subdivided in two categories:

  - Object-related. These special keys are assigned to :class:`json` objects by the :literal:`getValue` function. These keys must never be used in a JSON input file.

    - :literal:`#root` Key storing the contents of the root :class:`json` object.
    - :literal:`#keypath` Key storing the (absolute) key path from which a :class:`json` object was retrieved.

  - Path-related. These special keys are used only in :class:`KeyPath` objects.
    
    - :literal:`~` Known as root key. Used to denote that a key path is absolute (i.e. relative to the root :class:`json` object). Relative paths start with a key other than :literal:`~`.
    - :literal:`..` Known as up key. Used to navigate up one level in the key tree.
      
For instance, imagine that our :class:`Integrator` has an :literal:`initialTime` property. If we want this property to be retrieved from the :literal:`initialEpoch` key of the :literal:`mainJson` object in case it is not defined for the integrator, we can update its :literal:`from_json` like this:

.. code-block:: cpp
  
  inline void from_json( const json& jsonIntegrator, Integrator& integrator )
  {
      ...
      try
      {
          integrator.initialTime = getValue< double >( jsonIntegrator, "initialTime" );
      }
      catch ( UndefinedKeyError& e )
      {
          integrator.initialTime = getValue< double >( jsonIntegrator, ".." / "initialEpoch" );
      }
  }

In this case, we navigate one level up in the key by using the special key :literal:`..`, so we end up in the :literal:`mainJson` object, and then we *can* access the key :literal:`initialEpoch`. However, in case it is possible for integrator objects to be defined at different levels in the key tree (i.e. not always immediately under the :literal:`mainJson`), it is better to use an absolute key path:

.. code-block:: cpp
  
  integrator.initialTime = getValue< double >( jsonIntegrator, "~" / "initialEpoch" );


Multi-source properties
~~~~~~~~~~~~~~~~~~~~~~~

As illustrated in the previous section, some properties (referred to as multi-source properties), such as the integrator's initial time, can be retrieved from different key paths (:literal:`initialEpoch` or :literal:`integrator.initialTime`) This is so because the value at :literal:`initialEpoch` can also be used by other parts of the simulation (such as by Spice, to load the ephemeris from an initial time until a final time). In order to avoid repeating information, this value can be ommitted for the individual Spice and integrator objects and retrieved from the :literal:`mainJson` object instead. However, if Spice is not being used, it makes more sense to define inside the integrator object. Catching the :literal:`UndefinedKeyError`, as shown above, allows to use any of these values. Since this is done frequently in many parts of the :literal:`json_interface`, the :literal:`getValue` function has been overloaded to allow an :class:`std::vector< KeyPath >` as second argument. The value will be retrieved from the first defined key path in that list. If the none of the key paths in the list are defined, an :literal:`UndefinedKeyError` will be thrown printing the last key path in the key. Thus, the :literal:`try-catch` block shown above can be replaced by:

.. code-block:: cpp
  
  integrator.initialTime = getValue< double >(
    jsonIntegrator, { "initialTime", "~" / "initialEpoch" } );



Defaultable properties
~~~~~~~~~~~~~~~~~~~~~~

Imagine that the :literal:`type` property of our :class:`Integrator` class can take several values, e.g. :literal:`"euler"` and :literal:`"rungeKutta4"`. If we want the Runge-Kutta 4 integrator to be used when the key :literal:`type` is not defined by the user, we could modify the :literal:`from_json` method:

.. code-block:: cpp
  
  inline void from_json( const json& jsonIntegrator, Integrator& integrator )
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

This is also done frequently throughout the :literal:`json_interface`, so an overload for the :literal:`getValue` function with a third argument (the default value) is provided. Since the template type is inferred by the compiler from the third argument's type, in most cases we can safely remove the template arguments :literal:`< std::string >` from the function call. Thus, the previous :literal:`try-catch` block becomes:

.. code-block:: cpp
  
  integrator.type = getValue( jsonIntegrator, "type", "rungeKutta4" );

Note that, if the user *does* provide a value for the integrator's type, but it is not of the expected type (i.e. not a string), the default value will not be used, and an :literal:`IllegalValueError` will be thrown.


Other value-access functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In this subsection, a few functions widely used in the :literal:`json_interface`, all defined in :class:`Tudat/External/JsonInterface/Support/valueAccess.h`, are mentioned together with an example. For more information, see the Doxygen documentation [LINK].

- :literal:`bool isDefined( const json& jsonObject, const KeyPath& keyPath )`

  .. code-block:: cpp
    
    isDefined( mainJson, "integrator" / "initialTime" );   // false

- :literal:`ValueType getAs( const json& jsonObject )`

  .. code-block:: cpp
    
    Integrator integrator = getAs< Integrator >( jsonIntegrator );

- :literal:`json getRootObject( const json& jsonObject )`

  .. code-block:: cpp
    
    getRootObject( jsonIntegrator );   // mainJson

- :literal:`KeyPath getKeyPath( const json& jsonObject )`

  .. code-block:: cpp
    
    getKeyPath( jsonIntegrator );   // { "~", "integrator" }

- :literal:`std::string getParentKey( const json& jsonObject )`

  .. code-block:: cpp
    
    getParentKey( jsonIntegrator );   // "integrator"


Enhanced value access for arrays
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Imagine that our propagation requires the use of three different integrators for different periods of time. The user provides settings for the different integrators by defining the :literal:`integrator` key of :literal:`mainJson` to be an array with three elements, but forgets to define the third one:

.. code-block:: json

  [
    {
      "type": "rungeKutta4",
      "stepSize": 10,
    },
    {
      "type": "rungeKutta4",
      "stepSize": 20,
    }
  ]

Thus, when we try to access:

.. code-block:: cpp
  
  json integrators = mainJson.at( "integrator" );
  Integrator thirdIntegrator = integrators.at( 2 );

we get the following error:

.. code-block:: cpp
  
  libc++abi.dylib: terminating with uncaught exception of type nlohmann::detail::out_of_range:
  [json.exception.out_of_range.401] array index 2 is out of range

which is useless when trying to identify the source of the problem.

If we use enhanced value access:

.. code-block:: cpp
  
  json thirdIntegrator = getValue< json >( mainJson, "integrator" / 2 );

we do get a comprehensible error message:

.. code-block:: cpp
  
  libc++abi.dylib: terminating with uncaught exception of type tudat::json_interface::UndefinedKeyError:
  Undefined key: integrator.2

Although the functionality is identical for :class:`json` objects of value type :literal:`object` and :literal:`array`, the internal implementation is different and can have consequences for a developer extending the JSON interface. As explained in [REF], the comprehensible error messages are generated by defining the special key :literal:`#keypath` of the :class:`json` objects retrieved by using the :literal:`getValue` function. If the retrieved object is of value type :literal:`object`, defining this key is trivial. However, if the retrieved object is of value type :literal:`array`, the key cannot be defined directly, as string keys cannot be defined for :class:`json` arrays.

To overcome this problem, :class:`json` objects of value type :literal:`array` are converted first to :class:`json` objects of value type :literal:`object`. Then, the special keys :literal:`#keypath` and :literal:`#root` *can* be set. This means that, when calling :literal:`getValue< json >( ... )`, the returned :class:`json` will always be of value type :literal:`object` or :literal:`primitive`, but never :literal:`array`. If one wants to convert 


