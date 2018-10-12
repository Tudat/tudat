.. _extendingJSON_settingsClasses:

.. role:: jsontype
.. role:: jsonkey

Settings classes
================

Almost all classes containing settings for Tudat simulations are used in the form of shared pointers. Currently, :class:`std::shared_ptr` is being used. Support for converting shared pointers to and from :class:`nlohmann::json` could be added just by writing:

.. code-block:: cpp
  
  namespace std
  {
  
  template< typename T >
  void to_json( nlohmann::json& jsonObject, const shared_ptr< T >& sharedPointer )
  {
      if ( sharedPointer != NULL )
      {
          jsonObject = *sharedPointer;                    // T's to_json called
      }
  }

  template< typename T >
  void from_json( const nlohmann::json& jsonObject, shared_ptr< T >& sharedPointer )
  {
      if ( ! jsonObject.is_null( ) )
      {
          T object = jsonObject;                          // T's from_json called
          sharedPointer = make_shared< T >( object );
      }
  }
  
  }

Then, if class :class:`T` is convertible to/from :class:`nlohmann::json`, :literal:`std::shared_ptr< T >` would also be. However, this approach has some drawbacks, namely:

  - :class:`T` must be default-constructible for the code to compile.
  - Separate :literal:`to_json` and :literal:`from_json` functions have to be written for each derived class of :class:`T`, potentially leading to the duplication of code.
  
Thus, given that the settings classes used throughout Tudat are always used as shared pointers, rather than providing :literal:`to_json` and :literal:`from_json` functions for these classes, these functions have been written for shared pointers of these classes. Before introducing the best-practices to be followed when writing these functions, the way in which the keys to be used in these functions are defined in the JSON Interface library is described.

Definition of keys
~~~~~~~~~~~~~~~~~~

In all the example :literal:`to_json` and :literal:`from_json` functions presented so far, the keys were hard-coded, i.e. literal strings were used when using the :literal:`[ ]` operator of a :class:`nlohmann::json` object, calling the :literal:`getValue` function or constructing a key path (by concatenating several strings). However, this approach makes code-updating very complex. Image that, in the future, we want to update a key called :literal:`initialTime` to :literal:`initialEpoch`. Although a global search could do the trick, this may result in modifying parts of the code that should not be modified. If we want to update the name of the key :jsonkey:`type` to :literal:`modelType`, but only for rotation model settings, the only option is doing it manually to avoid changing also the :literal:`type` keys of other objects.

Thus, in the JSON Interface library, literal strings are never used inside :literal:`to_json` and :literal:`from_json` functions. Instead, all the keys that are recognised by the JSON Interface are declared in :class:`Tudat/JsonInterface/Support/keys.h`, and their string-value is defined in :class:`Tudat/JsonInterface/Support/keys.cpp`. This is done using a struct called :class:`Keys` containing several nested structs for each level. For instance:

.. code-block:: cpp
  :caption: :class:`Tudat/JsonInterface/Support/keys.h`
  :name: keys-h
  
  namespace tudat
  {

  namespace json_interface
  {

  struct Keys
  {
      static const std::string initialEpoch;
      static const std::string finalEpoch;
      
      ...

      static const std::string bodies;
      struct Body
      {
          ...

          static const std::string rotationModel;
          struct RotationModel
          {
              static const std::string type;
              static const std::string originalFrame;
              static const std::string targetFrame;
              static const std::string initialOrientation;
              static const std::string initialTime;
              static const std::string rotationRate;
          };
          
          ...
      };
      
      ...
  };
  
  }  // namespace json_interface
  
  }  // namespace tudat


.. code-block:: cpp
  :caption: :class:`Tudat/JsonInterface/Support/keys.cpp`
  :name: keys-cpp
  
  namespace tudat
  {

  namespace json_interface
  {

  const std::string Keys::initialEpoch = "initialEpoch";
  const std::string Keys::finalEpoch = "finalEpoch";
  
  ...
  
  //  Body
  const std::string Keys::bodies = "bodies";
  
  ...
  
  // //  Body::RotationModel
  const std::string Keys::Body::rotationModel = "rotationModel";
  const std::string Keys::Body::RotationModel::type = "type";
  const std::string Keys::Body::RotationModel::originalFrame = "originalFrame";
  const std::string Keys::Body::RotationModel::targetFrame = "targetFrame";
  const std::string Keys::Body::RotationModel::initialOrientation = "initialOrientation";
  const std::string Keys::Body::RotationModel::initialTime = "initialTime";
  const std::string Keys::Body::RotationModel::rotationRate = "rotationRate";

  ...
    
  }  // namespace json_interface
  
  }  // namespace tudat

Note that the keys for the different derived classes of :class:`RotationModelSettings` are all defined at the same level (i.e. a different struct is not created for each derived class). When going through a settings class and defining its keys, it is good practice to define also the keys for the derived classes that will not supported by the JSON Interface (initially), and commenting them out.

When one wants to modify a key, changing its string value in :class:`keys.cpp` should suffice. However, it is good practice to keep the name of the keys and the values of the keys consistent, so "Rename Symbol Under Cursor" should be used as well to replace all the occurrences of the key.

.. caution:: When debugging an :class:`UndefinedKeyError`, the following situation can arise when parsing, for instance, the following JSON file (only relevant section shown):

  .. code-block:: json
  
    "bodies": {
      "Earth": {
        "rotationModel": {
          "type": "simple",
          "originalFrame": "ECLIPJ2000",
          "targetFrame": "IAU_Earth",
          "initialTime": 0,
          "rotationRate": 7e-05
        }
      }
    }

  .. code-block:: txt

    libc++abi.dylib: terminating with uncaught exception of type tudat::json_interface::UndefinedKeyError:
    Undefined key: bodies.Earth.rotationModel.initialOrientation

  The key :jsonkey:`initialOrientation` stores a defaultable-property, i.e. if not provided the value can be inferred from the keys :literal:`originalFrame`, :literal:`targetFrame` and :literal:`initialTime`. Thus, the error should not be generated if the key is not defined. One would probably tend to start by looking at the :literal:`from_json` function of :class:`EphemerisSettings` when debugging this issue. However, the source of the problem can be in the :class:`keys.cpp`. Even if the :literal:`from_json` function is completely correct, the previous error would be printed if, in the :class:`keys.cpp`, we had:

  .. code-block:: cpp
    
    const std::string Keys::Body::RotationModel::initialOrientation = "initialOrientation";
    const std::string Keys::Body::RotationModel::initialTime = "initialOrientation";

  which can happen easily when copy-pasting. Thus, what is actually happening is that, when retrieving the value for :literal:`initialTime` (non-defaultable property), the key :jsonkey:`initialOrientation` is mistakenly accessed (and not found). To prevent these issues, a search for any given key inside the :class:`keys.cpp` file should always result in an even number of occurrences. In this way, we also make sure that e.g. the value stored at the key :jsonkey:`rotationRate` does not end up being used for the property :literal:`initialTime` of our :class:`RotationModelSettings`, in which case no error or warning would be generated during conversion to :class:`nlohmann::json` as both as non-defaultable properties and store values of the same type (:class:`double`).
  

Writing :literal:`from_json` functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Generally, the settings classes used in Tudat are not default-constructible. By providing :literal:`from_json` functions for their shared pointers, rather than for the class itself, we can get the code to compile without the need to provide default constructors because a shared pointer is default-constructible (it is :literal:`NULL` by default).

To illustrate the structure of a :literal:`from_json` function for a shared pointer to a settings class, the :class:`RotationModelSettings` example is described here:

.. code-block:: cpp
  :linenos:
  :caption: :class:`Tudat/JsonInterface/Environment/rotationModel.h`
  :name: rotationModel-h
  
  #include <Tudat/SimulationSetup/EnvironmentSetup/createRotationModel.h>

  #include "Tudat/JsonInterface/Support/valueAccess.h"
  #include "Tudat/JsonInterface/Support/valueConversions.h"

  namespace tudat
  {

  namespace simulation_setup
  {

  //! Map of `RotationModelType`s string representations.
  static std::map< RotationModelType, std::string > rotationModelTypes =
  {
      { simple_rotation_model, "simple" },
      { spice_rotation_model, "spice" }
  };

  //! `RotationModelType`s not supported by `json_interface`.
  static std::vector< RotationModelType > unsupportedRotationModelTypes = { };

  //! Convert `RotationModelType` to `nlohmann::json`.
  inline void to_json( nlohmann::json& jsonObject, const RotationModelType& rotationModelType )
  {
      jsonObject = json_interface::stringFromEnum( rotationModelType, rotationModelTypes );
  }

  //! Convert `nlohmann::json` to `RotationModelType`.
  inline void from_json( const nlohmann::json& jsonObject, RotationModelType& rotationModelType )
  {
      rotationModelType = json_interface::enumFromString( jsonObject, rotationModelTypes );
  }

  //! Create a `nlohmann::json` object from a shared pointer to a `RotationModelSettings` object.
  void to_json( nlohmann::json& jsonObject, const std::shared_ptr< RotationModelSettings >& rotationModelSettings );

  //! Create a shared pointer to a `RotationModelSettings` object from a `nlohmann::json` object.
  void from_json( const nlohmann::json& jsonObject, std::shared_ptr< RotationModelSettings >& rotationModelSettings );

  } // namespace simulation_setup

  } // namespace tudat

.. code-block:: cpp
  :linenos:
  :caption: :class:`Tudat/JsonInterface/Environment/rotationModel.cpp`
  :name: rotationModel-cpp-from-json
  
  namespace tudat
  {

  namespace simulation_setup
  {

  ...

  //! Create a shared pointer to a `RotationModelSettings` object from a `nlohmann::json` object.
  void from_json( const nlohmann::json& jsonObject, std::shared_ptr< RotationModelSettings >& rotationModelSettings )
  {
      using namespace json_interface;
      using K = Keys::Body::RotationModel;

      // Base class settings
      const RotationModelType rotationModelType = getValue< RotationModelType >( jsonObject, K::type );
      const std::string originalFrame = getValue< std::string >( jsonObject, K::originalFrame );
      const std::string targetFrame = getValue< std::string >( jsonObject, K::targetFrame );

      switch ( rotationModelType ) {
      case simple_rotation_model:
      {
          const double initialTime = getValue< double >( jsonObject, K::initialTime );

          // Get JSON object for initialOrientation (or create it if not defined)
          nlohmann::json jsonInitialOrientation;
          if ( isDefined( jsonObject, K::initialOrientation ) )
          {
              jsonInitialOrientation = getValue< nlohmann::json >( jsonObject, K::initialOrientation );
          }
          else
          {
              jsonInitialOrientation[ K::originalFrame ] = originalFrame;
              jsonInitialOrientation[ K::targetFrame ] = targetFrame;
              jsonInitialOrientation[ K::initialTime ] = initialTime;
          }

          rotationModelSettings = std::make_shared< SimpleRotationModelSettings >(
                      originalFrame,
                      targetFrame,
                      getAs< Eigen::Quaterniond >( jsonInitialOrientation ),
                      initialTime,
                      getValue< double >( jsonObject, K::rotationRate ) );
          return;
      }
      case spice_rotation_model:
      {
          rotationModelSettings = std::make_shared< RotationModelSettings >(
                      rotationModelType, originalFrame, targetFrame );
          return;
      }
      default:
          handleUnimplementedEnumValue( rotationModelType, rotationModelTypes, unsupportedRotationModelTypes );
      }
  }
  
  } // namespace simulation_setup

  } // namespace tudat

Note that the files :class:`Tudat/JsonInterface/Support/valueAccess.h` and :class:`Tudat/JsonInterface/Support/valueConversions.h` are always included. The former includes enhaced value access functions (:literal:`getValue`) and the latter overrides (and defines) :literal:`to_json` and :literal:`from_json` functions for frequently-used types, such as :class:`std::vector` or :class:`Eigen::Matrix`. The file :class:`Tudat/JsonInterface/Support/valueAccess.h` includes :class:`Tudat/JsonInterface/Support/keys.h`, so all the keys available in the JSON Interface are readily accessible.

Typically, the first lines of a :literal:`from_json` function are (for our example):

.. code-block:: cpp

  using namespace json_interface;
  using K = Keys::Body::RotationModel;

Generally, when creating a shared pointer to a settings class, only the keys for that class are needed (unless some keys of the root :literal:`mainJson` object have to be accessed). Thus we can use a shorter-name such as :literal:`K`. The other keys can still be accessed using the full name :literal:`Keys::...`. Do not write :literal:`using Keys = Keys::...`, as this would result in all the other keys being unaccessible.

Although it may be convenient to make the check on whether the provided :literal:`jsonObject` object is :literal:`null`, and return immediately a :literal:`NULL` shared pointer if it is, this situation will generally not happen in practice. When the user does not want to provide a rotation model, rather than writing :literal:`"rotationModel": null`, they leave the key :jsonkey:`rotationModel` undefined. The :literal:`from_json` function of :class:`BodySettings` is responsible for only calling the :literal:`from_json` function of :class:`RotationModelSettings` if the key :jsonkey:`rotationModel` is defined. If the user does provide :literal:`null` manually in their input file, this will result in an :literal:`UndefinedKeyError` for key :jsonkey:`bodies.Earth.rotationModel.type` (the first key to be accessed in the :literal:`from_json` function) will be thrown.

The settings classes used in Tudat typically have a type property that can be used to determine which derived class should be used when creating the shared pointer object. This is retrieved in line 15. Then, the settings for the base class (shared by all the derived classes) are retrieved in lines 16 and 17. The next step is to write a :class:`switch` that modifies (re-constructs) the :literal:`rotationModelSettings` (passed by reference) depending on the :literal:`rotationModelType`. Generally, a :literal:`return` is added at the end of each switch case.

Finally, the default case of the switch always calls the :literal:`handleUnimplementedEnumValue` function, with the first argument the type converted from JSON, the second argument the map of string representations for the enumeration, and the third argument the list of enumeration values not supported by the JSON interface. This will throw an :literal:`UnsupportedEnumError`, suggesting the user to write their own JSON-based C++ application, or an :literal:`UnimplementedEnumError`, if the provided enum value is not marked as unsupported, but we (the coders) forgot to write its implementation, printing a warning in which the user is kindly asked to open an issue on GitHub.

In most cases, the defaultable properties use the default value defined in the setting class constructors. For instance, consider the following constructor:

.. code-block:: cpp

  BasicSolidBodyGravityFieldVariationSettings(
          const std::vector< std::string > deformingBodies,
          const std::vector< std::vector< std::complex< double > > > loveNumbers,
          const double bodyReferenceRadius,
          const std::shared_ptr< InterpolatorSettings > interpolatorSettings = NULL ):

If the user does not provide the key :jsonkey:`interpolator`, the same default value defined in the constructor (:literal:`NULL`) should be used to create the settings object from the :class:`nlohmann::json` object. To keep the behaviour of C++ Tudat applications and JSON-based Tudat applications consistent, if in the future the default interpolator settings are changed from :literal:`NULL` to e.g. :literal:`std::make_shared< LagrangeInterpolatorSettings >( 6 )` in the constructor, this change should also be reflected in the JSON Interface. To make this happen automatically, the default values are not hard-coded in the :literal:`from_json` functions. Instead, an instance constructed only with the mandatory properties is used to create the actual shared pointer:

.. code-block:: cpp

  BasicSolidBodyGravityFieldVariationSettings defaults( { }, { }, TUDAT_NAN );
  gravityFieldVariationSettings = std::make_shared< BasicSolidBodyGravityFieldVariationSettings >(
              getValue< std::vector< std::string > >( jsonObject, K::deformingBodies ),
              getValue< std::vector< std::vector< std::complex< double > > > >( jsonObject, K::loveNumbers ),
              getValue< double >( jsonObject, K::referenceRadius ),
              getValue( jsonObject, K::interpolator, defaults.getInterpolatorSettings( ) ) );

Since the settings classes are only used to store information, their constructors are generally empty, so in most cases the temporary object from which the default properties is retrieved can be constructed with values such as :literal:`{ }`, :literal:`""` or :literal:`TUDAT_NAN` for the mandatory properties.



Writing :literal:`to_json` functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

An example of a :literal:`to_json` function is provided below:

.. code-block:: cpp
  :linenos:
  :caption: :class:`Tudat/JsonInterface/Environment/rotationModel.cpp`
  :name: rotationModel-cpp-to-json
  
  namespace tudat
  {

  namespace simulation_setup
  {

  //! Create a `nlohmann::json` object from a shared pointer to a `RotationModelSettings` object.
  void to_json( nlohmann::json& jsonObject, const std::shared_ptr< RotationModelSettings >& rotationModelSettings )
  {
      if ( ! rotationModelSettings )
      {
          return;
      }
      using namespace json_interface;
      using K = Keys::Body::RotationModel;

      const RotationModelType rotationModelType = rotationModelSettings->getRotationType( );
      jsonObject[ K::type ] = rotationModelType;
      jsonObject[ K::originalFrame ] = rotationModelSettings->getOriginalFrame( );
      jsonObject[ K::targetFrame ] = rotationModelSettings->getTargetFrame( );

      switch ( rotationModelType )
      {
      case simple_rotation_model:
      {
          std::shared_ptr< SimpleRotationModelSettings > simpleRotationModelSettings =
                  std::dynamic_pointer_cast< SimpleRotationModelSettings >( rotationModelSettings );
          assertNonNullPointer( simpleRotationModelSettings );
          jsonObject[ K::initialOrientation ] = simpleRotationModelSettings->getInitialOrientation( );
          jsonObject[ K::initialTime ] = simpleRotationModelSettings->getInitialTime( );
          jsonObject[ K::rotationRate ] = simpleRotationModelSettings->getRotationRate( );
          return;
      }
      case spice_rotation_model:
          return;
      default:
          handleUnimplementedEnumValue( rotationModelType, rotationModelTypes, unsupportedRotationModelTypes );
      }
  }
  
  ...
  
  }  // namespace simulation_setup
  
  }  // namespace tudat

First, a check on the nullity of the shared pointer is done. If it is :literal:`NULL`, the null :class:`nlohmann::json` object won't be modified. Otherwise, the settings object will be used to define the keys of the :class:`nlohmann::json` object.

The structure is similar to the one for :literal:`from_json` functions. The main difference is that the :literal:`[]` mutator operator is used to modify the :class:`nlohmann::json` object, instead of using the :literal:`getValue` function to access it. Additionally, in every switch case the original shared pointer has to be dynamically casted to the corresponding derived class. Then, the function :literal:`assertNonNullPointer` is called. This throws an :literal:`NullPointerError` when the settings derived class and the value of its type property do not match.

Some switch cases, such as :literal:`spice_rotation_model`, are empty because they do not contain additional information other than that of the original base class. Thus, the only needed statement is :literal:`return;`.

