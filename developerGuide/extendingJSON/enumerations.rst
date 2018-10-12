.. _extendingJSON_enumerations:

.. role:: jsontype
.. role:: jsonkey

Enumerations
============

Enumerations are widely used by Tudat to define the type of settings objects. In this section, it will be explained how to implement support for the enumeration :literal:`EphemerisType`. The same process can be followed to add support for other enumerations in the future.

.. code-block:: cpp
  :caption: :class:`Tudat/SimulationSetup/EnvironmentSetup/createEphemeris.h`
  :name: createEphemeris-h
  
  ...
  
  namespace tudat
  {

  namespace simulation_setup
  {

  ...
  enum EphemerisType
  {
      approximate_planet_positions,
      direct_spice_ephemeris,
      tabulated_ephemeris,
      interpolated_spice,
      constant_ephemeris,
      kepler_ephemeris,
      custom_ephemeris
  };
  
  ...
  

Built-in implementation
~~~~~~~~~~~~~~~~~~~~~~~

The JSON for modern C++ library has built-in support for enumerations. For instance:

.. code-block:: cpp
  
  nlohmann::json j = 1;
  EphemerisType ephemerisType = j;                      // direct_spice_ephemeris
  
.. code-block:: cpp
  
  EphemerisType ephemerisType = interpolated_spice;
  nlohmann::json j = ephemerisType;
  std::cout << ephemerisType << std::endl;              // 3

However, the built-in implementation uses the value of the enumeration rather than its name. Clearly, it is not convenient for the user having to specify the type of the ephemeris based on how they are listed in the C++ code. Additionally, if the order of the values in the enumeration list changes in the future, converting :literal:`1` to :class:`EphemerisType` may not result in :literal:`direct_spice_ephemeris` anymore. Thus, a custom implementation is needed, in which the JSON representation of an enumeration is the name of its possible values. For instance:

.. code-block:: cpp
  
  nlohmann::json j = "direct_spice_ephemeris";
  EphemerisType ephemerisType = j;
  
Without the custom implementation, this leads to a run-time error.


Name-based implementation
~~~~~~~~~~~~~~~~~~~~~~~~~

Custom implementations for the :literal:`to_json` and :literal:`from_json` functions of all the supported enumerations are provided in the JSON Interface library. In this way, it is possible to convert :class:`nlohmann::json` objects of value type :jsontype:`string` to :class:`enum` and vice versa. However, the string representation for each enum value has to be manually provided. Although it is possible to replace the enumeration value names by identical strings during compile time, this was deemed too complex and, additionally, the enumeration names used in Tudat are not always optimal. For instance, consider this JSON file:

.. code-block:: json
  :caption: :class:`bodies.h`
  :name: bodies-json
    
  {
    "Earth": {
      "ephemeris": {
        "type": "tabulated_ephemeris"
      }
    }
  }

In this case, the :literal:`_ephemeris` part is redundant. In Tudat this is necessary when other enumerations declared in the same namespace can share the same names (e.g. :literal:`tabulated_atmosphere`). However, in a JSON file, the string :literal:`tabulated_ephemeris` will only be used inside :literal:`ephemeris` objects of body objects, so the string :literal:`tabulated` is unambiguous. Thus, when defining the string representation for existing Tudat enumerations, the redundant information is removed. Additionally, the naming convention is to use :literal:`lowerCamelCase` strings.

The definition of the string representation of the enum values is done in a file in the JSON Interface directory but in the enumeration's namespace (not in the :literal:`json_interface` namespace). This decision was made taking into account that these variables are only used inside the :literal:`to_json` and :literal:`from_json` functions of the enumeration, which must be declared in the enumeration's namespace. A map is used to define the string representation of each enumeration type:

.. code-block:: cpp
  :caption: :class:`Tudat/JsonInterface/Environment/ephemeris.h`
  :name: ephemeris-h
  
  #include <Tudat/SimulationSetup/EnvironmentSetup/createEphemeris.h>
  #include "Tudat/JsonInterface/Support/valueAccess.h"
  #include "Tudat/JsonInterface/Support/valueConversions.h"

  ...
  
  namespace tudat
  {

  namespace simulation_setup
  {

  ...
  
  //! Map of `EphemerisType`s string representations.
  static std::map< EphemerisType, std::string > ephemerisTypes =
  {
      { approximate_planet_positions, "approximatePlanetPositions" },
      { direct_spice_ephemeris, "directSpice" },
      { tabulated_ephemeris, "tabulated" },
      { interpolated_spice, "interpolatedSpice" },
      { constant_ephemeris, "constant" },
      { kepler_ephemeris, "kepler" },
      { custom_ephemeris, "custom" }
  };

  //! `EphemerisType` not supported by `json_interface`.
  static std::vector< EphemerisType > unsupportedEphemerisTypes =
  {
      custom_ephemeris
  };
  
  ...

As you can see, the string representations are provided for **all** the enumeration values, even those that are not supported by the JSON interface. For instance, :literal:`custom_ephemeris` is not supported by the JSON interface, because a :class:`std::function` cannot be provided using JSON files. Thus, this enum value is marked as unsupported by adding it to the variable :literal:`unsupportedEphemerisTypes`. In this way, when the user provides the value :literal:`"custom"`, for the key :jsonkey:`ephemeris`, rather than getting an :class:`IllevalValueError`, an :class:`EphemerisType` with value :literal:`custom_ephemeris` will be created without printing any warning. Then, when the actual :class:`EphemerisSettings` are created, in its :literal:`from_json` function, the user will get an error in which it is said that custom ephemeris is not supported by the JSON interface but it *does* exist in Tudat, so if they want to use it they have to build their own custom JSON-based C++ application, in which the ephemeris function is defined manually (after reading the JSON input file containing the remainder of the settings).

Although a :literal:`to_json` and :literal:`from_json` function has to be written for each enumeration, they each are just a single line in which the functions :literal:`stringFromEnum` and :literal:`enumFromString` defined in :class:`Tudat/JsonInterface/Support/errorHandling.h` are called:

.. code-block:: cpp
  :caption: :class:`Tudat/JsonInterface/Environment/ephemeris.h`
  :name: ephemeris-h-to-from-json
  
  #include <Tudat/SimulationSetup/EnvironmentSetup/createEphemeris.h>
  #include "Tudat/JsonInterface/Support/valueAccess.h"
  #include "Tudat/JsonInterface/Support/valueConversions.h"
  
  ...
  
  namespace tudat
  {

  namespace simulation_setup
  {

  ...

  //! Convert `EphemerisType` to `nlohmann::json`.
  inline void to_json( nlohmann::json& jsonObject, const EphemerisType& ephemerisType )
  {
      jsonObject = json_interface::stringFromEnum( ephemerisType, ephemerisTypes );
  }

  //! Convert `nlohmann::json` to `EphemerisType`.
  inline void from_json( const nlohmann::json& jsonObject, EphemerisType& ephemerisType )
  {
      ephemerisType = json_interface::enumFromString( jsonObject, ephemerisTypes );
  }
  
  ...

If the string representation of :literal:`ephemerisType` is not known, the recognised strings will be printed and an :literal:`IllegalValueError` will be thrown:

.. code-block:: txt

  Unknown string "constatn" for enum tudat::simulation_setup::EphemerisType
  Recognized strings:
    approximatePlanetPositions
    directSpice
    tabulated
    interpolatedSpice
    constant
    kepler
    custom
  libc++abi.dylib: terminating with uncaught exception of type tudat::json_interface::IllegalValueError: 
  Illegal value for key: bodies.Earth.ephemeris.type
  Could not convert value to expected type tudat::simulation_setup::EphemerisType

.. note:: All the files in which the :literal:`to_json` and :literal:`from_json` functions of enumerations and settings file are defined must include the files :class:`Tudat/JsonInterface/Support/valueAccess.h` and :class:`Tudat/JsonInterface/Support/valueConversions.h`. The former includes :class:`Tudat/JsonInterface/Support/errorHandling.h`, exposing the functions :literal:`json_interface::stringFromEnum` and :literal:`json_interface::enumFromString`.

.. caution:: When converting an enumeration to or from a :class:`nlohmann::json` object, the file in which its custom :literal:`to_json` and :literal:`from_json` functions are defined must be included (if the conversion takes place in a different file), or the functions have to be declared before being used if defined in the same file. If one forgets to include this file, the code will compile without giving any errors or warnings and the default implementation will be used, leading to a run-time error in which it is said that an :class:`int` was expected when converting to the enumeration type.
