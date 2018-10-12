.. _extendingJSON_enhancedValueSet:

.. role:: jsontype
.. role:: jsonkey

Enhanced value set
==================

In the previous section, the enhanced value access functionality has been introduced. This functionality allows the user to easily validate their input JSON files by providing comprehensible error messages in case non-defaultable properties are missing or the provided values are not convertible to the expected type. This was achieved by using :class:`KeyPath` s, special keys and the :literal:`getValue` function.

However, there is nothing like a :literal:`setValue` function, and :class:`KeyPath` s are not used in :literal:`to_json` functions. This is so because there is no need for validation when converting Tudat objects to JSON. If the process is not successful, it will be the developer's fault and not the user's. Additionally, although it is not always possible to access an arbitrary key, since it may be undefined, it is always possible to define any arbitrary key for a :class:`nlohmann::json` object.

Nonetheless, there are a few functions provided by the :literal:`json_interface` that can make the writing of :literal:`to_json` functions easier, namely:

- :literal:`void assignIfNotNaN( nlohmann::json& j, const std::string& key, const EquatableType& value )`: updates or defines :literal:`j[ key ]` only if :literal:`value == value` (this comparison returns :literal:`false` for :literal:`TUDAT_NAN`). Can only be used when the comparison operator is defined for :literal:`EquatableType`.

- :literal:`void assignIfNotNull( nlohmann::json& j, const std::string& key, const std::shared_ptr< T >& object )`: updates or defines :literal:`j[ key ]` only if :literal:`object != NULL`.

- :literal:`void assignIfNotEmpty( nlohmann::json& j, const std::string& key, const ContainerType& container )`: updates or defines :literal:`j[ key ]` only if :literal:`object.empty( ) == false`. Can only be used when the method :literal:`empty` is defined for :literal:`ContainerType`.

Note that the following is not possible:

.. code-block:: cpp
  
  nlohmann::json mainJson;
  KeyPath keyPath = "integrator" / "type";
  mainJson[ keyPath ] = "rungeKutta4";             // compile error
  mainJson[ "integrator.type" ] = "rungeKutta4";   // wrong! {"integrator.type":"rungeKutta4"}

Instead, the default basic value set methods have to be used:

.. code-block:: cpp
  
  nlohmann::json mainJson;
  mainJson[ "integrator" ][ "type" ] = "rungeKutta4";

If one wants to use the properties of a given object to update the keys of a :class:`nlohmann::json` object that is above in the key tree, e.g.:

.. code-block:: cpp
  
  inline void to_json( nlohmann::json& j, const Integrator& integrator )
  {
      j[ "type" ] = integrator.type;
      j[ "stepSize" ] = integrator.stepSize;
      j[ "<-" / "initialEpoch" ] = integrator.initialTime;   // compile error
  }

this is not possible (yet). This kind of implementation is required only once in the whole JSON Interface library, when converting a :class:`Simulation` object to :class:`nlohmann::json`, as the information contained in the propagator settings is stored in three different places in the simulation object (part of its information is stored in the key :jsonkey:`propagators` of the :literal:`mainJson`, the termination settings are stored in the key :jsonkey:`termination` of the :literal:`mainJson`, and the print interval is stored in the :jsonkey:`options.printInterval` key path of the :literal:`mainJson`). To do so, the following implementation was chosen (code has been simplified):

.. code-block:: cpp
  :caption: :class:`simulation.h`
  :name: simulation-h
  
  void to_json( nlohmann::json& mainJson, const Simulation& simulation )
  {
      ...
      mainJson[ "integrator" ] = simulation.integrator;
      ...
      propagators::to_json( mainJson, simulation.propagator );
  }

.. code-block:: cpp
  :caption: :class:`propagator.h`
  :name: propagator-h
  
  namespace propagators
  {
      void to_json( nlohmann::json& mainJson, const Propagator& propagator )
      {
          mainJson[ "propagators" ][ 0 ][ "type" ] = propagator.type;
          mainJson[ "propagators" ][ 0 ][ "centralBodies" ] = propagator.centralBodies;
          ...
          mainJson[ "termination" ] = propagator.terminationSettings;
          mainJson[ "options" ][ "printInterval" ] = propagator.printInterval;
      }
  }

This is the only case in which a :literal:`to_json` function is manually called in the JSON Interface library. Note that, when passed to the :literal:`to_json` function, the :literal:`mainJson` object is not re-initialised, so the keys defined before this function call are kept.
