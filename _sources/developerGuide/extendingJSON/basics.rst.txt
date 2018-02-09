.. _extendingJSON_basics:

.. role:: jsontype
.. role:: jsonkey

Basics
======

The external library used by the :ref:`jsonInterface` is `JSON for Modern C++ <http://github.com/nlohmann/json>`_. The basic concepts necessary to understand the following sections are provided here.

Including the library
~~~~~~~~~~~~~~~~~~~~~

.. code-block:: cpp

    #include "json/src/json.hpp"


Serialisation and deserialisation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Parse a JSON-formatted :class:`std::string` as :class:`nlohmann::json` (deserialisation):

.. code-block:: cpp
    
    std::string str = "[0, 1, 2]";
    nlohmann::json j = nlohmann::json::parse( str );


Parse a JSON file as :class:`nlohmann::json` (deserialisation):

.. code-block:: cpp
    
    std::string filePath = ...
    std::ifstream stream( filePath );
    nlohmann::json j = nlohmann::json::parse( stream );
    
Convert a :class:`nlohmann::json` object to JSON-formatted :class:`std::string` (serialisation):

.. code-block:: cpp
    
    nlohmann::json j = ...
    std::string oneLineString = j.dump( );
    std::string multiLineStringWithIndent = j.dump( 2 );

where the (optional) argument of the :literal:`dump` method is the number of spaces used for each indentation level.

Print a :class:`nlohmann::json` object as JSON-formatted text (serialisation):

.. code-block:: cpp
    
    nlohmann::json j = ...
    std::cout << j << std::endl;
    std::cout << j.dump( 2 ) << std::endl;




.. _extendingJSON_basics_valueTypes:

Value types
~~~~~~~~~~~

A :class:`nlohmann::json` object can have 6 possible value types. The value type of a :class:`nlohmann::json` object can be inferred during construction, or later during its lifecycle. The value type of a :class:`nlohmann::json` object can change.

    - :jsontype:`null`
    
      .. code-block:: cpp
          
          nlohmann::json j;                            // j.is_null( ) -> true
          
          nlohmann::json k = ...                       // k.is_null( ) -> false
          k = json( );                                 // k.is_null( ) -> true

    - :jsontype:`object`: similar to :class:`std::map`, but with the possibility to mix different value types and with keys always of type string.
    
      .. code-block:: cpp
          
          nlohmann::json j = { { "a", "zero" }, { "b", 1 } };    // j.is_object( ) -> true
          j[ "b" ];                                              // 1
          j.size( );                                             // 2
          
          nlohmann::json k;                                      // k.is_object( ) -> false
          k[ "a" ] = 0;                                          // k.is_object( ) -> true

    - :jsontype:`array`: similar to :class:`std::vector`, but with the possibility to mix different value types.
    
      .. code-block:: cpp
          
          nlohmann::json j = { 0, "one" };             // j.is_array( ) -> true
          j[ 1 ];                                      // "one"
          j.empty( );                                  // false
          
          nlohmann::json k;                                      // k.is_array( ) -> false
          k[ 0 ] = 0.5;                                // k.is_array( ) -> true
          k.push_back( "pi" );                         // 
          // k[ "b" ] = 1;                             // run-time error
          k = json( );                                 // k.is_null( ) -> true
          k[ "b" ] = 1;                                // k.is_object( ) -> true

    - :jsontype:`string`: similar to :class:`std::string`.
          
    - :jsontype:`number`: similar to :class:`double`, :class:`float` or :class:`int`.
        
    - :jsontype:`boolean`: similar to :class:`bool`.

The types :jsontype:`object` and :jsontype:`array` are known as structured types (the method :literal:`is_structured` will return :literal:`true`), while the types :jsontype:`string`, :jsontype:`number` and :jsontype:`boolean` are known as primitive types (the method :literal:`is_primitive` will return :literal:`true`).

Structured :class:`nlohmann::json` objects are containers of :class:`nlohmann::json` objects (of any value type). This means that if we create the following :class:`nlohmann::json` object of array type:

.. code-block:: cpp
    
    double aDouble = 1.0;
    std::string aString = "two";
    nlohmann::json j = { aDouble, aString };       

then when we access any of its elements (e.g. :literal:`j[ 0 ]` or :literal:`j[ 1 ]`) we will not get a :class:`double` or :class:`std::string`, but a :class:`nlohmann::json` object (with value type :jsontype:`number` or :jsontype:`string` in this case).


