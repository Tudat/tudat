.. _extendingJSON_valueConversions:

.. role:: jsontype
.. role:: jsonkey

Value conversions
=================

Conversion to :class:`nlohmann::json`
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The value of (the elements of) a :class:`nlohmann::json` object can be set as follows:

.. code-block:: cpp

    int myInt = 1;
    nlohmann::json j = myInt;  // j.is_numeric( ) -> true

    double myDouble = 1.0;
    nlohmann::json k;          // k.is_null( ) -> true
    k[ 0 ] = myDouble;         // k.is_array( ) && k[ 0 ].is_numeric( ) && j == k[ 0 ] -> true
    
However, this does not work (yet):

.. code-block:: cpp

    MyClass myObject = ...
    nlohmann::json j = myObject;            // compile error

because the JSON representation of :class:`MyClass` objects is not known.

At first, one could think that the :class:`nlohmann::json` class has constructors taking as an argument a basic type such as :class:`double`, :class:`int`, :class:`bool`, :class:`std::string`, :class:`std::vector` or :class:`std::map`. However, this is not true. When writing:

.. code-block:: cpp
    
    std::string myString = "hi";
    nlohmann::json j = myString;

the following is happening behind the scenes:
    
.. code-block:: cpp

    ...
    nlohmann::json j;
    std::to_json( j, "hi" );
    return j;
    
Indeed, :literal:`to_json` functions are defined for many frequently used types in the file :literal:`"json/src/json.hpp"`. These functions are always defined in the namespace of the type that is on the right-hand side of the assignment. We can create :literal:`to_json` functions for our custom classes:

.. code-block:: cpp
    :caption: :class:`dog.h`
    :name: dog-h
    
    namespace animals
    {
        class Dog
        {
        public:
            std::string name;
            double weight;
            
            Dog( const std::string& name, const double weight ) :
              name( name ), weight( weight ) { }
        };
    }
    
.. code-block:: cpp
    :caption: :class:`jsonInterface.h`
    :name: jsonInterface-h
    
    #include "json/src/json.hpp"
    #include "dog.h"
    
    namespace animals
    {            
        void to_json( nlohmann::json& j, const Dog& dog )
        {
            j[ "name" ] = dog.name;      // std::to_json( nlohmann::json&, const string& ) will be called
            j[ "weight" ] = dog.weight;  // to_json( nlohmann::json&, const double& ) will be called
        }
    }
    
.. code-block:: cpp
    :caption: :class:`jsonInterface.cpp`
    :name: jsonInterface-cpp
    
    #include "jsonInterface.h"
    
    void printDog( )
    {
        nlohmann::json jsonDog = animals::Dog( "Senda", 20 );
        std::cout << jsonDog << std::endl;                      // {"name":"Senda","weight":20}
    }
    
    
Conversion from :class:`nlohmann::json`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The value of the elements of a :class:`nlohmann::json` object can be accessed as follows:

.. code-block:: cpp

    nlohmann::json j = { "no", "yes" };
    j[ 0 ];                                                    // "no"
    j.at( 1 );                                                 // "yes"

    nlohmann::json k = { { "half", 0.5 }, { "twice", 2.0 } };
    k[ "half" ];                                               // 0.5
    k.at( "twice" );                                           // 2.0

However, as discussed in :ref:`extendingJSON_basics_valueTypes`, the returned values are not numbers or strings, but :class:`nlohmann::json` objects. Thus:

.. code-block:: cpp

    std::string str = j[ 0 ] + ", thanks";          // compile error

    std::string no = j[ 0 ];                        // "no"
    std::string str = no + ", thanks";              // "no, thanks"

The implicit conversion is done by overloading the :literal:`=` operator. This means that the following won't work:

.. code-block:: cpp

    std::string str = std::string( j[ 0 ] ) + ", thanks";     // compile error

However, an explicit conversion is always possible by calling the templated :literal:`get` method:

.. code-block:: cpp

    std::string str = j[ 0 ].get< std::string >( ) + ", thanks";     // "no, thanks"

What is happening behind the scenes, either when an implicit conversion takes place or the :literal:`get` method is called, is:

.. code-block:: cpp
    
    ...
    std::string str;                // str -> ""
    std::from_json( j, str );       // str -> "no"
    return str;

The :literal:`from_json` functions must be defined in the namespace of the type we are converting to. Many are pre-defined for frequently-used types. But, for our :class:`Dog` class, this is not possible yet:

.. code-block:: cpp

    nlohmann::json jsonDog = { { "name", "Senda" }, { "weight", 20 } };
    Dog dog = jsonDog;                                                   // compile error

We have to define its :literal:`from_json` function:

.. code-block:: cpp
    :caption: :class:`jsonInterface.h`
    :name: jsonInterface-h-from-json
    
    #include "json/src/json.hpp"
    #include "dog.h"
    
    namespace animals
    {            
        ...
        
        void from_json( const nlohmann::json& j, Dog& dog )
        {
            dog = Dog( j.at( "name" ), j.at( "weight" ) );
        }
    }
        
However, we are still getting a compile error, because this is happening behind the scenes:

.. code-block:: cpp
    
    ...
    animals::Dog dog;                  // compile error: Dog is not default-constructible!
    animals::from_json( j, dog );      
    return dog;

Thus, we need to define a default constructor for the class :class:`Dog`. We can do this by providing default values for all the arguments in the constructor:

.. code-block:: cpp
    :caption: :class:`dog.h`
    :name: dog-h-default-constructible
    
    namespace animals
    {
        class Dog
        {
        public:
            std::string name;
            double weight;
                          
            Dog( const std::string& name = "", const double weight = 0.0 ) :
                name( name ), weight( weight ) { }
        };
    }
    
Now, we can do:

.. code-block:: cpp

    nlohmann::json jsonDog = { { "name", "Senda" }, { "weight", 20 } };
    Dog dog = jsonDog;                                                   // fine



