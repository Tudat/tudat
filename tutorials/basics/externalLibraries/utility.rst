.. _externalUtility:

Dynamic Memory and Function Objects
***********************************

Dynamic memory and function objects are two concepts that are widely used in Tudat. They allow for a much more streamlined and powerful implementation of the C++ language, and greatly simplify application and features implementation. 

These two features used to be part of the :literal:`boost` library (introduced in :ref:`externalBoost`), but their use became so extensive and their utility so considerable, that they were included in the standard C++ libraries, after C++11. 

On these page, the basic concepts of these C++ features will be explained, at a basic level. Together with the tutorial provided in :ref:`externalUtilityExamples`, you should get an idea of how these are used throughout Tudat and of how you can implement them in your own applications.

Dynamic Memory
==============

In Tudat, dynamic memory allocation is used to ensure that memory is allocated and destroyed correctly (and automatically), such that no manual calls to :literal:`delete` are required.

Pointers are a very common example of where these features come in handy. The allocation and descrution of pointers in C++ can be a hard task to understand, when first learning to use this language. Dynamic memory management helps reducing this gap by introducing the concept of smart pointers.

Dynamic memory management is defined in the header :literal:`<memory>`. You can find a complete description of the other elements available in this header, `here <https://en.cppreference.com/w/cpp/memory>`_.

Common Features
~~~~~~~~~~~~~~~

While going through virtually any Tudat application or function, you will routinely encounter a couple of elements related to dynamic memory. These are:

	.. class:: shared_ptr

		It is a template class that stores a smart pointer to a dynamically allocated object. The smart pointer behaves much like a normal C++ pointer, except that the object pointed to is guaranteed to be **automatically deleted** when the last :class:`shared_ptr` pointing to the object is destroyed or reset. This feature is used in Tudat because it is a safe way to use pointers to dynamically allocated objects. In Tudat it is frequently combined with :class:`make_shared`.

	.. class:: make_shared

		It is a factory function that creates an object of a given type and **returns a shared pointer** to it. This way the explicit use of new is avoided. This feature is frequently used in Tudat because it is an easy, efficient and safe way to use pointers to dynamically allocated objects.

Function Objects
================

In case of passing a function as an argument we use the :literal:`function` class, which allows for very generic interfaces and makes it easy to pass a function, where it doesn't matter whether the function is a free function or a class member function.

Function objects are defined in the header :literal:`<functional>`. You can find a complete description of the other elements available in this header, `here <https://en.cppreference.com/w/cpp/utility/functional>`_.

Common Features
~~~~~~~~~~~~~~~

   .. class:: function

      The :class:`function` function contains a family of class templates that are function object wrappers. This feature is used in Tudat because it is an easy way to pass a function as an input to another object or function. It is used instead of a normal C++ function pointer. The :class:`function` type describes what the function should 'look' like, in the sense that it defines both its input types (and order) and its output type. Any function that fits this profile is accepted when it is passed. In Tudat it is frequently combined with :class:`bind` to create a link to a class member function.

   .. class:: bind

      The :class:`bind` function is a template library that implements a simple and versatile function argument binding mechanism. Typically, it is used to create :class:`function` objects from either a free function or a member function. It supports arbitrary function objects, functions, function pointers, and member function pointers. It can bind any argument to a specific value or route input arguments into arbitrary positions. This feature is used in Tudat because it is a versatile method, and it is a good way to pass function objects.