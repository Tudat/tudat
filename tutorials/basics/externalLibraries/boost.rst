.. _externalBoost:

Boost: Basic Concepts
=====================
Boost is a collection of useful libraries providing a wide variety of functionality in C++ that we use. The boost library is included in the Tudat Bundle, pre-compiled and ready for use. Alternatively, you can download the library from the boost website itself, where you will find downloads and documentation for different operating systems. Note that when downloading the boost library, you will have to compile it yourself.

.. tip:: For Windows, pre-built binaries are available. In this case, compilation is not necessary by the user.

General information on boost can be found on the boost website. This website also contains documentation on all boost features. Be sure to look up the version number of the boost library that you are using when checking out this documentation. Note that the boost documentation can be quite complex, especially for people who are starting to learn and use C++ and Tudat.

Why use boost?
~~~~~~~~~~~~~~
Boost is a collection of generally well-coded C++ libraries that are a substantial addition to the standard libraries. Some boost features are particularly useful for Tudat and are frequently used for two reasons:

    1. It makes the use of some Tudat features easier.
    2. It makes the Tudat library more robust.

In Tudat, the main overall uses of the boost libraries are in cases where **dynamic memory allocation** is used and when a function is passed as an argument to another object. In case of dynamic memory allocation, use the boost shared pointer class to ensure that memory is allocated and destroyed correctly, so that no manual calls to 'delete' are required. In case of passing a function as an argument we use the boost function class, which allows for very generic interfaces and makes it easy to pass a function, where it doesn't matter whether the function is a free function or a class member function.

Although the above two cases represent the main boost features that are throughout much of the code, there are many other instances of usages of other specific boost functionalities in both the library and users' applications.

Common Boost Features
~~~~~~~~~~~~~~~~~~~~~
A couple of elements from boost are frequently used in Tudat. These are:

**Basic Features:**

    - :class:`shared_ptr`
        This is a part of the boost smart pointer library. It is a template class that stores a smart pointer to a dynamically allocated object. The smart pointer behaves much like a normal built-in C++ pointer, except that the object pointed to is guaranteed to be **automatically deleted** when the last :class:`shared_ptr` pointing to the object is destroyed or reset. This boost feature is used in Tudat because it is a safe way to use pointers to dynamically allocated objects. In Tudat it is frequently combined with :class:`make_shared`.

    - :class:`make_shared`
        This is a part of the boost smart pointer library. It is a factory function that creates an object of a given type and **returns a shared pointer** to it. This way the explicit use of new is avoided. This boost feature is frequently used in Tudat because it is an easy, efficient and safe way to use pointers to dynamically allocated objects.

**Advanced Features:**

    - :class:`function`
        The boost function library contains a family of class templates that are function object wrappers. This boost feature is used in Tudat because it is an easy way to pass a function as an input to another object or function. It is used instead of a normal C++ function pointer. The boost function type describes what the function should 'look' like, in the sense that it defines both its input types (and order) and its output type. Any function that fits this profile is accepted when it is passed. In Tudat it is frequently combined with :class:`bind` to create a link to a class member function.

    - :class:`bind`
        The boost bind library is a template library that implements a simple and versatile function argument binding mechanism. Typically, it is used to create boost function objects from either a free function or a member function. It supports arbitrary function objects, functions, function pointers, and member function pointers. It can bind any argument to a specific value or route input arguments into arbitrary positions. This boost feature is used in Tudat because it is a versatile method, and it is a good way to pass function objects.

Boost Shared Libraries (Dynamic Linking)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
During the build process of any application, multiple object files can be linked in a library. This library can be linked either statically or dynamically. Usually, we make use of the static boost libraries. These files have the extension ``.a`` and can be found in the ``boost/stage/lib directory``. However, it could be that you require the dynamically linked libraries (for instance, if you are making use of PaGMO). These files have the extension ``.dll`` and can be found in the ``boost/shared`` directory.

