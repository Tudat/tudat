.. _externalBoost:

Boost: Basic Concepts
=====================

Boost is a collection of useful libraries providing a wide variety of functionality in C++ that we use. The boost library is included in the Tudat Bundle, pre-compiled and ready for use. Alternatively, you can download the library from the boost website itself, where you will find downloads and documentation for different operating systems. Note that when downloading the boost library, you will have to compile it yourself.

.. tip:: For Windows, pre-built binaries are available. In this case, compilation is not necessary by the user.

General information on boost can be found on the boost website. This website also contains documentation on all boost features. Be sure to look up the version number of the boost library that you are using when checking out this documentation. Note that the boost documentation can be quite complex, especially for people who are starting to learn and use C++ and Tudat.

Why use Boost?
~~~~~~~~~~~~~~

Boost is a collection of generally well-coded C++ libraries that are a substantial addition to the standard libraries. Some boost features are particularly useful for Tudat and are frequently used for two reasons:

   1. It makes the use of some Tudat features easier.
   2. It makes the Tudat library more robust.

Although the above two cases represent the main boost features that are throughout much of the code, there are many other instances of usages of other specific boost functionalities in both the library and users' applications.

Boost Shared Libraries (Dynamic Linking)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

During the build process of any application, multiple object files can be linked in a library. This library can be linked either statically or dynamically. Usually, we make use of the static boost libraries. These files have the extension ``.a`` and can be found in the ``boost/stage/lib directory``. However, it could be that you require the dynamically linked libraries (for instance, if you are making use of PaGMO). These files have the extension ``.dll`` and can be found in the ``boost/shared`` directory.