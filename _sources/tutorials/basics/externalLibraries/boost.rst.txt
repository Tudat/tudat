.. _externalBoost:

Boost: Basic Concepts
=====================

Boost is a collection of useful libraries providing a wide variety of functionality in C++ that we use. The Boost library is included in the Tudat Bundle, pre-compiled and ready for use. Alternatively, you can download the library from the Boost website itself, where you will find downloads and documentation for different operating systems. Note that when downloading the Boost library, you will have to compile it yourself.

.. tip:: For Windows, pre-built binaries are available. In this case, compilation is not necessary by the user.

General information on Boost can be found on the Boost website. This website also contains documentation on all Boost features. Be sure to look up the version number of the Boost library that you are using when checking out this documentation. Note that the Boost documentation can be quite complex, especially for people who are starting to learn and use C++ and Tudat.

Why use Boost?
~~~~~~~~~~~~~~

Boost is a collection of generally well-coded C++ libraries that are a substantial addition to the standard libraries. Some Boost features are particularly useful for Tudat and are frequently used for two reasons:

   1. It makes the use of some Tudat features easier.
   2. It makes the Tudat library more robust.

Although the above two cases represent the main Boost features that are throughout much of the code, there are many other instances of usages of other specific Boost functionalities in both the library and users' applications.

Common Boost Features
~~~~~~~~~~~~~~~~~~~~~

In general, you will see the use of Boost in the following cases:

   - **Multi-dimensional Arrays**

      In C++, vectors and matrices are usually defined by using the :literal:`Eigen` library (see :ref:`externalEigen` for an introduction). Despite their usefulness, however, these linear algebra tools are only defined for up to 2-dimensional matrices. In case you wanted to store data with more dimensions, then it is common practice in Tudat to use :literal:`boost::multi_array`. You can find examples of the use of this objects, whenever files are read or written, and for the storage of aerodynamic coefficients, tabulated atmospheres, etc.

   - **Random Number Generators**

      Generation of random numbers in Tudat is done with the aid of Boost. You can see :ref:`tudatFeaturesProbabilityDistributions` for a description of the available distributions and how they are used. 

   - **Special Math Operators**

      The use of inverse hyperbolic, cumulative distribution and probability density functions are taken from the Boost libraries.

   - **Date Conversions**

      Similarly to the case above, conversion between various dates format (day/month/year to Gregorian dates) is done with Boost.

   - **File Management**

      File management with Boost allows for the creation of directories, to check if directories exist, manipulation of paths, etc.  

   - **Unit Tests**

      Testing of Tudat features is fully based on the unit testing suite provided by Boost. It allows for very structured and consistent verification of Tudat function and for simple addition of new tests.

   - **Other**

Boost Shared Libraries (Dynamic Linking)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

During the build process of any application, multiple object files can be linked in a library. This library can be linked either statically or dynamically. Usually, we make use of the static Boost libraries. These files have the extension ``.a`` and can be found in the ``boost/stage/lib directory``. However, it could be that you require the dynamically linked libraries (for instance, if you are making use of PaGMO). These files have the extension ``.dll`` and can be found in the ``boost/shared`` directory.