.. _preknowledgeIndex:

Tudat Basics
============

These pages of the wiki will help you build a strong knowledge basis to get started with Tudat. It treats some functions of external libraries often used in Tudat and concepts which are mandatory to understand the Tudat example codes discussed later. But also, it contains some additional concepts which are used within Tudat, but are not necessary to understand the Tudat example code. The following concepts should be understood before proceeding further:

- Eigen common types, as discussed in discussed :ref:`externalEigen`
- Boost features, as discussed in :ref:`externalBoost`
- Dynamic memory (such as :literal:`shared_ptr` and :literal:`make_shared`) and function objects (such as :literal:`function` and :literal:`bind`), as discussed in :ref:`externalUtility`

Additional information/concepts include basic C++ tutorials, the use iterators and more advanced Eigen and boost concepts. Don't forget to check these sections when encountering unclear used concepts in the Tudat code.


.. toctree::
   :numbered:
   :maxdepth: 1

   externalLibraries/eigen
   externalLibraries/eigenExamples
   externalLibraries/utility
   externalLibraries/utilityExamples
   externalLibraries/cppTutorials
   externalLibraries/boost
