.. _jsonInterface:

JSON Interface
==============

The JSON inteface makes it possible to use Tudat without knowledge of C++. In fact, the JSON interface creates the C++ settings classes, as described in :ref:`tudatFeaturesIndex`, by itself. These pages of the wiki detail how the JSON interface can be used to run propagations by providing JSON input files and provide a few tutorials showcasing the main features. 

The usage of JSON can be done in two ways. The first method, which is exaplained in :ref:`basicJsonInterfaces`, consists of only the use of :literal:`.json` files, and the :literal:`json_interface` library. For the second method, JSON is combined with C++, to extend its functionalities. The latter method is explained in :ref:`jsonInterface_customJsonApp`.

.. toctree::
   :numbered:
   :maxdepth: 2

   basicJsonInterfaces
   customJsonApp
   tutorials/index
   keys/rst/index
   modularFiles
   multicase
