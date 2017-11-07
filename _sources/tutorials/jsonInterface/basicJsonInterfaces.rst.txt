.. _basicJsonInterfaces:

Basic JSON interfaces
=====================

There are two ways of using the JSON interface with Tudat. Either full numerical simulation are done with a JSON file which provides propagation output. The second option is using a JSON file for loading simulation settings and extend it with C++ code yourself. This can for example be used for writing guidance. The second interface is explained in :ref:`jsonInterface_tutorials_lifetimeMaximisationCPP`.  The first option can be called in Terminal and Qt:


The :literal:`json_interface` application can be used to run Tudat propagations by providing the path to a JSON input file as command-line argument in Terminal::

  tudatBundle/tudat/bin/json_interface main.json
  
.. note:: The extension can be omitted, in which case :literal:`.json` will be assumed. The argument can also be the path to a directory containing a :literal:`main.json` file, and it can be omitted if a :literal:`main.json` file exists in the current directory.

It is also possible to use the :literal:`json_interface` application from Qt Creator. In the Projects > Run tab, choose :literal:`json_interface` in "Run configuration" and write the (absolute) path to your input file in "Command line arguments".

.. image:: qt.png
 
Note that details about the contents of the JSON file is explained in upcoming :ref:`jsonInterface_tutorials`.
