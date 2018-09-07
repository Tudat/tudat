.. _basicJsonInterfaces:

Basic JSON Interfaces
=====================

There are two ways of using the JSON interface with Tudat. In the first case, you can write a full numerical simulation in a JSON file which provides propagation output, whereas in the second option is to use a JSON file to load simulation settings and extend it with C++ code yourself. The latter is usefult if, for example, you need to write a custom guidance algorithm for your application. The second interface is explained in :ref:`jsonInterface_tutorials_lifetimeMaximisationCPP`.  The first option can be called in Terminal and Qt by following the procedure below.

The :literal:`json_interface` application can be used to run Tudat propagations by providing the path to a JSON input file as command-line argument in Terminal::

	tudatBundle/tudat/bin/json_interface main.json

.. note:: The extension can be omitted, in which case :literal:`.json` will be assumed. The argument can also be the path to a directory containing a :literal:`main.json` file, and it can be omitted if a :literal:`main.json` file exists in the current directory.

It is also possible to use the :literal:`json_interface` application from Qt Creator. In the Projects > Run tab, choose :literal:`json_interface` in "Run configuration" and write the (absolute) path to your input file in "Command line arguments", as shown in the image below.

.. image:: qt.png

Note that details about the contents of the JSON file is explained in the upcoming :ref:`jsonInterface_tutorials` page.