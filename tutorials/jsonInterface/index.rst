.. _jsonInterface:

JSON Interface
==============

These pages of the wiki detail how the JSON interface can be used to run propagations by providing JSON input files and provide a few tutorials showcasing the main features.

The :literal:`json_interface` application can be used to run Tudat propagations by providing the path to a JSON input file as command-line argument in Terminal:

.. code-block:: txt

  tudatBundle/tudat/bin/json_interface main.json
  
The extension can be omitted, in which case :literal:`.json` will be assumed. The argument can also be the path to a directory containing a :literal:`main.json` file, and it can be omitted if a :literal:`main.json` file exists in the current directory.

It is also possible to use the :literal:`json_interface` application from Qt Creator. In the Projects > Run tab, choose :literal:`json_interface` in "Run configuration" and write the (absolute) path to your input file in "Command line arguments".

.. image:: qt.png

.. toctree::
   :hidden:
   :maxdepth: 2

   keys/rst/index
   modularFiles
   multicase
   customJsonApp
   tutorials/index

