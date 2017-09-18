.. _jsonInterface:

JSON Interface
==============

These pages of the wiki detail how the JSON interface can be used to run propagations by providing JSON input files and provide a few tutorials showcasing the main features.

The application :literal:`json_interface` can be used to run Tudat propagations by providing the path to a JSON input file as command-line argument:

.. code-block:: txt

  tudatBundle/tudat/bin/json_interface main.json
  
The extension can be omitted, in which case :literal:`.json` will be assumed. The argument can also be the path to a directory containing a :literal:`main.json` file, and it can be omitted if a :literal:`main.json` file exists in the current directory.

.. toctree::
   :maxdepth: 2

   keys/index
   modularFiles
   multicase
   tutorials/index

