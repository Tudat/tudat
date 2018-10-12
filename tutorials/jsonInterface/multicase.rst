.. _jsonInterface_multicase:

.. role:: jsontype
.. role:: jsonkey

Multi-case Simulations
======================

The JSON Interface allows combining several input files to generate a single object that is used to set up a propagation. This is slightly different than the modular files described in the previous page. In the previous case, during the creation of the root input file, the special strings :literal:`"$(file.json)"` were replaced by the contents of the reference file. In this case, the contents of two files are combined, processing first one object and then (re-)defining the keys specified in the second object. This is done by writing, for instance:

.. code-block:: json
  :caption: :class:`mass8000.json`

  [
    "$(shared.json)",
    {
      "bodies.satellite.mass": 8000
    }
  ]

The first element of the array is (a reference to a file containing) the settings for the simulation (the file referred to as root or main file so far). Then, the second object is used to (re-)define the values of some keys by providing key paths. For instance, if the contents of the reference file are:

.. code-block:: json
  :caption: :class:`shared.json`

  {
    "bodies": {
      "Earth": {
        "useDefaultSettings": true
      },
      "satellite": {
        "mass": 5000
      }
    }
  }

parsing the file :class:`mass8000.json` will result in the following object:

.. code-block:: json

  {
    "bodies": {
      "Earth": {
        "useDefaultSettings": true
      },
      "satellite": {
        "mass": 8000
      }
    }
  }

Even when the keys :jsonkey:`satellite.mass` or :jsonkey:`satellite` are undefined in :class:`shared.json`, the same combined object would be obtained.

This feature is especially useful when running multi-case simulations (for solving optimisation problems), as it allows the user to define the shared settings in a file, that is referenced from mergeable files in which some keys (the optimisation variables) are (re-)defined. In this way, repeating the same information in every file is avoided, which can lead to large file sizes specially for complex propagations or when the number of simulations to be run is large.

For instance, the following tree structure can be used to run the same simulations with a variable value for the mass of the vehicle: ::

  root
  | 
  | inputs
  |      |
  |      | mass5000.json
  |      | mass6000.json
  |      | mass7000.json
  |      | mass8000.json
  |      
  | shared.json

Then, the :literal:`json_interface` application is called for each of the files inside the :class:`inputs` directory. Note that there is no need to call it for the file :class:`shared.json`. If one has `GNU Parallel <https://www.gnu.org/software/parallel/>`_ installed, it is possible to run the simulations in parallel by writing in Terminal: ::

  parallel json_interface ::: inputs/*.json


Note that the second element of the array to be merged can contain several keys to be (re-)defined, and even more than one object can be provided. For instance, consider the following tree structure: ::

  root
  | 
  | inputs
  |      | 
  |      | rk4
  |      |   | 
  |      |   | mass5000.json
  |      |   | mass6000.json
  |      |   | mass7000.json
  |      |   | mass8000.json
  |      | 
  |      | rk78
  |           | 
  |           | mass5000.json
  |           | mass6000.json
  |           | mass7000.json
  |           | mass8000.json
  |      
  | shared.json
  | rk4.json
  | rk78.json


.. code-block:: json
  :caption: :class:`rk4.json`
  
  {
    "type": "rungeKutta4",
    "stepSize": 20
  }


.. code-block:: json
  :caption: :class:`rk78.json`
  
  {
    "type": "rungeKuttaVariableStepSize",
    "rungeKuttaCoefficientSet": "rungeKuttaFehlberg78",
    "initialStepSize": 20,
    "minimumStepSize": 1,
    "maximumStepSize": 1e4
  }

Then, each of the mergeable files would look like this:

.. code-block:: json
  :caption: :class:`inputs/rk4/mass8000.json`

  [
    "$(../../shared.json)",
    {
      "integrator": "$(../../rk4.json)",
      "bodies.satellite.mass": 8000
    }
  ]


It is also possible to use the following integrator file:

.. code-block:: json
  :caption: :class:`rk4.json`
  
  [
    "$(shared.json)",
    {
      "integrator": {
        "type": "rungeKutta4",
        "stepSize": 20
      }
    }
  ]

and then

.. code-block:: json
  :caption: :class:`inputs/rk4/mass8000.json`

  [
    "$(../../rk4.json)",
    {
      "bodies.satellite.mass": 8000
    }
  ]

In this case, the :class:`mass8000.json` file is a mergeable file that references another mergeable file. The mergeable :class:`rk4.json` loads first the data defined in :class:`shared.json` and then defines the key :jsonkey:`integrator` to be equal to the provided object (i.e. :literal:`{ "type": "rungeKutta4", "stepSize": 20 }`). In both cases, the final merged object used to actually set up the simulation will be identical.

.. caution:: Note that the following file would result in a different behaviour:

  .. code-block:: json
    :caption: :class:`inputs/rk4/mass8000.json`
  
    [
      "$(../../rk4.json)",
      {
        "bodies": {
          "satellite": {
            "mass": 8000
          }
        }
      }
    ]

  In this case, the contents of :class:`rk4.json` would be loaded first, and then the property of the key :jsonkey:`bodies` would be re-defined to be equal to :literal:`{ "satellite": { "mass": 8000 } }`. This would result in the loss of other keys defined inside :literal:`bodies.satellite` and :literal:`bodies` in the file :class:`shared.json`.

Generally, we will want to save the results of each simulation (e.g. the epochs and states) to a different file, so that we end up with the following file tree: ::

  root
  | 
  | inputs
  |      | 
  |      | mass5000.json
  |      | mass6000.json
  |      | mass7000.json
  |      | mass8000.json
  |      
  | outputs
  |       | 
  |       | mass5000.txt
  |       | mass6000.txt
  |       | mass7000.txt
  |       | mass8000.txt
  |      
  | shared.json

We can do this by defining the key :jsonkey:`export` in the :class:`shared.json` file:

.. code-block:: json
  :caption: :class:`shared.json`
  
  {
    "export": {
      "variables": [
        {
          "type": "independent"
        },
        {
          "type": "state"
        }
      ]
    }
  }

and then, in each file inside the :class:`inputs` directory, we would have to define the file to which the results of that simulation should be saved:

.. code-block:: json
  :caption: :class:`inputs/mass8000.json`

  [
    "$(../shared.json)",
    {
      "bodies.satellite.mass": 8000,
      "export.file": "@path(../outputs/mass8000.txt)"
    }
  ]

However, there is a way to avoid having to include this additional line in each of the input files. Before showing how this can be done, it is necessary to define the following concepts:

  - **Declaration file**: file in which a JSON key and corresponding value are defined.
  - **Parent file**: file from which the declaration file is referenced.
  - **Root file**: file provided as input argument to the :literal:`json_interface` application.

Then, the following special variables can be used inside strings anywhere in input files, which will be replaced by:

  - :literal:`${FILE_DIR}`: absolute path of the directory where the declaration file is located.
  - :literal:`${FILE_STEM}`: filename (without extension) of the declaration file.
  - :literal:`${FILE_NAME}`: filename (with extension) of the declaration file.
  - :literal:`${PARENT_FILE_DIR}`: absolute path of the directory where the parent file is located.
  - :literal:`${PARENT_FILE_STEM}`: filename (without extension) of the parent file.
  - :literal:`${PARENT_FILE_NAME}`: filename (with extension) of the parent file.
  - :literal:`${ROOT_FILE_DIR}`: absolute path of the directory where the root file is located.
  - :literal:`${ROOT_FILE_STEM}`: filename (without extension) of the root file.
  - :literal:`${ROOT_FILE_NAME}`: filename (with extension) of the root file.

For instance, we can remove the line :literal:`export.file = ...` from each of the individual input files by writing in the shared file:

.. code-block:: json
  :caption: :class:`shared.json`
  
  {
    "export": {
      "variables": [
        {
          "type": "independent"
        },
        {
          "type": "state"
        }
      ],
      "file": "@path(outputs/${ROOT_FILE_STEM}.txt)"
    }
  }

When running :literal:`json_interface mass8000.json`, the string :literal:`"@path(outputs/${ROOT_FILE_STEM}.txt)"` will be resolved to :literal:`"../outputs/mass8000.txt"`. Note that the path, defined relative to the input file in which it was declared, is converted to a path that is relative to the root input file provided to the application as command-line argument (i.e. :class:`mass8000.json`). This is only possible by using :literal:`"@path(relPath)"`. If we do not use the :literal:`@path` keyword, the JSON parser cannot tell regular strings and paths apart, and in this case that would have resulted in the creation of an :class:`outputs` directory inside the :class:`inputs` directory. In summary:

  - When a path is provided as a plain string (e.g. :literal:`"relativePath"`), it must be relative to the root file.
  - When a path is provided using the :literal:`@path` keyword (e.g. :literal:`"@path(relativePath)"`), it must be relative to the declaration file.
  - When a special string is used to include (parts of) the contents of another JSON file (e.g. :literal:`"$(shared.json)"`), the :literal:`@path` keyword is not used and it must be relative to the declaration file.

It is recommended to never provide relative paths as plain strings, and to always use either :literal:`@path( )` or :literal:`$( )`, so that the paths are always specified relative to the declaration file. When no modular or mergeable files are used, the root file and the definition file are always the same, so using the :literal:`@path` keyword makes it unnecessary but still recommended, as the project could be modularised in the future or parts of it may be end up being used in other projects, potentially requiring the use of the :literal:`@path` keyword.
