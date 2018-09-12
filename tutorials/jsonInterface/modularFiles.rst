.. _jsonInterface_modularFiles:

.. role:: jsontype
.. role:: jsonkey

Modular Files
=============

The JSON Interface allows the use of modular JSON input files. This means that a simulation can be set up by using different files. Since the list of settings to be provided can be very long, in some cases it is convenient to specify those settings in different files. Additionally, modularity allows re-using parts from other projects, and a set of files containing settings for e.g. several integrators can be provided and then only one of them be used by referencing it from the root input file (i.e. the file provided as command-line argument to the :literal:`json_interface` application).

Modularity can be nested, i.e. not only the root input file, but any file, can contain references to other files. The contents of a file can be included by using the following special string:

.. code-block:: json
  :caption: :class:`main.json`

  {
    "integrator": "$(rk78.json)"
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

After parsing by the JSON Interface, this will result in:

.. code-block:: json
  :caption: :class:`main.json`

  {
    "integrator": {
      "type": "rungeKuttaVariableStepSize",
      "rungeKuttaCoefficientSet": "rungeKuttaFehlberg78",
      "initialStepSize": 20,
      "minimumStepSize": 1,
      "maximumStepSize": 1e4
    }
  }

The path can be relative to the directory of the file in which it is specified, or absolute. The extension can be omitted, in which case :literal:`.json` will be assumed. The path to a directory can also be provided if that directory contains a :literal:`main.json` file.

It is also possible to refer to parts of a file, rather than importing all of its contents. This is done by using the following syntax:

.. code-block:: json

  "$(file.json){variables}"

where :literal:`variables` is a comma-separated list of key paths to variables defined in :literal:`file.json`. For instance, imagine that we define a set of variables in a separate file:

.. code-block:: txt

  root
  | 
  | main.json
  | variables.json



.. code-block:: json
  :caption: :class:`variables.json`

  {
    "satellite": {
      "altitude": {
        "body": "satellite",
        "dependentVariableType": "altitude",
        "relativeToBody": "Earth"
      },
      "drag": {
        "body": "satellite",
        "dependentVariableType": "acceleration",
        "accelerationType": "aerodynamic",
        "bodyExertingAcceleration": "Earth"
      },
      "srp": {
        "body": "satellite",
        "dependentVariableType": "acceleration",
        "accelerationType": "cannonBallRadiationPressure",
        "bodyExertingAcceleration": "Sun"
      }
    }
  }

Since these variables are used in several places (for defining termination conditions and export settings), rather than repeating the objects we can reference parts of this file. For instance, in the root file we can write:

.. code-block:: json
  :caption: :class:`main.json`

  {
    "termination": {
      "variable": ["$(variables.json){satellite.altitude}"],
      "lowerLimit": 110E+3
    }
    "export": [
      {
        "file": "finalAltitude.txt",
        "variables": ["$(variables.json){satellite.altitude}"],
        "onlyFinalStep": true,
        "epochsInFirstColumn": false
      },
      {
        "file": "accelerations.txt",
        "variables": "$(variables.json){satellite.drag,satellite.srp}"
      }
    ]
  }

The special string :literal:`"$(variables.json){satellite.altitude}"` is replaced by the object:

.. code-block:: json

  {
    "body": "satellite",
    "dependentVariableType": "altitude",
    "relativeToBody": "Earth"
  }
  
while the string :literal:`"$(variables.json){satellite.drag,satellite.srp}"` is replaced by an array of objects:

.. code-block:: json

  [
    {
      "body": "satellite",
      "dependentVariableType": "acceleration",
      "accelerationType": "aerodynamic",
      "bodyExertingAcceleration": "Earth"
    },
    {
      "body": "satellite",
      "dependentVariableType": "acceleration",
      "accelerationType": "cannonBallRadiationPressure",
      "bodyExertingAcceleration": "Sun"
    } 
  ]

The key :literal:`export.variables` expects an array of objects. Then, we need to explicitly specify that we want to create an array with one element by writing either :literal:`[ "$(variables.json){satellite.altitude}" ]` or adding a comma at the end of the special string's list of variables, as in :literal:`"$(variables.json){satellite.altitude,}"`.

Now imagine we provide a file containing pre-defined settings for celestial bodies:

.. code-block:: txt

  root
  | 
  | bodies.json
  | main.json

.. code-block:: json
  :caption: :class:`bodies.json`

  {
    "Sun": { ... },
    "Mercury": { ... },
    "Venus": { ... },
    "Earth": { ... },
    "Moon": { ... },
    "Mars": { ... },
    "Jupiter": { ... },
    "Saturn": { ... },
    "Uranus": { ... },
    "Neptune": { ... }
  }

Then, in the main file, if we want to propagate the dynamics of the inner solar system, we can write:

.. code-block:: json
  :caption: :class:`main.json`

  {
    "bodies": {
      "Sun": "$(bodies.json){Sun}",
      "Mercury": "$(bodies.json){Mercury}",
      "Venus": "$(bodies.json){Venus}",
      "Earth": "$(bodies.json){Earth}",
      "Moon": "$(bodies.json){Moon}",
      "Mars": "$(bodies.json){Mars}"
    }
  }

However, we can also obtain the same :jsonkey:`bodies` object using a compact expression:

.. code-block:: json
  :caption: :class:`main.json`

  {
    "bodies": "$(bodies.json){Sun:Sun,Mercury:Mercury,Venus:Venus,Earth:Earth,Moon:Moon,Mars:Mars}"
  }
