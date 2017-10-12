.. _matlabInterface_tutorialsInputOutput_modularInputFiles:

Modular input files
===================

This tutorial describes how to use the MATLAB Interface to generate multiple JSON files that are then combined into a single JSON object by the :literal:`json_interface` application. The code for this example can be found at:

.. code-block:: txt

  tudatBundle/matlabInterface/Examples/InputOutput/ModularInputFiles/generateInput.m


The first part of the script sets up the :class:`Simulation` object. This is done in a similar way as explained for previous examples.

Once that all the settings have been defined, instead of calling :literal:`simulatinon.run()`, we export parts of the :class:`Simulation` object to different JSON files. Before that, we change the current working directory in MATLAB to the directory in which the root JSON file (:literal:`main.json`) will be generated, so that relatives paths in the MATLAB script and in the input files coincide:

.. code-block:: matlab

  if exist('input','dir') ~= 7
      mkdir('input');
  end
  originalWorkingDirectory = cd;
  cd('input');

We also store the original working directory to set it back once that the file generation is done, and create the directory if it does not exist.

We want to export different parts of the settings to separated files. For instance, we export the body settings to a file :literal:`bodies.json` inside the :literal:`input` directory:

.. code-block:: matlab

  simulation.bodies = json.modular(simulation.bodies,'bodies.json');

The :literal:`json.modular` function exports the object specified as first argument to a file specified by the second argument, and returns the string :literal:`$(bodies.json)` in this case, which is assigned to :literal:`simulation.bodies`. In this way, the simulation object does not store the bodies anymore, but a reference to the file containing that information, by means of the special string :literal:`$()`.

We can do the same for other settings:

.. code-block:: matlab

  simulation.propagators = {json.modular(propagator,'translationalPropagator.json')};
  simulation.integrator = json.modular(integrator,'rk4.json');
  simulation.export = json.modular(simulation.export,'export.json');

Finally, we export the main JSON file:

.. code-block:: matlab

  json.export(simulation,'main.json');

which looks like this:

.. code-block:: json
  :caption: :class:`matlabInterface/Examples/InputOutput/ModularInputFiles/input/main.json`
  
  {
    "initialEpoch": 4.77171E+8,
    "finalEpoch": 4.7760480000001341E+8,
    "globalFrameOrigin": "SSB",
    "globalFrameOrientation": "J2000",
    "spice": {
      "useStandardKernels": true
    },
    "bodies": "$(bodies.json)",
    "propagators": [
      "$(translationalPropagator.json)"
    ],
    "integrator": "$(rk4.json)",
    "export": "$(export.json)",
    "options": {
      "printInterval": 86400,
      "fullSettingsFile": "@path(fullSettings.json)"
    }
  }

After running :literal:`json_interface main.json`, the directory tree looks like this:

.. code-block:: txt

  ModularInputFiles
  | 
  | input
  |     |
  |     | bodies.json
  |     | export.json
  |     | fullSettings.json
  |     | main.json
  |     | rk4.json
  |     | translationalPropagator.json
  |      
  | output 
  |      |
  |      | epochs.txt
  |      | states.txt
  |      
  | generateInput.m
  | processOutput.m

The script :literal:`processOutput.m` can be used to generate a plot using the output files.
