.. _matlabInterface_tutorialsInputOutput_basicPropagation:

Unperturbed Earth-orbiting Satellite
====================================

This tutorial describes how to propagate the orbit of an unperturbed satellite about Earth using the MATLAB Interface to generate an input file that is then provided to the :literal:`json_interface` application. This example is similar (but not identical) to the example :ref:`walkthroughsUnperturbedEarthOrbitingSatellite` written in C++. The code for this example can be found at:

.. code-block:: txt

  tudatBundle/matlabInterface/Examples/InputOutput/BasicPropagation/generateInput.m


The first part of the script sets up the :class:`Simulation` object. This is identical to the initial part of :ref:`matlabInterface_tutorialsSeamless_basicPropagation`, so it is not repeated here.

Once that all the settings have been defined, instead of calling :literal:`simulatinon.run()`, we export the :class:`Simulation` object to a JSON file. However, we have to specify first the results to be exported to output files:

.. code-block:: matlab

  inputFile = fullfile('input','main.json');
  simulation.addResultsToExport(fullfile('..','output','results.txt'),{'independent','state'});

Here, the :literal:`fullfile` function is used to build a platform-dependent path (e.g. :literal:`input/main.json` on Unix systems, :literal:`input\\main.json` on Windows). In the first line, we declare the path to which the input file will be saved (will be used later). Then, in the second line, we request the variables :literal:`independent` and :literal:`state` to be saved to a file at :literal:`../output/results.txt`.

Then, we request the application to generate a file containing the full settings actually used for the propagation (including all the default values):

.. code-block:: matlab

  simulation.options.fullSettingsFile = 'fullSettings.json';

This file will be generated inside the :literal:`input` directory, because paths are specified relative to the root input file (:literal:`input/main.json` in this case), and not relative to the directory in which the MATLAB script is located.

The last step is to generate the JSON file:

.. code-block:: matlab

  json.export(simulation,inputFile);

After running :literal:`json_interface main.json`, the resulting directory tree looks like this:

.. code-block:: txt

  BasicPropagation
  | 
  | input
  |     |
  |     | fullSettings.json
  |     | main.json
  |      
  | output 
  |      |
  |      | results.txt
  |      
  | generateInput.m
  | processOutput.m

The script :literal:`processOutput.m` can be used to generate a plot using the results file.
