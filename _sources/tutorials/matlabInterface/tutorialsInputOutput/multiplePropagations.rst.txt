.. _matlabInterface_tutorialsInputOutput_multiplePropagations:

Multiple propagations
=====================

This tutorial describes how to use the MATLAB Interface to generate multiple JSON files that are then provided to the :literal:`json_interface` application to solve an optimisation process (or run several cases sharing a large part of the settings). The code for this example can be found at:

.. code-block:: txt

  tudatBundle/matlabInterface/Examples/InputOutput/MultiplePropagations/generateInput.m


The first part of the script sets up the :class:`Simulation` object. This is identical to the initial part of :ref:`matlabInterface_tutorialsSeamless_multiplePropagations`, so it is not repeated here.

Once that all the settings have been defined, instead of calling :literal:`simulatinon.run()`, we export the :class:`Simulation` object to a JSON file. However, we have to specify first the results to be exported to output files:

.. code-block:: matlab

  simulation.addResultsToExport(fullfile('..','output','${ROOT_FILE_STEM}.txt'),{'independent','state',altitude});

Here, the :literal:`fullfile` function is used to build a platform-dependent path (using :literal:`/` on Unix systems and :literal:`\\` on Windows). We defined request the variables :literal:`independent`, :literal:`state` and :literal:`altitude` to be saved to a file located at :literal:`../output/results.txt`. The JSON files will be located inside a :literal:`input` directory, so the output files will be generated inside :literal:`matlabInterface/Examples/InputOutput/MultiplePropagations/output` since the provided file is relative to the file in which it will be declared. We use the special string :literal:`${ROOT_FILE_STEM}`, which, during run-time, will be replaced by the filename (without extension) of the JSON file that is being used as input argument for the :literal:`json_interface` application (e.g. when running :literal:`json_interface mass200.json`, the :literal:`${ROOT_FILE_STEM}` of is :literal:`mass200`).

Then, since we are not running the propagations from MATLAB, we cannot use :literal:`fprintf` to print a message when each propagation starts/terminate. Thus, we set the option:

.. code-block:: matlab

  simulation.options.notifyOnPropagationTermination = true;

to print a message on the Terminal window when each propagation ends.

Now we have defined all the settings that are shared amongst all the propagations. Thus, we can generate a JSON file containing these shared settings:

.. code-block:: matlab

  json.export(simulation,fullfile('input','shared.json'));

Now, we generate the input files in which the values of the mass and reference area of the body :literal:`satellite` are defined, after including the contents of :literal:`shared.json`.

.. code-block:: matlab

  masses = 200:100:800;
  referenceAreas = 8:-1:2;
  for i = 1:length(masses)
      m = masses(i);
      A = referenceAreas(i);
      filePath = fullfile('input',sprintf('mass%i.json',m));
      fileContents = json.merge('$(shared.json)','bodies.satellite.mass',m,'bodies.satellite.referenceArea',A);
      json.export(fileContents,filePath);
  end

The :liteal:`json.merge` function returns an array in which the first element is the object passed as first argument to the function, and in which the second element is a map/object containing a set of key-value pairs. For instance, one the file generated in the first iteration of the loop is:

.. code-block:: json
  :caption: :class:`matlabInterface/Examples/InputOutput/MultiplePropagations/output/mass200.json`

  [
    "$(shared.json)",
    {
      "bodies.satellite.mass": 200,
      "bodies.satellite.referenceArea": 8
    }
  ]
  
Now that all the input files have been generated, we can run each of them, e.g. :literal:`json_interface mass200.json` (we must not run the file :literal:`shared.json`). We can also use `GNU Parallel <https://www.gnu.org/software/parallel/>`_ and write:

.. code-block:: txt

  parallel json_interface ::: input/mass*.json

We will get the following output in Terminal:

.. code-block:: txt

  FAILURE: propagation of file "~/tudatBundle/matlabInterface/Examples/InputOutput/MultiplePropagations/input/mass200.json" terminated with errors.
  Minimum step size exceeded.
  Error, propagation terminated at t=141762.97465721529, returning propagation data up to current time
  SUCCESS: propagation of file "~/tudatBundle/matlabInterface/Examples/InputOutput/MultiplePropagations/input/mass300.json" terminated with no errors.
  SUCCESS: propagation of file "~/tudatBundle/matlabInterface/Examples/InputOutput/MultiplePropagations/input/mass400.json" terminated with no errors.
  SUCCESS: propagation of file "~/tudatBundle/matlabInterface/Examples/InputOutput/MultiplePropagations/input/mass500.json" terminated with no errors.
  SUCCESS: propagation of file "~/tudatBundle/matlabInterface/Examples/InputOutput/MultiplePropagations/input/mass600.json" terminated with no errors.
  SUCCESS: propagation of file "~/tudatBundle/matlabInterface/Examples/InputOutput/MultiplePropagations/input/mass700.json" terminated with no errors.
  SUCCESS: propagation of file "~/tudatBundle/matlabInterface/Examples/InputOutput/MultiplePropagations/input/mass800.json" terminated with no errors.

As you can see, the first propagation fails, so it is terminated before reaching the termination condition. This was done on purpose to showcase the propagation failure handling capabilities of the MATLAB and JSON interfaces.

The directory tree looks like this after generating the results:

.. code-block:: txt

  MultiplePropagations
  | 
  | input
  |     |
  |     | mass200.json
  |     | ...
  |     | mass800.json
  |     | shated.json
  |      
  | output 
  |      |
  |      | mass200.txt
  |      | ...
  |      | mass800.txt
  |      
  | generateInput.m
  | processOutput.m

The script :literal:`processOutput.m` can be used to generate a plot using the results file. This script uses the :literal:`[results,failure] = import.results` function from the MATLAB Interface. In this way, it can be known whether one of the generated output files corresponds to a propagation that terminated before reaching the termination conditions, due to a propagation error. This is possible because the file :literal:`mass200.txt` has a header indicating that a failure took place:

.. code-block:: txt
  :caption: :class:`matlabInterface/Examples/InputOutput/MultiplePropagations/output/mass200.txt`

  FAILURE
    0                 4667617.77873914  2333808.88936957  4042275.57154398  -5494.76025754408  2747.38012877204  4758.6019707383   230000.00000205 
    20                4556442.06802713  2388104.54009974  4136317.49487142  -5622.29470046153  2681.9360346687   4645.15892336866  230001.407005601
  ...
  
In order to disable this functionality (i.e. tagging output files with the :literal:`FAILURE` header when the propagation fails), one can write:

.. code-block:: matlab

  simulation.options.tagOutputFilesIfPropagationFails = false;
