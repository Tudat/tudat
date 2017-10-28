.. _guideToDebugging:

Guide To Debugging
==================
This page contains some valuable tips for debugging your application. These tips can help narrow down where your error occurs and thereby helps you solve many issues by yourself. 


Debuggers
~~~~~~~~~

Qt debugger
***********
This build in debugger is available on all platforms. It can be started from within Qt by selecting ``Debug -> Start Debugging -> Start Debugging``. 

In the debugger console the output of the debugger can be seen. This often gives information on where your problem arises. The following is an example of the output of the Qt built in debugger for an error occuring in the FlightConditions update function. 

.. figure:: images/QtDebugger.png

gdb debugger
************
On linux use can be made of the gdb debugger. It is installed from the terminal by using (for debian-based distributions)::

   $ sudo apt-get update
   $ sudo apt-get install gdb

Then navigate to the folder with your application usually something like::

   $ cd tudatBundle/tudatApplication/myApplicationsFolder/bin/applications/

Then run gdb with the following command::

   $ gdb ./myApplicationName

Followed by::

   (gdb) run

This will run the application and stops when it stumbles upon debugging symbols. If your application crashes the error can be traced back at this point by using::

   (gdb) backtrace

When applied to the same example as for the Qt build-in debugger, gdb provide the following output in the terminal:

.. literalinclude:: DebuggingOutput.txt


Searching Tudat project files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Users can search the complete Tudat project from within Qt. This can be usefull e.g. for finding the input and output types of functions, for finding the derived classes of settings or searching for error messages given during runtime. This can help locate the problem. 

Searching the complete project can be done by using ``ctrl+shift+f`` which opens the following:

.. figure:: images/searchCompleteProject.png

As an example the resulting output when searching for ``"Error in dynamics simulator, propagator settings not defined"`` is shown below:

.. figure:: images/searchResults.png

.. _debuggingListOfKnownErrors:

List of known errors
~~~~~~~~~~~~~~~~~~~~
Here some known errors are described, if available it also provide the link to the relevant github issue for more information and context. This section is under development and will be updated with new content in the future. Past issues can be found on the `github issue page <https://github.com/Tudat/tudat/issues>`_. Don't forget to check the closed issues.

.. _debuggingOpeningCMake:

Opening CMake
*************

.. _debuggingdDownloadingCompilingBoost:

Downloading/compiling boost
***************************

- :literal:`file COPY cannot find boost_.../stage`

   This error was fixed by removing the ~ from the TEMP and TMP directory. See `github issue 259 <https://github.com/Tudat/tudat/issues/259>`_ for more details on the issue.
   
.. _debuggingCompilationLinkingCode:

Compilation/linking of code
***************************

.. _debuggingFailedUnitTests:

Failed unit tests
*****************
 

Frequently made coding mistakes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This section sums op some of the most made mistakes leading to failures either in runtime or during compilation. 

   - Accessing vector elements out of vector bounds::
      
      Eigen::Vector2d vectorOfLengthTwo(2.0, 3.0);
      double intermediateResult = vectorOfLengthTwo(2); 

   This will actually work, but the value obtained ``intermediateResult`` is the value located just next to the ``vectorOfLengthTwo`` in your memory and could be anything. 






