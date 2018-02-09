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

- .. figure:: images/compilerNotFound.png

This indicates that your C and or C++ compiler are not set correctly. This screenshot indicates an error for the C compiler, for the C++ compiler case :literal:`CMAKE_CXX_COMPILER` instead of :literal:`CMAKE_C_COMPILER` will give the error. To fix, open Qt Creator and go to Tools/Options. Under the Build & Run section, choose the Kits tab and select Desktop. Please verify the following:

.. figure:: ../../installation/images/compilerCheck.png


This screenshot is for Windows, on which the MinGW compiler is used. For your system, the C/C++ compiler may be a version of GCC or Clang. Make Sure that both compiler (indicated by the red box) are both set, and set to the same compiler version (see :ref:`troubleshootingConfig` Step 2.6).

- Cmake tool incorrect

.. code-block:: cpp

   Running "/usr/bin/cmake /home/dominic/Software/tudatBundle -GNinja -DCMAKE_CXX_COMPILER:STRING=/usr/bin/clang++ -DCMAKE_C_COMPILER:STRING=/usr/bin/clang '-DCMAKE_PREFIX_PATH:STRING=%  ]
   {Qt:QT_INSTALL_PREFIX}' -DQT_QMAKE_EXECUTABLE:STRING=" in /home/dominic/Software/build-tudatBundle-Desktop-Default.
   -- Configuring incomplete, errors occurred!
   See also "/home/dominic/Software/build-tudatBundle-Desktop-Default/CMakeFiles/CMakeOutput.log".
   See also "/home/dominic/Software/build-tudatBundle-Desktop-Default/CMakeFiles/CMakeError.log".
   CMake Error: CMake was unable to find a build program corresponding to "Ninja".  CMAKE_MAKE_PROGRAM is not set.  You probably need to select a different build tool.   
   *** cmake process exited with exit code 1.

This indicates that your CMake generator is not set correctly. See :ref:`troubleshootingConfig` Step 2.6 for how to fix this.


- If nothing happens (no output, no error message) when opening CMake, your CMake executable is set wrongly. See :ref:`troubleshootingConfig` Step 2.6 for how to fix this.

- Cmake GUI opens:

.. figure:: images/cmakeGui.png

This occurs if the cmake executable is set to cmake-gui , instead of cmake. See :ref:`troubleshootingConfig` Step 2.6 for how to fix this.


.. _debuggingdDownloadingCompilingBoost:

Downloading/compiling boost
***************************

- Errors occur when building boost. See :ref:`troubleshootingConfig` Step 2.1. If this does not correct the issue, open a Github issue, with the files indicates in :ref:`troubleshootingConfig` Step 2.2.

- :literal:`file COPY cannot find boost_.../stage`

   This error was fixed by removing the ~ from the TEMP and TMP directory. See `github issue 259 <https://github.com/Tudat/tudat/issues/259>`_ for more details on the iss
   
.. _debuggingCompilationLinkingCode:

Compilation/linking of code
***************************

-   Undefined reference to ..../libsofa... or ..../libsofa... This may occur if you are using multi-core compilation. Save the error output, and start the compilation again, if the same error occurs again, open an issue on Github to report the issue.

-   Undefined reference to ..../libtudat... This indicates that the required Tudat libraries cannot be found when compiling. If this happens for the Tudat libraries, copy the compile output to a text file, and open an issue on Github. If it is your own program, first check if you have added the required link libraries to your CMake file.

-   Undefined reference to ......libboost/. This indicates that no compatible version of boost can be found when compiling the code. Most likely boost was not compiled correctly due to an incompatibility with your compiler. Check :ref:`troubleshootingConfig`, step 2.1 for compatibility between your compiler and boost version. If this does not resolve your problem, copy the compile output to a text file, retrieve the files listed under step 2.2 of :ref:`troubleshootingConfig`, and open an issue on Github.

-   Error "out of memory allocating XXXX bytes" during compilation (with XXX some number). This indicates that your compiler is using too much RAM, and your system cannot allocate it. First, copy your full compile output (tab at bottom of Qt Creator) to a text file. Then, change the :literal:`COMPILE_HIGH_ACCURACY_ESTIMATION_TESTS` CMake argument to OFF, and recompile. *Whether this fixes the error or not, open a Github issue*. This problem should have been corrected, and any occurence should be communicated. When 


-   .. code-block:: cpp

      Linking CXX executable ..... C:/PROGRA2/Qt/Tools/MINGW41/bin/../lib/gcc/i686-w64-mingw32/4.9.2/../../../../i686-w64-mingw32/lib/../lib/libmingw32.a(lib32_libmingw32_a-crt0_c.o):crt0_c.c:  
      (.text.startup+0x39): undefined reference to `WinMain@16'. 
     
   This error (or something similar) can occur (on Windows) if your compiler is in a directory containing a space (and possibly other non-standard character). Make sure that Qt and your compiler are installed in a directory like C:/Qt, C:/mingw, etc. Avoid the C:/Program Files directory.

.. _debuggingFailedUnitTests:

Failed unit tests
*****************

See :ref:`configureTudatLibraries`, end of step 6.

Frequently made coding mistakes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This section sums op some of the most made mistakes leading to failures either in runtime or during compilation. 

   - Accessing vector elements out of vector bounds::
      
      Eigen::Vector2d vectorOfLengthTwo(2.0, 3.0);
      double intermediateResult = vectorOfLengthTwo(2); 

   This will actually work, but the value obtained ``intermediateResult`` is the value located just next to the ``vectorOfLengthTwo`` in your memory and could be anything. 






