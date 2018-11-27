.. _debuggingOpeningCMake:

Installation errors
===================

This page contains often encountered errors during installation of Tudat. If your problem is not listed, please visit the projects `issue page <https://github.com/Tudat/tudat/issues>`_ and make sure to check the already closed issues. If this still doesn't resolve your issue, please open an `issue <https://github.com/Tudat/tudat/issues/new>`_ on Github yourself.

Opening CMake
~~~~~~~~~~~~~

- Compiler not found::

      -- The C compiler identification is unknown
      -- The CXX compiler identification is unknown
      CMake Error at CMakeLists.txt:14 (project):
       The CMAKE_C_COMPILER:

        cl

       is not a full path and was not found in the PATH.

       Tell CMake where to find the compiler by setting either the environment
       variable "CC" or the CMake cache entry CMAKE_C_COMPILER to the full path to
       the compiler, or to the compiler name if it is in the PATH. 


      CMake Error at CMakeLists.txt:14 (project):
       The CMAKE_CXX_COMPILER:
   
        cl

       is not a full path and was not found in the PATH

       Tell CMake where to find the compiler by setting either the environment
       variable "CXX" or the CMake cache entry CMAKE_CXX_COMPILER to the full path
       to the compiler, or to the compiler name if it is in the PATH. 
   
 This indicates that your C and or C++ compiler are not set correctly. This screenshot indicates an error for the C compiler, for the C++ compiler case :literal:`CMAKE_CXX_COMPILER` instead of :literal:`CMAKE_C_COMPILER` will give the error. This can be fixed by following the steps in QtCreator Kits in :ref:`verifyKitsAndCMake` 
   
- Cmake tool incorrect::

    Running "/usr/bin/cmake /home/dominic/Software/tudatBundle -GNinja -DCMAKE_CXX_COMPILER:STRING=/usr/bin/clang++ -DCMAKE_C_COMPILER:STRING=/usr/bin/clang '-DCMAKE_PREFIX_PATH:STRING=%  ]
    {Qt:QT_INSTALL_PREFIX}' -DQT_QMAKE_EXECUTABLE:STRING=" in /home/dominic/Software/build-tudatBundle-Desktop-Default.
    -- Configuring incomplete, errors occurred!
    See also "/home/dominic/Software/build-tudatBundle-Desktop-Default/CMakeFiles/CMakeOutput.log".
    See also "/home/dominic/Software/build-tudatBundle-Desktop-Default/CMakeFiles/CMakeError.log".
    CMake Error: CMake was unable to find a build program corresponding to "Ninja".  CMAKE_MAKE_PROGRAM is not set.  You probably need to select a different build tool.   
    *** cmake process exited with exit code 1.

 This indicates that your CMake generator is not set correctly. This can be fixed by following the steps in QtCreator Kits in :ref:`verifyKitsAndCMake`. 

- Nothing happens (no output, no error message) when opening CMake

   Your CMake executable is set wrongly. This can be fixed by following the steps in QtCreator CMake binary in :ref:`verifyKitsAndCMake` 

- Cmake GUI opens:

   .. figure:: images/cmakeGui.png

   This occurs if the cmake executable is set to cmake-gui , instead of cmake. This can be fixed by following the steps in QtCreator CMake binary in :ref:`verifyKitsAndCMake`. 

.. _verifyKitsAndCMake:

Verify Build & Run options
**************************

- QtCreator CMake binary

   1. Open QtCreator, go to ``Preferences/Options``, select the ``Build & Run`` section and switch to the ``CMake`` tab
   2. It is very important that the QtCreator is pointed to the correct cmake binary. CMake ships with multiple binaries and often the wrong one is selected. The correct binaries are:

      - ``C:\Program Files (x86)\CMake\bin\cmake.exe`` instead of ``C:\Program Files (x86)\CMake\bin\cmake-gui.exe``
      - ``/usr/bin/cmake or /usr/local/bin/cmake`` instead of ``/usr/local/bin/cmake-gui``

   3. ``/Applications/CMake.app/Contents/bin/cmake`` instead of ``/Applications/CMake.app/Contents/MacOS/CMake`` or ``/Applications/CMake.app``.

   .. note:: Make a screenshot of the CMake tab if problems persist.

- QtCreator Kits

   1. Open QtCreator, go to ``Preferences/Options``, select the ``Kits`` tab (on older versions of Qt, choose the ``Build & Run`` -> ``Kit``):
   2. Verify settings:

      - Generator:
         - Unix Makefiles on Linux or Mac OS X
         - MinGW Makefiles on Windows
      - Extra generator: CodeBlocks
      - Device type: Desktop
      - Compiler C/C++
         - MinGW >= 7.3.0 
         - GCC or Clang on Linux or Mac OS X
         - C compiler should be non-empty, it is needed for certain libraries.

   It should look like:

   .. figure:: images/compilerCheck.png

   This screenshot is for Windows, on which the MinGW compiler is used. For your system, the C/C++ compiler may be a version of GCC or Clang. Make Sure that both compiler (indicated by the red box) are both set, and set to the same compiler version.


   .. note:: Make a screenshot of the Kits tab if problems persist.


.. _debuggingdDownloadingCompilingBoost:

Downloading/compiling boost
~~~~~~~~~~~~~~~~~~~~~~~~~~~

- Boost and compiler are incompatible
   Not all versions of Boost are compatible with each system. Please refer to the compatibility table below and pick a different version (by commenting out/in lines in ``tudatBundle/CMakeLists.txt``). The installed version can be found as described in :ref:`verifyInstallationCmakeAndCompiler`. 

   +------------------+--------------+--------------+--------------+
   |**Compiler/Boost**|**Boost 1.53**|**Boost 1.57**|**Boost 1.60**|
   +------------------+--------------+--------------+--------------+
   |MinGW 4.9.1       |      ✓       |       ?      |        ?     |
   +------------------+--------------+--------------+--------------+
   |MinGW 4.9.2       |      ✓       |       ✓      |        ✓     |
   +------------------+--------------+--------------+--------------+
   |MinGW 5.3         |      ✗       |       ✓      |        ✓     |
   +------------------+--------------+--------------+--------------+


- Boost download has failed

   The following error might sometimes occur when trying to download Boost::

      Downloading boost 1.64.0 to C:/tudatBundle/boost/build
      [download 100% complete]
      CMake Error at external/CMake/add_boost.cmake:214 (file):
        file DOWNLOAD HASH mismatch
         
   This error is likely to be due to an issue originating from the Boost server itself. It will usually be solved by waiting some time and re-running the configuration process (which will make a new attempt to download Boost).
   
   If it does not solve the error, then please check your anti-virus settings which might block the Boost download on Windows machines.

- Boost build failed
   - Go to the ``tudatBundle/boost/stage/lib`` folder and verify all the libraries you require are present.
   - Make note of all files in this folder.
   - Go to the ``tudatBundle/boost/boost`` folder and locate ``version.hpp`` and verify with the compatibility table above.
   - Copy this file along with your report.

   If Boost still fails, go to your build directory and locate the following four files::

      build-*/boost_1_XX_*/build_bootstrap.log
      build-*/boost_1_XX_*/build_b2.log
      build-*/boost_1_XX_*/cmake-config.jam
      build-*/boost_1_XX_*/project-config.jam

   Copy all four files along with your report and create an `issue <https://github.com/Tudat/tudat/issues/new>`_ on the Github project page.


- :literal:`file COPY cannot find boost_.../stage`

   This error was fixed by removing the ``~`` from the ``TEMP`` and ``TMP`` directory. See `github issue 259 <https://github.com/Tudat/tudat/issues/259>`_ for more details on the issue.

- System ``boost`` library used on linux

   When boost is already installed on the system your compiler might prefer the system ``boost`` over the ``boost`` in your tudatBundle as discussed in `this issue <https://github.com/Tudat/tudat/issues/203>`_. This could lead to incompatibility issues. 
   
.. _debuggingCompilationLinkingCode:

Compilation/linking of code
~~~~~~~~~~~~~~~~~~~~~~~~~~~

- Errors during compilation of Tudat on Windows

   Some error message could occur during the compilation of Tudat on Windows. However, if Qt is continuing compilation there is no need to worry.

- Many warnings (yellow triangle with exclemation mark) during compilation.

   These warnings can be safely ignored. To mute warnings and only display errors please check :ref:`qtBasics`.

- ``Undefined reference to ..../libsofa... or ..../libsofa...`` 
   This may occur if you are using multi-core compilation. Save the error output, and start the compilation again, if the same error occurs again, open an `issue <https://github.com/Tudat/tudat/issues/new>`_ on Github to report the issue.

- ``Undefined reference to ..../libtudat...``
   This indicates that the required Tudat libraries cannot be found when compiling. If this happens for the Tudat libraries, copy the compile output to a text file, and open an `issue <https://github.com/Tudat/tudat/issues/new>`_ on Github. If it is your own program, first check if you have added the required link libraries to your CMake file.

- ``Undefined reference to ......libboost/``
   This indicates that no compatible version of boost can be found when compiling the code. Most likely boost was not compiled correctly due to an incompatibility with your compiler. Check :ref:`debuggingdDownloadingCompilingBoost`, for compatibility between your compiler and boost version. If this does not resolve your problem, copy the compile output to a text file, retrieve the files listed under boost build failed, and open an `issue <https://github.com/Tudat/tudat/issues/new>`_ on Github.

- ``Error "out of memory allocating XXXX bytes"``
   This indicates that your compiler is using too much RAM, and your system cannot allocate it. First, copy your full compile output (tab at bottom of Qt Creator) to a text file. Then, change the :literal:`COMPILE_HIGH_ACCURACY_ESTIMATION_TESTS` CMake argument to OFF, and recompile. *Whether this fixes the error or not, open a Github issue*. This problem should have been corrected, and any occurence should be communicated. 

- ``Undefined reference to  `WinMain@16'``::

   Linking CXX executable ..... C:/PROGRA2/Qt/Tools/MINGW41/bin/../lib/gcc/i686-w64-mingw32/4.9.2/../../../../i686-w64-mingw32/lib/../lib/libmingw32.a(lib32_libmingw32_a-crt0_c.o):crt0_c.c:  
   (.text.startup+0x39): undefined reference to `WinMain@16'. 
     
  This error (or something similar) can occur (on Windows) if your compiler is in a directory containing a space (and possibly other non-standard character). Make sure that Qt and your compiler are installed in a directory like ``C:/Qt``, ``C:/mingw``, etc. Avoid the ``C:/Program Files`` directory.

- ``libbacktrace could not find executable to open``
   This error is due to multi-core compilation on Windows. Restarting the compile process fixes the issue (multiple restarts could be required), or compile with a singly thread.

- ``unlink .../libcspice.a: Permission denied``::

   mingw32-make.exe[2]: *** Deleting file 'C:/tudatBundle/cspice/lib/libcspice.a'
   mingw32-make.exe[2]: unlink: C:/tudatBundle/cspice/lib/libcspice.a: Permission denied

  This error may occur on Windows machines due to anti-virus settings blocking the Tudat compilation. You might want to check those settings and modify them, so that your anti-virus does not interfere with the building of Tudat. This should prevent the issue from happening. Please be aware that the modification of your anti-virus settings is your own responsibility.


.. _debuggingFailedUnitTests:

Failed unit tests
~~~~~~~~~~~~~~~~~~
There is a possibility of one or more unit tests failing. Usually, there is no cause for alarm, as this just means that your computer is rounding some variables a bit differently. To be sure, `open an issue on Github <https://github.com/Tudat/tudat/issues/new>`_. In this issue, attach the file ``LastTest.log``, which should be in the ``/Testing/Temporary/`` directory in your build folder. In the issue description and title, note that it concerns failed unit test(s) and mention your operating system. We'll get back to you with a fix for the failure ASAP.

