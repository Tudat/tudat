.. _configureTudatLibraries:

Compile and Test Tudat Libraries
================================

To make use of the Tudat libraries, you need to compile them. Below we will guide you through this process. In case if any issues, open a Github issue explaining your problem, stating your operating system, compiler, and any output of the process you have that may help us to solve the issue.

    .. note:: When attaching output, always put it into a text file, not a screenshot, making sure to copy all information in the ``Compile Output`` (if the error occurs during compilation) or the ``Application Output`` (if the error occurs when running an application).

**Step 1: Open Qt Creator**
    Launch the Qt Creator application that you installed. If any windows pop-up offering to help you get started, launch tutorials etc., hit Cancel. This should bring you into the editor.

**Step 2: Open the project**
    The next step is to open the CMake project of the tudatBundle. Click on ``Open File or Project...`` from the File drop-down menu. Navigate to where you extracted your Tudat Bundle, and navigate to the ``tudatBundle`` folder. Within this directory, you will see a file called ``CMakeLists.txt``. This is the main project file for any CMake project. Click on ``Open``, after selecting the CMakeLists.txt file.

    .. note:: Please note that you can safely ignore any git-related errors/warnings that Qt Creator throws. Example: Cannot run "git rev-parse --git-dir" in "C:\tudatBundle".

**Step 3: Condigure project**
    You will now get a 'Configure Project' screen. Leave all settings to default, and click ``Configure Project``. 

    The process of configuring the Tudat project and the required libraries will now be started. You will see output generated in the ``General messages`` box at the bottom of your screen, that will look something similar to::

        -- The C compiler identification is GNU 4.8.4
        -- The CXX compiler identification is GNU 4.8.4
        -- Check for working C compiler: /usr/bin/cc
        -- Check for working C compiler: /usr/bin/cc -- works
        -- Detecting C compiler ABI info
        -- Detecting C compiler ABI info - done
        -- Detecting C compile features
        -- Detecting C compile features - done
        -- Check for working CXX compiler: /usr/bin/c++
        -- Check for working CXX compiler: /usr/bin/c++ -- works
        -- Detecting CXX compiler ABI info
        -- Detecting CXX compiler ABI info - done
        -- Detecting CXX compile features
        -- Detecting CXX compile features - done
        -- /home/dominicdirkx/Software/tudatClean/tudatBundle/tudat/Tudat/External/CMake/
        -- /home/dominicdirkx/Software/tudatClean/tudatBundle
        -- BOOST: Using gnu.
        -- Downloading boost 1.60.0 to /home/dominicdirkx/Software/tudatClean/build-tudatBundle-Desktop-Default
        -- [download 0% complete]
        -- .......
        -- [download 100% complete]
        -- Extracting boost 1.60.0 to /home/dominicdirkx/Software/tudatClean/build-tudatBundle-Desktop-Default/boost_unzip
        -- Building b2 (bjam)
        -- ./bootstrap.sh;--with-toolset=gcc
        -- Build boost (note that this may take a while, please sit back)
        -- ./b2;link=static;threading=multi;runtime-link=shared;--build-dir=Build;stage;-d+2;--hash;--ignore-site-config;variant=release;cxxflags=-fPIC;cxxflags=-std=c++11;--layout=tagged;toolset=gcc;-sNO_BZIP2=1;--with-filesystem;--with-system;--with-thread;--with-regex;--with-date_time;--with-test
        -- Building CSpice from within TudatBundle.
        -- WARNING: building release version!
        -- JsonCpp Version: 1.6.5
        -- Building NRLMSISE00 from within TudatBundle.
        -- WARNING: building release version!
        -- Building Tudat from within TudatBundle.
        -- Tudat Relative path (wrt to project): /tudat/Tudat
        -- WARNING: building release version!
        -- Using gnucxx compiler.
        -- Performing Test CXX_SUPPORTS_CXX11
        -- Performing Test CXX_SUPPORTS_CXX11 - Success
        -- Found Eigen3: /usr/include/eigen3 (Required is at least version "2.91.0")
        -- Boost version: 1.60.0
        -- Found the following Boost libraries:
        -- date_time
        -- system
        -- unit_test_framework
        -- filesystem
        -- regex
        -- SPICE disabled!
        -- NRLMSISE-00 disabled!
        -- Building SatellitePropagatorExamples from within TudatBundle.
        -- Relative path (wrt to project): /tudatExampleApplications/satellitePropagatorExamples/SatellitePropagatorExamples
        -- WARNING: building release version!
        -- Using gnucxx compiler.
        -- Boost version: 1.60.0
        -- Found the following Boost libraries:
        -- thread
        -- date_time
        -- system
        -- unit_test_framework
        -- filesystem
        -- regex
        -- Found Tudat: /home/dominicdirkx/Software/tudatClean/tudatBundle/tudat/Tudat/.. (Required is at least version "2.0")
        -- Building SpiceAndJSON from within TudatBundle.
        -- Relative path (wrt to project): /tudatExampleApplications/libraryExamples/SpiceAndJSON
        -- WARNING: building release version!
        -- Using gnucxx compiler.
        -- Boost version: 1.60.0
        -- Found the following Boost libraries:
        -- thread
        -- date_time
        -- system
        -- unit_test_framework
        -- filesystem
        -- regex
        -- Relative path to Tudat found: /tudat/Tudat
        -- SPICE_LIBRARIES: cspice
        -- Found SPICE: /home/dominicdirkx/Software/tudatClean/tudatBundle/cspice/include/../..
        -- JSONCPP_LIBRARIES: jsoncpp
        -- Found JSONCPP: /home/dominicdirkx/Software/tudatClean/tudatBundle/jsoncpp/include/json/../../include
        -- Building TemplateApplication from within TudatBundle.
        -- Relative path (wrt to project): /tudatExampleApplications/templateApplication/TemplateApplication
        -- WARNING: building release version!
        -- Using gnucxx compiler.
        -- Boost version: 1.60.0
        -- Found the following Boost libraries:
        -- thread
        -- date_time
        -- system
        -- unit_test_framework
        -- filesystem
        -- regex
        -- Configuring done
        -- Generating done
        -- Build files have been written to: /home/dominicdirkx/Software/tudatClean/build-tudatBundle-Desktop-Default

    Depending on your system, boost may or may not be downloaded and compiled by CMake (it typically is). Depending on the speed of your computer and internet connection, this may take anywhere from several to 30 minutes. You can safely **ignore CMake warnings** about unused variables, specifically manually-specified variables were not used by the project, and warnings about relative paths, such as the one below::

      CMake Warning (dev) in tudat/Tudat/Astrodynamics/OrbitDetermination/RotationalDynamicsPartials/CMakeLists.txt:
        Policy CMP0081 is not set: Relative paths not allowed in LINK_DIRECTORIES
        target property.  Run "cmake --help-policy CMP0081" for policy details.
        Use the cmake_policy command to set the policy and suppress this warning.

        Found relative path while evaluating link directories of
        "tudat_torque_partials":

          "SOFA_BASE_PATH-NOTFOUND/../lib"

      This warning is for project developers.  Use -Wno-dev to suppress it.

    In case an **error occurs** during this portion of the installation, copy the full contents of the ``General Messages`` tab from Qt (bottom of screen) into a text file and post this with your Github issue.

**Step 4: Updating the Settings (optional)**
   Before building the libraries, you can modify some of the CMake settings to suit your needs. If you are not sure what your needs are (yet), leave all settings as they are, and proceed to the following step. If you are installing Tudat for the AE4868 course, do not modify the settings. To change the CMake settings, go to Projects->Build->CMake, see screenshot below (note that it may look slightly different, depending on your Qt version/operating system):

   .. figure:: images/cmakeSettings.png

   Below, we give some examples of changes that you may wish to make:

      * If you do not plan on using the estimation (variational equations propagation, observation models, acceleration partials *etc.*), you can switch ``BUILD_WITH_ESTIMATION_TOOLS`` to ``OFF``.
      * If you plan on using extended precision (*e.g.* more than 16 significant digits) for either state or time representation, you should switch ``BUILD_EXTENDED_PRECISION_PROPAGATION_TOOLS`` to ``ON``.
      * If you plan to use the Pagmo toolbox for optimization, set ``USE_PAGMO`` to ``ON``. Note that this will also trigger the compilation several example applications on how to use Pagmo with Tudat.

   .. note:: You can modify the CMake settings at any later point in time, but this may require a rebuild of a significant part of the libraries.

**Step 5: Build the libraries**
    Now all that remains to be done is to build the libraries. Typically, Tudat is compiled using a single core on your system. The compiler can be instructed to use multiple cores for compilation, check the FAQ :ref:`faqCompilationInstallation` for details. 

    .. note:: When encountering a compilation error during multi-core compilation, try reinintializing the compilation (clicking the hammer again). If the same errors occurs, and you wish to open an issue, rerun with a SINGLE thread, and post the output that this produces. Multi-core compilation output can often be garbled, and difficult to interpret.

    To compile all the libraries, simply click on the "hammer" build icon at the bottom-left of your screen (or use the menu ``Build`` at the top and select ``Build all``). You will see a ``Compile Output`` console window pop-up, showing the status of the build process, as the compiler walks through all the project files, and generates the libraries that we need. The entire build process could take anywhere from 15 minutes (Linux/Mac modern workstation; 12 threads) to 3-8 hours (Windows; single core), depending on the specifications of your computer. Have patience! It will all work out in the end. Once the building is complete, you're done! You have now successfully built Tudat and all required libraries on your computer.
    
**Step 6: Running the unit tests**
   For each part of the code in Tudat, we have written unit tests, which are included in the repository. Before moving on with using Tudat, you should run all the unit tests to ensure that your installation is functioning as it should. To run all unit tests, go to the project tab, and go to the ``Build Steps`` block. Write "test" in the ``Tool Arguments`` (may be called ``Additional Arguments``) line, as shown below.

   .. figure:: images/testSettings.png

   Now, go back to your code by clicking on the ``Edit`` tab, and click the ``Compile`` (hammer) button on the lower left. In the ``Compile Output`` console window at the bottom of your screen, you should see all the unit tests being run, with output as follows::

      15:15:48: Running steps for project TudatBundle...
      15:15:48: Starting: "/usr/bin/make" test
      Running tests...
      Test project /home/dominicdirkx/Software/tudat/build-tudatBundle-Desktop-Default 
      Start   1: sofa-test
      1/249 Test   #1: sofa-test ............................................................   Passed    0.01 sec
            Start   2: test_AerodynamicMomentAndAerodynamicForce
      2/249 Test   #2: test_AerodynamicMomentAndAerodynamicForce ............................   Passed    3.06 sec
            Start   3: test_AerodynamicsNamespace
      3/249 Test   #3: test_AerodynamicsNamespace ...........................................   Passed    0.00 sec
            Start   4: test_AerodynamicCoefficientGenerator
      4/249 Test   #4: test_AerodynamicCoefficientGenerator .................................   Passed    0.03 sec
            Start   5: test_ExponentialAtmosphere
      5/249 Test   #5: test_ExponentialAtmosphere ...........................................   Passed    0.00 sec
            Start   6: test_CustomConstantTemperatureAtmosphere
      6/249 Test   #6: test_CustomConstantTemperatureAtmosphere .............................   Passed    0.00 sec
            Start   7: test_TabulatedAtmosphere
      7/249 Test   #7: test_TabulatedAtmosphere .............................................   Passed   26.81 sec
            Start   8: test_TabulatedAerodynamicCoefficients
      8/249 Test   #8: test_TabulatedAerodynamicCoefficients ................................   Passed    1.37 sec
      ...
      ...
      ...
      243/249 Test #243: test_JsonInterfaceTermination ........................................   Passed    0.02 sec
              Start 244: test_JsonInterfaceThrust
      244/249 Test #244: test_JsonInterfaceThrust .............................................   Passed    0.01 sec
              Start 245: test_JsonInterfaceTorque
      245/249 Test #245: test_JsonInterfaceTorque .............................................   Passed    0.00 sec
              Start 246: test_JsonInterfaceVariable
      246/249 Test #246: test_JsonInterfaceVariable ...........................................   Passed    0.01 sec
              Start 247: test_JsonInterfaceObservation
      247/249 Test #247: test_JsonInterfaceObservation ........................................   Passed    0.09 sec
              Start 248: test_JsonInterfaceParameter
      248/249 Test #248: test_JsonInterfaceParameter ..........................................   Passed    0.05 sec
              Start 249: test_JsonInterfaceSimulationSingleSatelliteVariational
      249/249 Test #249: test_JsonInterfaceSimulationSingleSatelliteVariational ...............   Passed    0.09 sec

      100% tests passed, 0 tests failed out of 249
      Total Test time (real) = 623.61 sec
      15:16:48: The process "/usr/bin/make" exited normally.

   Depending on your exact compilation settings, and the speed of your system, running the unit tests may take anywhere from several to 30 minutes. Also, depending on your settings, and version of the code, you will run a different number of unit tests.

   If the output ends with ``100% tests passed, 0 tests failed``, all is well and you do not need to take any further action. After running the unit tests, make sure to remove the 'test' text that you've typed in the project tab. If any tests fail the reader is refered to :ref:`debuggingFailedUnitTests`. 

    .. note:: After running the unit tests, make sure to remove the ``test`` text that you've typed in the ``Build Steps``, Qt will not compile the code as long as it is there.

   So, welcome to Tudat. You are now ready to run one of the many example applications that came bundled with Tudat, and get started on setting up your won application. The applications are explained in detail in the tutorials at Tutorials and Documentation. The next and last (optional) part explains you how to set-up a new application or add existing ones to your Tudat Bundle.

**Step 7: Run An Application**

   Before moving on to using Tudat for the example applications (or your own application), modify the ``Build Settings`` to build only the current application. Now that the unit tests are built and run, there is no need to recompile everything everytime. Only the portions relevant for the specific application under consideration need to be compiled. See the screenshot below for the option to tick that enforces this behaviour. 

    .. note:: When (re)running the unit tests, always first recompile the code with the targets set to ``all``. If there is no ``all`` tick box, uncheck all other boxes, it will automatically revert to ``all`` as default.

   .. figure:: images/runSettings.png

   For your convenience, we have shipped some example applications for you to play around with. The structure of these applications is discussed in detail in the :ref:`walkthroughsIndex`. 

   To select a specific application to run, click on the ``Build and Run Settings`` (computer) icon and select your application. For starters, select :literal:`application_SingleSatellitePropagator`. The exact contents and results of this simulation are show in the :ref:`walkthroughsUnperturbedEarthOrbitingSatellite` tutorial By clicking the ``Run`` button (play icon in bottom left), the code will be compiled and the selected application will be executed. The output of your application is displayed in the ``Application Output`` box at the bottom of your screen. In addition, a folder 'SimulationOutput' will have been created in your :literal:`/tudatExampleApplications/satellitePropagatorExamples/SatellitePropagatorExamples/` directory, containing the propagation output.
   
   .. code-block:: cpp

      Starting .../tudatBundle/tudatExampleApplications/satellitePropagatorExamples/bin/applications/application_SingleSatellitePropagator...
      Single Earth-Orbiting Satellite Example.
      The initial position vector of Asterix is [km]:
      7037.48
      3238.06
      2150.72
      The initial velocity vector of Asterix is [km/s]:
      -1.46566
      -0.0409584
      6.6228
      After 86400 seconds, the position vector of Asterix is [km]:
      -4560.45
      -1438.32
       5973.99
      And the velocity vector of Asterix is [km/s]:
      -4.55021   
      -2.41254
      -4.95063
      .../tudatBundle/tudatExampleApplications/satellitePropagatorExamples/bin/applications/application_SingleSatellitePropagator exited with code 0

      	
