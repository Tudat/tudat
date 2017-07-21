.. _configureTudatLibraries:

Configure Tudat Libraries
=========================

.. warning:: This step requires that you have CMake, GCC compiler, Qt Creator installed and the Tudat Bundle (incl. libraries) on your computer.

To make use of the Tudat libraries, you need to compile them. The steps that follow will take you through the process of compiling Tudat. ``Eigen`` is the linear algebra library we use, and does not have to be compiled. ``Boost`` is a collection of useful libraries that are pre-compiled by your system, or are compiled when first running Tudat.

**Step 1: Open Qt Creator**
    Launch the Qt Creator application that you installed. If any windows pop-up offering to help you get started, launch tutorials etc., hit Cancel. This should bring you into the editor.

**Step 2: Open the project**
    The next step is to open the CMake project that contains cspice, Tudat, etc. This will allow us to compile the libraries. Click on ``Open File`` or ``Project...`` from the File drop-down menu. Navigate to where you extracted your Tudat Bundle, and navigate to the ``tudaBundle`` folder. Within this directory, you will see a file called ``CMakeLists.txt``. This is the main project file for any CMake project. Click on ``Open``, after selecting the CMakeLists.txt file.

    .. note:: Please note that you can safely ignore any git-related errors/warnings that Qt Creator throws. Example: Cannot run "git rev-parse --git-dir" in "C:\Users\Me\tudatBundle".

**Step 3: Set build-directory location**
    You will now see the ``CMake Wizard``, which will guide you through the process of configuring your library compilation. Note that the look of this interface can be a bit different, depending on your system. Hit ``Continue``, or ``Configure``, depending on which option is shown shown.

**Step 4: Run CMake**
    This will now bring you to the ``Run CMake`` screen (note, in some cases this screen may be skipped altogether, don't worry). If you are on the Run CMake screen, click ``Run CMake``.

    .. note:: The process of configuring the Tudat project and the required libraries will now be started. You will see output generated in the wizard, or the ``General messages`` box at the bottom of your screen, that will look something like::

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

    Depending on your system, boost may or may not be downloaded and compiled by CMake. Depending on the speed of your computer and internet connection, this may take anywhere from several to 15 minutes. You can safely ignore CMake warnings about unused variables, specifically manually-specified variables were not used by the project. In case of any problems, please send an e-mail to tudat-AE@tudelft.nl

**Step 5: Build the libraries**
    Now all that remains to be done is to build the libraries. To do this, simply click on the "hammer" build icon at the bottom-left of your screen (or use the menu ``Build`` at the top and select ``Build all``). You will see a ``Compile Output`` console window pop-up, showing the status of the build process, as the compiler walks through all the project files, and generates the libraries that we need. The entire build process could take anywhere from 5 to 30 minutes, depending on the specifications of your computer. Have patience! It will all work out in the end :). Once the building is complete, you're done! You have now successfully built Tudat and all required libraries on your computer. In case of any problems, please contact tudat_ae@tudelft.nl.

**Step 6: Run Template Application**
    For your convenience, we have shipped some example applications for you to play around with. As the basis for your future applications, your Tudat Bundle is shipped with a magnificent template application. All that remains to be done is to run it. To select a specific application to run, click on the ``Build and Run Settings`` (computer) icon and select the ``application_HelloWorld`` application. By clicking the ``Run`` button (play icon in bottom left), the code will be compiled and the selected application will be executed. However, this will also recompile all off the applications in your current project. Assuming that you have made no changes to the code, this process should be quite quick, but can take up to several minutes on a Windows machine. To tell Qt Creator to only build a single executable, click the project tab on the left. Subsequently, click on ``Details`` under ``Build Steps``. You will see a list of all applications and static libraries in the project. Select the one(s) you want to compile. Note that all dependencies of a given application will automatically be compiled as well. Now go back to your coding window by hitting ``Edit``. Click the ``Run`` button again. The output of your application is displayed in the green box.

    .. tip:: The output in the greenbox will not match your results exactly, since the Template Application generates random floating-point numbers to provide an example of using the Boost libraries. SO DON'T WORRY SO LONG AS THE OUTPUT IS SIMILAR!.

Congratulations! You've built and run your first Tudat application :)

**Step 7: Running the unit tests**
For each part of the code in Tudat, we have written unit tests, which are included in the repository. Before moving on with using Tudat, you should run all the unit tests to ensure that your installation is functioning as it should. To run all unit tests, go back to the project tab, and again go to the ``Build Steps`` block. In this block, uncheck the ``application_HelloWorld`` from the previous part and write "test" in the ``Additional Arguments`` line, as shown below. After running the unit tests, make sure to remove the "test" text that you've typed in here, Qt will not compile the code as long as it is there. Now, go back to your code by clicking on the ``Edit`` tab, and click the ``Compile`` (hammer) button on the lower left. In the ``Compile Output`` console window at the bottom of your screen, you should see all the unit tests being run, with output as follows::

    15:15:48: Running steps for project TudatBundle...
    15:15:48: Starting: "/usr/bin/make" test
    Running tests...
    Test project /home/dominicdirkx/Software/tudat/build-tudatBundle-Desktop-Default
    Start 1: test_Sofa
    1/132 Test 1: test_Sofa ................................................ Passed 0.03 sec
    Start 2: test_AerodynamicMomentAndAerodynamicForce
    2/132 Test 2: test_AerodynamicMomentAndAerodynamicForce ................ Passed 0.22 sec
    Start 3: test_AerodynamicsNamespace
    3/132 Test 3: test_AerodynamicsNamespace ............................... Passed 0.00 sec
    Start 4: test_AerodynamicCoefficientGenerator
    4/132 Test 4: test_AerodynamicCoefficientGenerator ..................... Passed 0.03 sec
    Start 5: test_ExponentialAtmosphere
    5/132 Test 5: test_ExponentialAtmosphere ............................... Passed 0.00 sec
    Start 6: test_TabulatedAtmosphere
    6/132 Test 6: test_TabulatedAtmosphere ................................. Passed 0.04 sec
    Start 7: test_TabulatedAerodynamicCoefficients
    7/132 Test 7: test_TabulatedAerodynamicCoefficients .................... Passed 1.61 sec
    Start 8: test_NRLMSISE00Atmosphere
    8/132 Test 8: test_NRLMSISE00Atmosphere ................................ Passed 0.01 sec
    Start 9: test_AstrodynamicsFunctions
    9/132 Test 9: test_AstrodynamicsFunctions .............................. Passed 0.00 sec
    Start 10: test_OrbitalElementConversions
    ...
    ...
    ...
    130/132 Test 130: test_SpiceInterface ...................................... Passed 0.05 sec
    Start 131: test_EnvironmentSetup
    131/132 Test 131: test_EnvironmentSetup .................................... Passed 2.90 sec
    Start 132: test_AccelerationModelSetup
    132/132 Test 132: test_AccelerationModelSetup .............................. Passed 0.16 sec
    100% tests passed, 0 tests failed out of 132
    Total Test time (real) = 59.57 sec
    15:16:48: The process "/usr/bin/make" exited normally.
    15:16:48: Elapsed time: 01:00.

    If the output ends with ``100% tests passed, 0 tests failed``, all is well and you do not need to take any further action. After running the unit tests , make sure to remove the 'test' text that you've typed in here, Qt will not compile the code as long as it is there. There is a possibility of one or more unit tests failing, though. Usually, there is no cause for alarm, as this just means that your computer is rounding some variables a bit differently, so that the 15th or 16th digit is different from what we expect. Just to be sure, however, you should "open an issue on Github": https://github.com/Tudat/tudat/issues/new In this issue, attach the file ``LastTest.log``, which should be in the ``/Testing/Temporary/`` directory in your build folder (which you specified in Step 2). In the issue description and title, note that it concerns failed unit test(s) and mention your operating system. We'll get back to you with a fix for the failure ASAP.

So, welcome to the Tudat universe :). The fun has just started though. You are now ready to run one of the many example applications that came bundled with Tudat, and this time it involves real simulations. The applications are explained in detail in the tutorials at Tutorials and Documentation. The next and last (optional) part explains you how to set-up a new application or add existing ones to your Tudat Bundle.
