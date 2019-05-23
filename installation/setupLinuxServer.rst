.. _setupLinuxServer:

Tudat on a remote Linux server
==============================

This section of the guide refers to the use of Tudat on a remote Linux server.  

Accessing the server
~~~~~~~~~~~~~~~~~~~~

This section aims at detailing how to access the server. Specific examples will be given with the TU Delft server (only accessible to TU Delft students and staff, please contact the server administrator to ask for server access). The global directions are however applicable to every server you may want to run Tudat on. 
Login to the server is usually done via SSH (secure shell). For Linux or Mac OS X environments, use the following commands from your terminal::

	ssh -Y <userID>@serverAdress (<userID>@eudoxos.tudelft.nl for the TU Delft server

.. note:: Adding :literal:`-Y` after the :literal:`ssh` command activates the interface that will allow you to open graphical software installed on the server on your local machine. 

From a Windows machine, you have to follow a slightly different procedure as you will need an SSH client  Xclient (for X11). One of the most commonly used SSH client is Putty, which you can download from https://www.putty.org. Once downloaded, run the installer (you can keep the default install options) and open it. Fill the tab 'session' as follows, with the adress of your server (HostName) and the port number you will use to access it, and then select the SSH option. The example below corresponds to the TU Delft server:

    .. figure:: images/puttySessionTab.PNG

Then, you have to go to the :literal:`Data` tab under :literal:`Connection` and fill in your user account.

    .. figure:: images/puttyDataTab.PNG

Finally, go to the :literal:`SSH` tab (still under :literal:`Connection`) and then into :literal:`X1`' and tick :literal:`Enable X11 forwarding`. Going back to the :literal:`Session` tab, enter a name for the connection and save. The configuration you have just set up is then saved. When you will later open Putty, simply double-clicking on the connection name you have chosen will directly initiate the connection. 

   .. note:: If you want to open graphical software from the server on your local machine, you will also need to install an Xclient (for X11) (VcXsrv or Exceed are possible solutions).

To move files back and forth from the server, you will need a FTP client. The recommended one is FileZilla (you can download it from https://filezilla-project.org, and run the installer keeping the default options). Open FileZilla, and enter the hostname of your server, port number (same as before), user account and associated password and click :literal:`Quickconnect` (in the blanks highlighted in blue in the figure below). For later connections, note that you can directly ckick on the small arrow next to the :literal:`Quickconnect` button (in red below) and find your connection details already saved in the drop-down menu.

   .. figure:: images/fileZilla.png

You can then simply drag-and-drop your files from the server (right column) to your local machine (left column) and vice-versa.

Installing and building Tudat Bundle
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Most of the steps described below are very similar to the install guide for Linux environments (:ref:`setupDevelopmentEnvironmentLinux`) but will be reminded here for the sake of clarity.

Once logged in your server account, navigate to the folder in which you want to download the Tudat code. To print the content of the folder you are in, use the :literal:`ls` command. To navigate from your current folder to the 'parent' one, use :literal:`cd ..` (so that if you are in :literal:`folder1/folder2/` and type :literal:`cd ..`, you will end up in :literal:`folder1/`). On the contrary, use :literal:`cd folderName` to do the opposite (enter :literal:`cd folder2` to get to :literal:`folder1/folder2/` from :literal:`folder1/`). Once you are in your chosen folder, clone the tudatBundle repository with the following command::

	git clone https://github.com/tudat/tudatBundle.git
	
When the tudatBundle download is completed, a folder :literal:`tudatBundle` has been created in your chosen directory (you can verify it by typing the :literal:`ls` command, as mentioned above). Go to this tudatBundle directory (:literal:`cd tudatBundle`) and enter the following command:: 

	git submodule update --init --recursive
	
This step aims at updating the required modules and might take some time.

Once all the previous steps are done, you might want to modify the default Tudat settings before building the whole tudat code. This can be done by opening the :literal:`CMakeLists.txt` file directly located in the tudatBundle folder. This file can be opened with any available editor (:literal:`vim` and :literal:`nano` are the most common ones (both are available on the TU Delft server)). Typing either :literal:`vim CMakeLists.txt` or :literal:`nano CMakeLists.txt` respectively will open the file in the terminal. Depending on the settings you want, you might set the different options (namely :literal:`USE_CSPICE`, :literal:`USE_JSON`, :literal:`USE_NRLMSISE00`, :literal:`NRLMSISE support`, :literal:`USE_SOFA`, :literal:`USE_PAGMO`, :literal:`USE_PYGMO`, :literal:`BUILD_WITH_ESTIMATION_TOOLS`) to either :literal:`ON` or :literal:`OFF`. As an example, the option :literal:`USE_PAGMO` is set to :literal:`OFF` by default but should be turned on if you are planning on using the PAGMO toolbox for optimisation. 

  .. warning:: As a general advice and if you do not need it in particular, please switch off the :literal:`USE_JSON` option on the TU Delft server. Indeed, it will most probably yield compilation issues because the version of the gcc compiler currently available on the server is not supported by the JSON interface. In case you really need the JSON interface, please contact the server administrator to discuss about installing more recent gcc versions.

If you have done any settings modification, and before building the Tudat code, make sure to check whether the file :literal:`CMakeCache.txt` exists in the tudatBundle directory. If so, please delete it (with the command :literal:`rm CMakeCache.txt`) to ensure that your settings changes will be effective when building the code.

The tudatBundle libraries are now ready to be built. Just entering the command :literal:`make` in your terminal will initiate the compilation (make sure you are in the tudatBundle folder before doing it). This will build the tudat code and might take a while (from several dozens of minutes up to several hours). When the process is over, Tudat has been successfully build in your server account! 
The only remaining step is to run all the unit tests to ensure Tudat is working properly. It can be done from the tudatBundle directory by typing the following command::

	cmake --build . --target all -- test
	
You should then be able to see the unit tests being run in your terminal, the output looking as follows::

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

Depending on your compilation settings, this step can take from several minutes to more than one hour (the number of unit tests also depends on your settings). If the output ends with :literal:`100% tests passed, 0 tests failed`, everything worked out and you do not need to take any further action. If any tests fail the reader is refered to :ref:`debuggingFailedUnitTests`. 


Running Tudat applications on the server
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You are now ready to play around with the sample applications already available in Tudat or to run your own applications on the server. The different applications can be accessed and run from the following folders, depending on the type of application.
	
	- If you want to re-run an unit test independently, go to :literal:`tudatBundle/tudat/bin/unit_tests/`.
	- If you want to run an example application, go to :literal:`tudatBundle/tudatExampleApplications/satellitePropagatorExamples/bin/applications/`.
	- If you want to run an PAGMO optimisation application, go to :literal:`tudatBundle/tudatExampleApplications/libraryExamples/bin/applications/`.

Once you are in the proper directory, use the 'ls' command to get the list of all the built executables. Enter :literal:`./` directly followed by the name of your targeted executable to run it (:literal:`./executable_name`). The executable outputs will then appear in your terminal. 

It is of course also possible to write your own applications and run them on the server. The guidelines to write your application are presented in :ref:`createNewApps`. As mentioned there, new applications are typically added to the :literal:`tudatBundle/tudatApplications/` folder. If so, then the executables created after building your new applications can be found in :literal:`tudatBundle/tudatApplications/bin/applications/` and run with the same :literal:`./` command which has been described above.
For more details about getting a new application from an existing github repository or creating a totally new one, the reader is referred to :ref:`createNewApps`.

When you modify either an application file or some parts of the tudat code, you need to redo the building process (with the command :literal:`make`, from the folder containing the files you have modified). This will automatically generate new executables corresponding to the updated version of the code.

As previously mentioned, the application outputs appear in the terminal when running the associated executable. However, it might be advantageous to run the application in the background to be allowed to log out from the server while keeping your application running. Several options are possible here:

	- If you have already run your application as described above (so with a simple :literal:`./executable_name` command) but want to put it in the background to be able to log out, then press :literal:`ctrl+Z` to pause it and enter the command :literal:`bg` to put it in the background. Then you can type :literal:`exit` to log out and your application will keep running.
	
	- You can also directly start the process in the background by using the command :literal:`./executable_name &`. However, if you do so, the application outputs will not be accessible. You can choose to store them in a log file so that you can still read them when the running is over. This can be done with the following command (:literal:`log_run_date` being the name of the log file)::
	
		./executable_name > log_run_date 2>&1 &   	
