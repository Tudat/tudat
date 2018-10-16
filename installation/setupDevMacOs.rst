.. _setupDevelopmentEnvironmentMacOs:

Install on Mac OS X
-------------------
This section of the guide refers to the installation of the IDE under Mac OS. These instructions should work for all versions of Mac OS X. You may need to give the Administrator password when installing some programs, depending on your user settings.

CMake
~~~~~
So, we are now going to walk through the process of installing CMake on your computer. We use CMake to manage software project(s), and to instruct your compiler how to generate and link libraries.

**Step 1: Download CMake**
    Download CMake from the `CMake website <https://cmake.org/download/>`_.

**Step 2: Launch Apple disk image**
    Once you have successfully downloaded the ``.dmg`` file, simply double-click on it.

**Step 3: Install CMake application**
    To complete the installation, drag and drop CMake onto the Applications folder.

**Step 4: Launch CMake**
    Launch CMake from Launchpad or Spotlight. The interface that pops up is *only* used to test whether CMake has been correctly installed, and to perform step 5 of this guide. You do not need to fill anything in.

    .. error:: If you get the warning: "'CMake' can't be opened because it is from an unidentified developer.", please make sure to enable "Allow apps dowloaded from: Anywhere" under the "System Preferences" > "Security & Privacy".
 
**Step 5: Add CMake to path**
    From the "Tools" menu select "How to Install For Command Line Use". From the dialog that pops up, note the ``cmake-gui path``, this may be required later. Open a terminal by executing Cmd+Space, typing terminal and confirming with Enter. Type::

        sudo mkdir -p /usr/local/bin
        sudo /Applications/CMake.app/Contents/bin/cmake-gui --install=/usr/local/bin

    Hopefully no errors occured. 

    .. error:: If an error occurs, check whether the path ``/Applications/CMake.app/Contents/bin/cmake-gui`` corresponds to the ``cmake-gui path`` that you noted down earlier. If not, change the above command so that the two match.

    Now, Verify that it has been correctly installed to PATH by executing::

        cmake --version

    .. error:: If cmake can't be found, even after succesfully installing CMake for command-line use, you first need to verify that the symbolic links were properly made, execute::

            ln -al /usr/local/bin/cmake

     You should see something like::

        lrwxr-xr-x 1 root wheel 42 5 Feb 10:49 /usr/local/bin/cmake -> /Applications/CMake.app/Contents/bin/cmake

     This means that ``/usr/local/bin`` is not added to the list of paths. Do this now manually by editing ``/etc/paths``::

        sudo pico /etc/paths

     and adding ``/usr/local/bin`` on a new line. Exit (Cmd+X) and save (Y) your changes. Close the terminal and restart the terminal::

        cmake --version

     should now give you the desired result, showing some details on your CMake installation.

git
~~~
We are now going to walk through the process of installing git on your computer. We use git to download the software of the project(s) and to make sure that you can always be up-to-date on the latest modifications to the code.

**Step 1: Download git**
    The install process for git is very similar to that for CMake: go to the `git website <https://git-scm.com/downloads/>`_ to download the installer ``.dmg``.

**Step 2: Install git**
    Run the ``.dmg`` and open the enclosed ``.pkg`` to install git on your system. Step through the installation and provide you administrator password when prompted.

XCode
~~~~~
For the compilation of Tudat and its libraries, XCode (or the command-line tools for XCode) version 7.3 or newer is required. You can upgrade XCode through the AppStore or by downloading a new version, to replace the old one.

**Step 1: Download XCode**
    Download XCode (command-line tools alone suffices) from the `Apple developer downloads <https://developer.apple.com/download/more/>`_ (ADC account and Apple ID required) or through the Mac App Store. Note installing only the "Command Line Tools OSX 10.XX for Xcode 7.X" offers a significant reduction in size (download size of 157MB vs 4.7GB).

**Step 2: Install XCode**
Open the downloaded ``.dmg`` and execute the enclosed ``.pkg`` to start the installation. Complete the installation.

Qt Creator
~~~~~~~~~~
**Step 1: Download Qt Creator**
    Download QtCreator from the `Qt website <https://www1.qt.io/download-open-source/>`_.

**Step 2: Execute the installer**
    Open the downloaded ``.dmg`` and execute the enclosed installer.

**Step 3: Skip account**
    You can safely skip logging into your Qt account. Press "Skip" and "Next", the online installer will prepare the sources. Click continue. The installer will now prepare the installation (this will take a short while). You might be prompted by an Xcode warning, even though you have Xcode or the Xcode command-line tools installed. If you encounter this, click away the warning by pressing "Ok", three times, the installation will continue as normal.

**Step 4: Choose a location**
    Specify your preferred installation directory (or leave it at default).

**Step 5: Select components**
    Click "Continue" until you get to the "Select Components" step. Here you get the option to select which parts of the Qt SDK you wish to install, shown below. Only QtCreator (default, can not be unchecked) from the Tools section is necessary. Finish the installation.

**Step 6: Check Settings**
   Once Qt Creator is installed, you will need to verify that the various compilation settings have been defined correctly. Make sure to check ``Qt Creator Kits`` on the :ref:`verifyKitsAndCMake` page.

