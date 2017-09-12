.. _troubleshooting:

Troubleshooting
==================
Welcome to the troubleshooting guide. You have most likely been asked to come here because you have a problem with your Tudat installation or its configuration. Hopefully the steps discussed here will help you get up and running!

Troubleshooting installation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Please complete all checks and copy the output (even if a check seems succesfull). Some checks will provide additional tests and solutions in case the check initially fails. Please complete these instructions as well and redo the original test!

**Step 1.1: Check CMake**
    Open your terminal emulator and type::

        cmake --version

    Report the full output of this command. If your output looks like the example below, you can assume that CMake has been installed correctly::

        cmake version 3.4.3  
        CMake suite maintained and supported by Kitware (kitware.com/cmake).

**Step 1.2: Check Compiler**
    Open your terminal and check the version of your C++ compiler using the command below. This command is followed by two example outputs. Both outputs are fine, but if you get compilation errors, please report the versions.

    Test Command 1::

        gcc -v
        mingw32-make -v

    Test Command 2::

        gcc -v
        make -v

    Test Command 1, Example Output 1 (GCC)::

         Using built-in specs.
         COLLECT_GCC=g++
         COLLECT_LTO_WRAPPER=/usr/lib/gcc/x86_64-unknown-linux-gnu/5.3.0/lto-wrapper
         Target: x86_64-unknown-linux-gnu
         Configured with: /build/gcc/src/gcc-5-20160209/configure --prefix=/usr [...]
         Thread model: posix
         gcc version 5.3.0 (GCC)


    Test Command 1. Example Output 2 (Clang)::

         clang version 3.7.1 (tags/RELEASE_371/final)
         Target: x86_64-unknown-linux-gnu
         Thread model: posix

    Test Command 2::

         GNU Make 4.1
         Built for x86_64-unknown-linux-gnu
         Copyright (C) 1988-2014 Free Software Foundation, Inc.
         License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>
         This is free software: you are free to change and redistribute it.
         There is NO WARRANTY, to the extent permitted by law.

    **Solution A:** Make sure the software is installed (and up-to-date).

    - Please make sure CMake has been installed on your system.
    - If so try uninstall and reinstalling the latest version, please follow the Installation Guide to make sure you don't miss any steps.
    - In case of CMake version CMake 3.4 or earlier, it needs to be fully uninstalled before installing 3.5 or later.

    **Solution B:** Make sure the software is added to the path.

    - On both Windows and OS X, some software is not added to the path by default.
    - Usually this is an option inside the installation (not for CMake under OS X), reinstalling and selecting this option should fix the issue.
    - If not you need to manually add the directory that contains the executables to you path, see instructions below.

Troubleshooting configuration
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The configuration of your project includes downloading the Tudat Bundle source and running CMake.

**Step 2.1: Boost compiler compatibility**

    Not all version of Boost are compatible with each system. Please refer to the compatibility table below and pick a different version (by commenting out/in lines in ``tudatBundle/CMakeLists.txt``).

    +------------------+--------------+--------------+--------------+
    |**Compiler/Boost**|**Boost 1.53**|**Boost 1.57**|**Boost 1.60**|
    +------------------+--------------+--------------+--------------+
    |MinGW 4.9.1       |      ✓       |       ?      |        ?     |
    +------------------+--------------+--------------+--------------+
    |MinGW 4.9.2       |      ✓       |       ✓      |        ✓     |
    +------------------+--------------+--------------+--------------+
    |MinGW 5.3         |      ✗       |       ✓      |        ✓     |
    +------------------+--------------+--------------+--------------+

**Step 2.2: Check if Boost build succesfully**

    - Go to the ``tudatBundle/boost/stage/lib`` folder and verify all the libraries are present.

    - Make note of all files in this folder.
    - Go to the ``tudatBundle/boost/boost`` folder and locate ``version.hpp``.
    - Copy this file along with your report.

    If Boost still fails, go to your build directory and locate the following four files::

        build-*/boost_1_XX_*/build_bootstrap.log
        build-*/boost_1_XX_*/build_b2.log
        build-*/boost_1_XX_*/cmake-config.jam
        build-*/boost_1_XX_*/project-config.jam

    Copy all four files along with your report.

**Step 2.3: Clean your CMake cache**

    If you have made changes to your configuration it is important to clean your cache.
    
    - Navigate to your build folder.
    - Delete ``CMakeCache.txt``.

**Step 2.4: Manually build Boost yourself**
    
    The top-level CMakeLists.txt of Tudat Bundle, downloads, extracts, configures, builds and installs Boost for you. Although this process is completely automated it can happen that it fails somewhere doing the former.

    .. note:: If the automated process fails, it is necessary to take note of where it fails please copy the output of CMake for specifics.

    1. Find the Boost version TudatBundle is trying to build.
        1. Open tudatBundle/CMakeLists.txt
        2. Look for the uncommented (without a # in front) instance of set(BoostVersion 1.XX.0).
    2. https://sourceforge.net/projects/boost/files/boost/
        1. Pick the version corresponding to your version. Do not select beta.
        2. It doesn't matter which archive type you select, generally pick .tar.bz2 for Linux and OS X and .zip for Windows.
    3. Unpack the folder somewhere, for instance /home/user/boost or c:\boost.
    4. Open terminal emulator and go to the Boost folder.
    5. Run bootstrap:
        1. ./bootstrap.sh --with-toolset=gcc
        2. .\bootstrap.bat gcc
    6. If successful, run bjam2:
        1. ./b2 toolset=gcc link=static threading=multi --build-dir=Build stage variant=release --layout=tagged cxxflags=-std=c++11 --with-filesystem --with-system --with-thread --with-regex --with-date_time --with-test
        2. .\b2.exe toolset=gcc link=static threading=multi --build-dir=Build stage variant=release --layout=tagged cxxflags=-std=c++11 --with-filesystem --with-system --with-thread --with-regex --with-date_time --with-test
    7. In case of errors try to identify if bjam fails for each module or only for select modules.
    8. Rerun the b2 command several times, each time with only one and a different --with-[module] argument.

    .. note::Make logs of bootstrap and b2 command

**Step 2.5: Check QtCreator CMake binary**

    1. To verify the toolchain is correctly configured, open QtCreator, go to Preferences/Options, select the Build & Run section and switch to the CMake tab
    2. It is very important that the QtCreator is pointed to the correct cmake binary. CMake ships with multiple binaries and often the wrong one is selected. The correct binaries are:
        1. C:\Program Files (x86)\CMake\bin\cmake.exe instead of C:\Program Files (x86)\CMake\bin\cmake-gui.exe
        2. /usr/bin/cmake or /usr/local/bin/cmake instead of /usr/local/bin/cmake-gui
    3. /Applications/CMake.app/Contents/bin/cmake instead of /Applications/CMake.app/Contents/MacOS/CMake or /Applications/CMake.app.

    .. note:: Make a screenshot of the CMake tab if problems persist.

**Step 2.6: Check QtCreator Kits**

    1. To further verify the toolchain is correctly configured, open QtCreator, go to Preferences/Options, select the Build & Run section and switch to the Kits tab
    2. Main generator:
        1. Unix Makefiles
        2. MinGW Makefiles
        3. Common mistakes are Ninja or NMake
        4. Extra generator: CodeBlocks
    3. Device type: Desktop
    4. Compiler C/C++:
        1. MinGW >= 4.9.2
        2. GCC or Clang on Linux or Mac OS X
        3. Also the C compiler should be non-empty, it is needed for certain libraries.

    .. note:: Make a screenshot of the Kits tab if problems persist.

Troubleshooting tools
~~~~~~~~~~~~~~~~~~~~~
Below are several tools and tricks that are essential for troubleshooting your installation.

**3.1: Terminal**
    The terminal (also console) is king in troubleshooting and the great unifier across systems. Some benefits:

        - Although the output can be harsh, plentiful and unformatted, but is often very complete.
        - Usually additional features and options can be accessed from the commandline that can not be gained through a regular GUI.
        - Commands are often very similar across different platforms and version of software.

    To open the terminal emulator on your system follow the instruction below:

        Linux::
        
            Super + T

        Windows::

            Win + R

        Mac OS X::

            Cmd + Space
            terminal

    To go from one folder to another you can use the cd (change directory) command:

        Linux::
        
            cd /home/John
            cd tudatBundle

        Windows::

            C:
            cd \Users\John
            cd tudatBundle

        Mac OS X::

            cd /Users/John
            cd tudatBundle

    This should be everything you need to get to the right place.

    .. note:: Note that for Windows you can't cd to another drive (e.g. from C:\\ to D:\\). To switch drives, simply enter D:. Also note that Windows directories are seperated by a backslash (\\) and UNIX-like systems follow POSIX's forward-slash seperator ("/"). However, Windows is POSIX compliant, meaning that cd \Users\John and cd /Users/John are both correct on Windows.

**Step 3.2: Creating log files**

    In several checks of the troubleshooting guide you are asked to note down the output of several commands. Most terminal emulators support copy pasting, but even this can be very bothersome for long output, rather you can direct the output from a command to a file::

        gcc -v                       (all output to screen)
        gcc -v > gcc_log.txt         (normal output to logfile, errors to screen) 
        gcc -v > gcc_log.txt 2>&1    (normal and error output to logfile)

    Unfortunately while creating a log file you can not see the output as it is directed to the file instead of the screen. On Linux and OS X there is a command called tee, which does both::

        gcc -v | tee gcc_log.txt         (normal output to file and screen, errors screen only)
        gcc -v  2>&1 | tee gcc_log.txt   (normal and error output to logfile and screen, FREFERRED)

    On Windows you can run the command twice (once without logging and once with) or checkout the logfile using a text editor.

**Step 3.3 Command not found**

    A number of causes can result in a "command not found" error while trying to execute something. The most common are (in no specific order):

    - The command was mistyped: even getting the case wrong makes a diffirence on some UNIX-like systems (e.g. Cmake and cmake are not the same).
    - Special character were not escaped: this problem starts with improper file and directory names, but My Program will be interpreted as two seperate command. Solutions are to escape the offending characters (the space in this case) My\ Program or quote all parts "My Program"
    - The program is simply not installed: the obvious solution would be to install the program.
    - The executable are installed, but not are not added to the PATH environment variable: the next section will focus on this problem in particular.

**Step 3.4: Modifying the PATH variable**

    All systems look for matching programs into all the directories mentioned in the PATH variable when a command is typed and executed. You can check which directories are added to your PATH using the following command::

        echo %Path%
        echo $PATH

    Each user has a user-specific PATH directories and ones that it inherits from system. Usually we recommend installing a new directory entry for the PATH variable for the system (a.k.a. for all users). Usually during the setup of the application that we want to use from the command line there are some installation options that we can set such that the installer takes care of this. If you haven't done so it will be easiest to go back to the installer for the program and check if such an options exists, before attempting to do so manually.

    .. warning:: This can break your system 
 
    Follow these three essential steps carefully before attempting:
 
        1. Make sure that the directory your want to add is not already in the list.
        2. Backup current PATH var before modifying (instructions included below).
        3. Quadruple check for typos in the directory name you're adding. Example below uses MyApp

    ::

        # Making a backup of the current path
        set path_backup %path%
        echo %path% > %userprofile%\path_backup.txt

        # Changing the local PATH variable 
        # NOTE (1) quotes around everything 
        #      (2) semi-colon ; between old ond new 
        #      (3) Percentage signs % around %path%
        setx path "%path%;C:\MyApp\bin\" 

        # Verify that the operation was succesful
        echo %path%

        # Otherwise restore from backup with the following command
        setx path %path_backup%

    ::

        # Making a backup of the current path
        cp /etc/paths ~/paths_backup.txt

        # Changing the local PATH variable
        echo "/Applications/MyApp.app/Contents/bin/" | sudo tee -a /etc/paths

        # Verify that the operation was succesful" 
        cat /etc/paths

        # Otherwise we can restore from backup with the command
        sudo cp ~/paths_backup.txt /etc/paths

