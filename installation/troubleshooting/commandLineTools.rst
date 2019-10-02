.. _troubleShootingInstallationCommandLineTools:

Command line tools
==================

This page contains some command line tools which are usefull for troubleshooting your installation. First some command line basics are described. Then the verification of your CMake and compiler installations are described.

Troubleshooting tools
~~~~~~~~~~~~~~~~~~~~~
Below are several tools and tricks that are essential for troubleshooting your installation.

**Terminal**
    
   The terminal (also console) is king in troubleshooting and the great unifier across systems. Some benefits:

       - Although the output can be harsh, plentiful and unformatted, but is often very complete.
       - Usually additional features and options can be accessed from the commandline that can not be gained through a regular GUI.
       - Commands are often very similar across different platforms and version of software.

   To open the terminal emulator on your system follow the instruction below:

        Linux::
        
            Super + T

        Windows::

            Win + R

        macOS::

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

        macOS::

            cd /Users/John
            cd tudatBundle

    This should be everything you need to get to the right place.

    .. note:: Note that for Windows you can't cd to another drive (e.g. from C:\\ to D:\\). To switch drives, simply enter D:. Also note that Windows directories are seperated by a backslash (\\) and UNIX-like systems follow POSIX's forward-slash seperator ("/"). However, Windows is POSIX compliant, meaning that cd \Users\John and cd /Users/John are both correct on Windows.

**Creating log files**

    In several checks of the troubleshooting guide you are asked to note down the output of several commands. Most terminal emulators support copy pasting, but even this can be very bothersome for long output, rather you can direct the output from a command to a file::

        gcc -v                       (all output to screen)
        gcc -v > gcc_log.txt         (normal output to logfile, errors to screen) 
        gcc -v > gcc_log.txt 2>&1    (normal and error output to logfile)

    Unfortunately while creating a log file you can not see the output as it is directed to the file instead of the screen. On Linux and macOS there is a command called tee, which does both::

        gcc -v | tee gcc_log.txt         (normal output to file and screen, errors screen only)
        gcc -v  2>&1 | tee gcc_log.txt   (normal and error output to logfile and screen, FREFERRED)

    On Windows you can run the command twice (once without logging and once with) or checkout the logfile using a text editor.

**Command not found**

    A number of causes can result in a "command not found" error while trying to execute something. The most common are (in no specific order):

    - The command was mistyped: even getting the case wrong makes a difference on some UNIX-like systems (e.g. Cmake and cmake are not the same).
    - Special character were not escaped: this problem starts with improper file and directory names, but My Program will be interpreted as two seperate command. Solutions are to escape the offending characters (the space in this case) My\ Program or quote all parts "My Program"
    - The program is simply not installed: the obvious solution would be to install the program.
    - The executables are installed, but not are not added to the PATH environment variable: the next section will focus on this problem in particular.

**Modifying the PATH variable**

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


.. _verifyInstallationCmakeAndCompiler: 

Verify installation of CMake and Compiler
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Please complete all checks and copy the output (even if a check seems succesfull). Some checks will provide additional tests and solutions in case the check initially fails. Please complete these instructions as well and redo the original test!


**Check CMake installation**
   Open your terminal (emulator) and type::

      cmake --version

   Report the full output of this command. If your output looks like the example below, you can assume that CMake has been installed correctly::

      cmake version 3.4.3  
      CMake suite maintained and supported by Kitware (kitware.com/cmake).

**Check Compiler installation**
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
    - If so try uninstall and reinstalling the latest version, please follow the Installation Guide to make sure you don't miss any steps.
    - In case of CMake version CMake 3.4 or earlier, it needs to be fully uninstalled before installing 3.5 or later.

   **Solution B:** Make sure the software is added to the path.

    - On both Windows and macOS, some software is not added to the path by default.
    - Usually this is an option inside the installation (not for CMake under macOS), reinstalling and selecting this option should fix the issue.
    - If not you need to manually add the directory that contains the executables to you path, see instructions below.


