.. _setupDevelopmentEnvironment:

1. Setup Development Environment
================================

The first step is to set up the Integrated Development Environment (IDE) that will be used to write, compile and run your code. We recommend to use Qt Creator as your IDE due to its multi-platform capability and its popularity among students and staff. Note that the compiler and the precise installation steps will vary depending on your OS, so please refer to the section that is relevant to you. 

Install on Linux
----------------
This section of the guide refers to the installation of the IDE under a Linux-based OS.

CMake & git
~~~~~~~~~~~
This step depends on your distribution and more specifically its package manager, the most common are treated below. Please refer to the documentation of your distribution.
Fire up a terminal and enter the following command appropriate for your distribution:

**Debian-based distributions (Debian, Ubuntu, Mint, etc.)**::

    sudo apt-get install build-essential cmake
    sudo apt-get install git

**RedHat-based distributions (Red Hat, Fedora, CentOS, etc.)**::

    sudo yum install make gcc gcc-c++ kernel-devel cmake 
    sudo yum install git

**OpenSUSE**::
 
    sudo zypper install --no-recommends gcc gcc-c++ make cmake
    sudo zypper install --no-recommends git

**ArchLinux**:: 

    sudo pacman -S base-devel cmake
    sudo pacman -S git

Note that the compiler and make tools are maintained by the manager. Please make sure your versions are up to date, see your package manager documentation for more information. In the very unlikely scenario that your manager does not have the previously mentioned packages or they are out of date, you can always install these packages manually or through experimental / user repositories. In this case refer to the documentation of GCC and CMake.

Qt Creator
~~~~~~~~~~
**Step 1: Download Qt SDK**
    Go to http://www.qt.io/download-open-source/ and download the Qt SDK for Linux, either the 32- or 64-bit install, depending on your system and preferences.

**Step 2: Install Qt Creator**
    Next, head to the download folder and make the installer executable by typing the following in a terminal::

        chmod +x qt-unified-linux-x64-[version]-online.run

    Then launch the installer using::

        ./qt-unified-linux-x64-[version]-online.run.

    Please follow the instructions given throughout the installer. Note that you can safely "Skip" creating an account. When you get to the component selection, you can deselect all optional components, save for Qt Creator. Once the installation is finished, you can run Qt Creator by typing /home/user/Qt/Tools/QtCreator/bin/qtcreator from the command-line or use your favourite application launcher.

You are now ready to download the Tudat bundle, use the navigation below to go to the next step.

Install on Windows
------------------
This section of the guide refers to the installation of the IDE under Windows. These instructions should work for all versions of Windows from XP to Windows 10.

CMake
~~~~~
So, we are now going to walk through the process of installing CMake on your computer, after which we will install Qt Creator, which comes bundled with the compiler and related tools. We use CMake to manage software project(s), and to instruct your compiler how to generate and link libraries.

**Step 1: Download CMake**
    To download CMake, go to the CMake website and look for the Windows Installer ``cmake-[version]-win32-x86.msi`` at the top of table for Binary distributions. Click the link, and choose Run or Save (if available) in order to get the installer on your computer.

**Step 2: Install CMake**
    If your installer does not run automatically once the download has been completed, go to your Downloads directory and start the installer by double-clicking on it. Click on Run if your operating system tells you it does not recognize the publisher. You now get to the welcome screen for the CMake installer. You will have to press Next and accept the License Agreement (if you choose) and Next again. Now you will get to the Install Options. Make sure to select either Add CMake to the system PATH for *all users* (preferred) or *current user* (alternatively).

    .. note:: In rare cases, the installer will be unable to add CMake to your system path, because it is already too long (mind you the PATH list, not the CMake path, which is limited at 4095 characters). In this case choose "Add CMake to the PATH for current user" instead or try to clean up your system path manually.

**Step 3: Complete the installation**
    Finish the remainder of the installation. Leave the rest default (location and shortcuts).

Qt Creator & MinGW
~~~~~~~~~~~~~~~~~~
Now we are going to walk through the process of installing Qt Creator on your computer. This will also install the MinGW, a Minimalist GNU Environmentn for Windows, containing (amongst others) GNU C++ compilers and Make programs.

**Step 1: Download the installer**
    To begin the installation, use the following link to download the installer for the Qt SDK: http://www.qt.io/download-open-source/ As before it is best to choose to Run the instsaller directly.

**Step 2: Install Qt Creator**
    If you have selected Run, the installer should start automatically. If it does not start automatically, go to your downloads directory and start the installer by double-clicking on it.

    You will now see a welcome screen. You can "Skip" creating an account. Next the online installer will do some checks, this may take a few moments. You should get the option to specify the install directory, see the step below.

**Step 3: Set installation directory**
   Do not install in Program Files, instead install Qt in your C:/ directory (or D:/ ...). Do not install it in the Program Files (x86) directory (or any directory with spaces), since this is know to cause issues later on.

**Step 4: Select custom installation components**
    Click Next until you get to the Select Components step. Here you get the option to select which parts of the Qt SDK you wish to install, shown below. Only Qt Creator (default, can not be unchecked) and MinGW (latest version) from the Tools section are necessary.

    .. note:: It is necessary to have MinGW 4.9.1 (or greater), otherwise you will not be able to compile the Tudat libraries.

**Step 5: Complete the installation**
    Click Next to go to the License Agreement, where you find the license attached to the use of Qt. Scroll through it if you wish to know what you are accepting. Agree to the conditions and continue to the next step. You can choose to have shortcuts installed or not. Now you can click Install to start the installation. QtCreator and MinGW will be installed on your computer.

**Step 6: Configure Qt Creator**
    Open Qt Creator and go to Tools/Options. Under the Build & Run section, choose the Kits tab and select Desktop. Please verify the following:
    - The compiler is set to MinGW >= 4.9.1.
    - The CMake tool is present.
    - The primary generator is MinGW and secondary generator is CodeBlocks.

The next step is to download the Tudat bundle. Click next to go there.

Install on Mac OS X
-------------------
This section of the guide refers to the installation of the IDE under Mac OS. These instructions should work for all versions of Mac OS X. You may need to give the Administrator password when installing some programs, depending on your user settings.

CMake
~~~~~
So, we are now going to walk through the process of installing CMake on your computer. We use CMake to manage software project(s), and to instruct your compiler how to generate and link libraries.

**Step 1: Download CMake**
    Download CMake from the CMake website.

**Step 2: Launch Apple disk image**
    Once you have successfully downloaded the ``.dmg`` file, simply double-click on it.

**Step 3: Install CMake application**
    To complete the installation you must solve the puzzle (no worries it's kindergarten stuff): drag and drop CMake onto the Applications folder.

**Step 4: Launch CMake**
    Launch CMake from Launchpad or Spotlight.

    .. note:: If you get the warning: "'CMake' can't be opened because it is from an unidentified developer.", please make sure to enable "Allow apps dowloaded from: Anywhere" under the "System Preferences" > "Security & Privacy".
 
**Step 5: Add CMake to path**
    From the "Tools" menu select "How to Install For Command Line Use". From the dialog that pops up, note the ``cmake-gui path``. Open a terminal by executing Cmd+Space, typing terminal and confirming with Enter. Type::

        sudo mkdir -p /usr/local/bin
        sudo /Applications/CMake.app/Contents/bin/cmake-gui --install=/usr/local/bin

    Hopefully no errors occured. Verify that it has been correctly installed to PATH by executing::

        cmake --version

    .. note:: If cmake can't be found, even after succesfully installing CMake for command-line use, you first need to verify that the symbolic links were properly made, execute::

            ln -al /usr/local/bin/cmake

     You should see something like::

        lrwxr-xr-x 1 root wheel 42 5 Feb 10:49 /usr/local/bin/cmake -> /Applications/CMake.app/Contents/bin/cmake

     This means that ``/usr/local/bin`` is not added to the list of paths. Do this now manually by editing ``/etc/paths``::

        sudo pico /etc/paths

     and adding ``/usr/local/bin`` on a new line. Exit (Cmd+X) and save (Y) your changes. Close the terminal and restart the terminal::

        cmake --version

     should now give you the desired result.

git
~~~
We are now going to walk through the process of installing git on your computer. We use git to download the software of the project(s) and to make sure that you can always be up-to-date on the latest modifications to the code.

**Step 1: Download git**
    The install process for git is very similar to that for CMake: go to the git website to download the installer ``.dmg``.

**Step 2: Install git**
    Run the ``.dmg`` and open the enclosed ``.pkg`` to install git on your system. Step through the installation and provide you administrator password when prompted.

XCode
~~~~~
For the compilation of Tudat and its libraries, XCode (or the command-line tools for XCode) version 7.3 or newer is required. You can upgrade XCode through the AppStore or by downloading a new version, to replace the old one.

**Step 1: Download XCode**
    Download XCode (command-line tools alone suffices) from the Apple developer downloads (ADC account and Apple ID required) or through the Mac App Store. Note installing only the "Command Line Tools OSX 10.XX for Xcode 7.X" offers a significant reduction in size (download size of 157MB vs 4.7GB).

**Step 2: Install XCode**
Open the downloaded ``.dmg`` and execute the enclosed ``.pkg`` to start the installation. Complete the installation.

Qt Creator
~~~~~~~~~~
**Step 1: Download Qt Creator**
    Download QtCreator from the Qt website.

**Step 2: Execute the installer**
    Open the downloaded ``.dmg`` and execute the enclosed installer.

**Step 3: Skip account**
    You can safely skip logging into your Qt account. Press "Skip" and "Next", the online installer will prepare the sources. Click continue. The installer will now prepare the installation (this will take a short while). You might be prompted by an Xcode warning, even though you have Xcode or the Xcode command-line tools installed. If you encounter this, click away the warning by pressing "Ok", three times, the installation will continue as normal.

**Step 4: Choose a location**
    Specify your preferred installation directory (or leave it at default).

**Step 5: Select components**
    Click "Continue" until you get to the "Select Components" step. Here you get the option to select which parts of the Qt SDK you wish to install, shown below. Only QtCreator (default, can not be unchecked) from the Tools section is necessary.

**Step 6: Complete installation**

