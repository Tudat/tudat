.. _setupDevelopmentEnvironmentLinux:

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

**Step 3: Check Settings**
   Once Qt Creator is installed, you will need to verify that the various compilation settings have been defined correctly. Make sure to check ``Qt Creator Kits`` on the :ref:`verifyKitsAndCMake` page.

You are now ready to download the Tudat bundle, use the navigation below to go to the next step.
