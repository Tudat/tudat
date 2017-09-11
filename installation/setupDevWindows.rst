.. _setupDevelopmentEnvironmentWindows:

Install on Windows
------------------
This section of the guide refers to the installation of the IDE under Windows. These instructions should work for all versions of Windows from XP to Windows 10.

CMake
~~~~~
So, we are now going to walk through the process of installing CMake on your computer, after which we will install Qt Creator, which comes bundled with the compiler and related tools. We use CMake to manage software project(s), and to instruct your compiler how to generate and link libraries.

**Step 1: Download CMake**
    To download CMake, go to the `CMake website <https://cmake.org/download/>`_ and look for the Windows Installer ``cmake-[version]-win32-x86.msi`` at the top of table for Binary distributions. Click the link, and choose Run or Save (if available) in order to get the installer on your computer.

**Step 2: Install CMake**
    If your installer does not run automatically once the download has been completed, go to your Downloads directory and start the installer by double-clicking on it. Click on Run if your operating system tells you it does not recognize the publisher. You now get to the welcome screen for the CMake installer. You will have to press Next and accept the License Agreement (if you choose) and Next again. Now you will get to the Install Options. Make sure to select either Add CMake to the system PATH for *all users* (preferred) or *current user* (alternatively).

    .. note:: In rare cases, the installer will be unable to add CMake to your system path, because it is already too long (mind you the PATH list, not the CMake path, which is limited at 4095 characters). In this case choose "Add CMake to the PATH for current user" instead or try to clean up your system path manually.

**Step 3: Complete the installation**
    Finish the remainder of the installation. Leave the rest default (location and shortcuts).

Qt Creator & MinGW
~~~~~~~~~~~~~~~~~~
Now we are going to walk through the process of installing Qt Creator on your computer. This will also install the MinGW, a Minimalist GNU Environmentn for Windows, containing (amongst others) GNU C++ compilers and Make programs.

**Step 1: Download the installer**
    To begin the installation, use the following link to download the installer for the Qt SDK: http://www.qt.io/download-open-source/ As before it is best to choose to Run the installer directly.

**Step 2: Install Qt Creator**
    If you have selected Run, the installer should start automatically. If it does not start automatically, go to your downloads directory and start the installer by double-clicking on it.

    You will now see a welcome screen. You can "Skip" creating an account. Next the online installer will do some checks, this may take a few moments. You should get the option to specify the install directory, see the step below.

**Step 3: Set installation directory**
   Do not install in Program Files, instead install Qt in your C:/ directory (or D:/ ...). Do not install it in the Program Files (x86) directory (or any directory with spaces), since this is know to cause issues later on.

**Step 4: Select custom installation components**
    Click Next until you get to the Select Components step. Here you get the option to select which parts of the Qt SDK you wish to install, shown below. Only Qt Creator (default, can not be unchecked) and MinGW (latest version) from the Tools section are necessary.

    .. warning:: It is necessary to have MinGW 4.9.1 (or greater), otherwise you will not be able to compile the Tudat libraries.

**Step 5: Complete the installation**
    Click Next to go to the License Agreement, where you find the license attached to the use of Qt. Scroll through it if you wish to know what you are accepting. Agree to the conditions and continue to the next step. You can choose to have shortcuts installed or not. Now you can click Install to start the installation. QtCreator and MinGW will be installed on your computer.

**Step 6: Configure Qt Creator**
    Open Qt Creator and go to Tools/Options. Under the Build & Run section, choose the Kits tab and select Desktop. Please verify the following:

    - The compiler is set to MinGW 4.9.1 (or newer).
    - The CMake tool is present.
    - The primary generator is MinGW and secondary generator is CodeBlocks.

The next step is to download the Tudat bundle. Click next to go there.
