.. _downloadTudatBundleCommandLine:

Download using Command-line
---------------------------
This part of the installation guide will allow you to download the Tudat Bundle by means of the command-line. This assumes that you have succesfully installed CMake and git as explained in the :ref:`setupDevelopmentEnvironment`.

**Step 1: Open terminal**
    See the note below to learn how to open the terminal depending on your OS.

**Step 2: Navigate to a desired parent directory and clone**::

    git clone https://github.com/tudat/tudatBundle.git

**Step 3: Enter directory, initiate and fetch submodules**::

    cd tudatBundle
    git submodule update --init --recursive

**Step 4: Add your remotes (Optional)**::

    git remote add myaccount https://github.com/myaccount/tudatBundle.git
    git fetch myaccount
    cd tudat
    git remote add myaccount https://github.com/myaccount/tudat.git
    git fetch myaccount

**Step 5: Navigate to tudatBundle and prepare build**
    In terminal, navigate to the bundle using::

        cd

    Create a new directory using::

        mkdir build

    Enter the new build directory::

        cd build

.. note:: Depending on your OS, opening the Terminal will require a different set of steps.

    Linux::

        Super + T
    
    Windows::
        
        Win + R
        cmd

    macOS::

        Cmd + Space
        terminal
