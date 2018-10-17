.. _downloadTudatBundle:

Download Tudat Bundle
=====================
The next step in the installation guide is to download the Tudat Bundle. This is a pre-configured setup that will allow us to get the libraries up and running in no time. The Tudat Bundle contains the following libraries:

- Tudat
- Boost (automatically downloaded and configured for your system when first building Tudat)
- Eigen
- cspice
- jsoncpp
- nrlmsise-00
- pagmo2

Our recommended approach to download Tudat bundle is by using a terminal. Alternatively, a git-client with a GUI (such as SmartGit) can be used. 

Cloning Tudat
~~~~~~~~~~~~~

For Linux or Mac, the regular terminal can be used directly. Windows users are recommended to use the :literal:`tudat_shell.bat`, located in :literal:`tudatBundle/external/tools`, `here <https://github.com/Tudat/tudatBundle/blob/master/external/tools/tudat_shell.bat>`_. Copy the contents into a plain text file, and save it as :literal:`tudat_shell.bat`.  Running this script will open a Windows terminal, with the git, make and cmake commands available.

After opening a terminal, navigate to the directpry where you want to download tudatBundle. You can use the :literal:`ls` command to get a listing of the current directory. To navigate 'up' one level (from ::literal`..../directoryA/directoryB/` to ::literal`..../directoryA/`), use the :literal:`cd ..` command. For the inverse (from ::literal`..../directoryA/` to ::literal`..../directoryA/directoryB/`) use :literal:`cd directoryB`.

  .. warning:: If you are working on a Windows system, it is highly recommended that you place the Tudat Bundle directly in your C:\ (or D:\, E:\, ... ) drive, as failing to do so has been know to cause compilation errors resulting from excesive path length. Mac and Linux users are nor constrained in where the place the bundle.

In the directory where tudatBundle is to be downloaded, enter the command::


    git clone https://github.com/tudat/tudatBundle.git

This will download the basic folder structure of the tudat bundle. Once the process is complete, use the following commands::

    cd tudatBundle
    git submodule update --init --recursive

Note that the last command in particular can take some time, depending on the speed of your internet connection.

  .. note:: This whole process is known to be decidedly slower on Windows than on Mac/Linux. 

Create a GitHub account
~~~~~~~~~~~~~~~~~~~~~~~
For those users who want to push their Tudat code to Github (including students taking the AE4866 and AE4868 courses), a Github account is required. Create an account `here: <https://github.com/join?source=header-home>`_ if you haven't done so yet. 

A GitHub account is completely free and provides you with hosting for your software projects. It also allows you to fork other repositories to create your own improved version of a piece of code.

As a student you can upgrade your free account to a personal account (normally $7/month). With a personal account you can have unlimited private repositories. Private repositories can be very helpful when developing new unreleased programs and/or writing your thesis. Applying for the free upgrade is easily done via `this <https://education.github.com/discount_requests/new>`_ link and tick boxes: Student and Individual account. Then enter the requested information.
