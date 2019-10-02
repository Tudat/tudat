.. _downloadTudatBundleSmartgitUser:

Download using SmartGit (User)
------------------------------
This guide is intended for Tudat users. A more elaborate guide for developers on how to setup SmartGit and clone Tudat Bundle go to the :ref:`downloadTudatBundleSmartgitDeveloper` guide.

Setup SmartGit
~~~~~~~~~~~~~~~~~
**Step 1.1: Download SmartGit**
    Download SmartGit for your operating system. Note that for Linux the package might be available through your package manager, see Linux/Debian instructions. Note that SmartGit comes bundled with JRE (Java Runtime Environment) on Windows and macOS. On Linux you need to install ``jre-openjdk``.

**Step 1.2: Install SmartGit**
    The installation should be straight-forward and the default installation is fine.

**Step 1.3: Launch SmartGit**
    At the end of the installation you will be prompted to launch SmartGit. A shortcut should have been added to your program list.

**Step 1.4: Choose a license type**
    SmartGit is free for non-commercial use. On first launch SmartGit will ask you to accept EULA and choose a license type. Choose non-commercial if you are planning to use SmartGit non-commercially (like for Tudat). You will be asked to confirm your choice.

**Step 1.5: Fill in your user information**
    Please fill in your username and e-mail address. This information will be used to attribute you when you commit code changes.

**Step 1.6: Complete the configuration**
    The next two steps can be left as default. If you already have a GitHub account you can link this to SmartGit. This can however be easily done later too.

Add the Tudat Bundle repository
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**Step 2.1: Clone an existing repository**
    At first launch you will be welcomed and asked to create/clone/reopen a repository, choose ``Clone`` and click ``Ok``.

**Step 2.2: Set remote address**
    As the remote address fill in ``https://github.com/tudat/tudatBundle.git`` and choose ``Next``. The Clone URL can generally be found at the top of a GitHub page.

**Step 2.3: Include submodules**
    Make sure submodules are include (default).

**Step 2.4: Set local directory**
   Pick a good location for your local directory. End the path with tudatBundle or something similar. The files will be downloaded directly in the folder you choose.

   .. tip:: It is good practice to avoid long paths and names with special characters (preferably no spaces). Pick something short and sensible. Especially if you are working on a Windows system, it is recommended that you place the Tudat Bundle directly in your C:\ drive, to avoid problems with filepath lengths.

Congratulations! You have now downloaded the tudatBundle. You can now head to the :ref:`configureTudatLibraries` guide to configure the bundled libraries correctly and build them.

