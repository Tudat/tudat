.. _downloadTudatBundleSmartgitDeveloper:

Download using SmartGit (Developer)
-----------------------------------

.. tip:: Before following this guide, you are highly advised to go through the :ref:`githubBasics` section to get up to speed with the usage of GitHub.

Setup SmartGit
~~~~~~~~~~~~~~~~~
**Step 1.1: Download SmartGit**
    Download SmartGit for your operating system `here <http://www.syntevo.com/smartgit/>`_ . Note that for Linux the package might be available through your package manager, see Linux/Debian instructions. Note that SmartGit comes bundled with JRE (java runtime environment) on Windows and macOS. On Linux you need to install ``jre-openjdk``.

**Step 1.2: Install SmartGit**
    The installation should be straight-forward and the default installation is fine.

**Step 1.3: Launch SmartGit**
    At the end of the installation you will be prompted to launch SmartGit. A shortcut should have been added to your program list.

**Step 1.4: Choose a license type**
    SmartGit is free for non-commercial use. On first launch SmartGit will ask you to accept EULA and choose a license type. Choose non-commercial if you are planning to use SmartGit non-commercially (like for Tudat). You will be asked to confirm your choice.

**Step 1.5: Fill in your user information**
    Please fill in your username and e-mail address. This information will be used to attribute you when you commit code changes.

**Step 1.6: Complete the configuration**
    The next two steps can be left as default. If you already have a GitHub account you can link this to SmartGit. This can however be easily done later too.


Create a GitHub account
~~~~~~~~~~~~~~~~~~~~~~~~~~
`Create a GitHub account <https://github.com/join?source=header-home>`_ if you haven't done so yet. A GitHub account is completely free and provides you with hosting for your software projects. It also allows you to fork other repositories to create your own improved version of a piece of code.

**Step 2.1: Sign-up for a GitHub account**
    As a student you can upgrade your free account to a personal account (normally $7/month). With a personal account you can have unlimited private repositories. Private repositories can be very helpful when developing new unreleased programs and/or writing your thesis.

**Step 2.2: Apply for a free account upgrade**

   Applying for the free upgrade is easily done via `this <https://education.github.com/discount_requests/new>`_ link and tick boxes: Student and Individual account. Then enter the requested information.

Clone the Tudat Bundle
~~~~~~~~~~~~~~~~~~~~~~
Open SmartGit and click ``Repository`` and ``Clone``. Set the remote url to the Tudat Bundle repository: https://github.com/tudat/tudatBundle.git . Leave the other options default, unless you know what you're doing. 

  .. tip:: It is good practice to avoid long paths and names with special characters (preferably also no spaces). Pick something short and sensible. Especially If you are working on a Windows system, it is recommended that you place the Tudat Bundle directly in your C:\ drive, to avoid problems with filepath lengths.
  
Note that this step can take some time (up to tens of minutes), be patient and don't break off the cloning. 

Congratulations! You have now downloaded the tudatBundle. You can now head to the :ref:`configureTudatLibraries` guide to configure the bundled libraries correctly and build them.



