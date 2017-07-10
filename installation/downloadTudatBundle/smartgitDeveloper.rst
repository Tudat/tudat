.. _downloadTudatBundleSmartgitDeveloper:

Download using SmartGit (Developer)
-----------------------------------

.. tip:: Before following this guide, you are highly advised to go through the :ref:`githubBasics` section to get up to speed with the usage of GitHub.

1. Setup SmartGit
~~~~~~~~~~~~~~~~~
**Step 1.1: Download SmartGit**
    Download SmartGit for your operating system. Note that for Linux the package might be available through your package manager, see Linux/Debian instructions. Note that SmartGit comes bundled with JRE (java runtime environment) on Windows and Mac OS X. On Linux you need to install ``jre-openjdk``.

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


2. Create a GitHub account
~~~~~~~~~~~~~~~~~~~~~~~~~~
Create a GitHub account if you haven't done so yet. A GitHub account is completely free and provides you with hosting for your software projects. It also allows you to fork other repositories to create your own improved version of a piece of code.

**Step 2.1: Sign-up for a GitHub account**
    As a student you can upgrade your free account to a personal account (normally $7/month). With a personal account you can have unlimited private repositories. Private repositories can be very helpful when developing new unreleased programs and/or writing your thesis.

**Step 2.2: Apply for a free account upgrade**

3. Fork both the Tudat Bundle and the Tudat repositories and Clone them locally
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**Step 3.1: Fork Tudat Bundle**
    Create your own copy of Tudat and Tudat Bundle under your GitHub account. Go to https://github.com/Tudat/tudatBundle and click on ``Fork``.

**Step 3.2: Fork Tudat**
    Go to https://github.com/Tudat/tudat and click on ``Fork``.

**Step 3.3: Clone your Tudat Bundle**
    Open SmartGit and click ``Repository`` and ``Clone``. Set the remote url to your Tudat Bundle repository, like so https://github.com/username/tudatBundle.git . Leave these options default, unless you know what you're doing. 

    .. tip:: It is good practice to avoid long paths and names with special characters (preferably also no spaces). Pick something short and sensible. Especially If you are working on a Windows system, it is recommended that you place the Tudat Bundle directly in your C:\ drive, to avoid problems with filepath lengths.


4. Synchronize with official Tudat Bundle
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**Step 4A: Set-up official Tudat Bundle (Do this once!)**
   Choose ``Remote`` and ``Add``. Set the url as the official Tudat Bundle repository url: https://github.com/Tudat/tudatBundle.git and call it "upstream". In your favorite text editor open ``tudatBundle/.gitmodules`` and change the url, like so::

    [submodule "tudat"]
        path = tudat
        url = https://github.com/[YOURUSERNAME]/tudat.git

   .. warning:: On Mac OS X this file is hidden, you can right-click tudatBundle in the repository tree in SmartGit and select ``Open`` in Terminal. Type `open .gitmodules` followed by enter. Lastly, you need to sync your submodules, choose ``Remote``, ``Submodule`` and ``Synchronize``.

**Step 4B: Update to the latest Tudat Bundle (Do this every time!)**
    Right-click on ``upstream`` and select ``Pull``. If a dialog pops up asking you to Rebase or Merge, select ``Rebase`` and click ``Configure``. Choose ``Fetch only``. Right-click again on the ``upstream`` entry and select ``Merge`` this time. You will be asked how to incorporate any changes. In almost all cases ``Fast forward`` is the best option. In case there are conflicts you should do a merge. You have now succesfully updated your local Tudat Bundle to the latest version. However, you need to synchronize this change to your remote version as well.


5. Synchronize with official Tudat
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**Step 5A: Set-up official Tudat (Do this once!)**
    In the ``Repository`` window select Tudat and right-click the ``origin`` entry in the ``Branches`` window and select ``Properties``. Change the official remote Tudat repository url to your forked repository. Now, all that remains is to re-add the official remote Tudat repository url like we did with the Tudat Bundle. First, make sure that Tudat is still selected in the ``Repositories`` window. Then, click ``Remote Add`` from the top menu. Like before fill in the official address for the remote url: https://github.com/Tudat/tudat.git. Again choose "upstream" as the name.

**Step 5B: Update to the latest Tudat (Do this every time!)**
    Right-click on upstream and select ``Pull``. Choose ``Fetch only``. Right-click again on the upstream entry and select ``Merge`` this time. You will be asked how to incorporate any changes. In almost all cases ``Fast forward`` is the best option. In case there are conflicts you should do a merge. You have now succesfully updated your local Tudat to the latest version. However, you need to synchronize this change to your remote version as well. It could be that you have to checkout the local ``master`` branch first if an error message pops up. Double click the local branches, ``master`` branch. Try ``Sync`` again.

Congratulations! You have now downloaded the tudatBundel. You can now head to the :ref:`configureTudatLibraries` guide to configure the bundled libraries correctly and build them.
