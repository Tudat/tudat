.. _setupNewApps:

Setup New Applications
======================

This guide will show how to start add new and existing applications to the Tudat Bundle environment. Existing applications (by another user) are considered new to the system, thus the process is the same. The guide assumes that you you have already set-up the build environment, and downloaded, configured and build the Tudat libraries.

New applications are added to the ``tudatBundle/tudatApplications`` folder, which is empty by default. Moreover, new applications are expected to have their own repository, since the tudatBundle git repository is set to not track any changes inside the tudatApplications folder by design.

As mentioned previously, you can start a new (empty) repository or continue developing an existing repository. You also have the options to make the repository private or public. This leads to the following three slightly different approaches:

- **A:** Create a new empty repository (private or public).
- **B:** Continue an existing public repository and make it private (for Space Project students).
- **C:** Continue an existing private/public repository and keep it as is.

Note this same **A**, **B**, **C** notation is used further on in the guide as well. You can skip sections not applicable to your use case.

**Step 1: Setting up a repository for your app**
    **(A + B)** Go to https://github.com/new and create a new repository for your application. It is recommended to set the scope to private. If the private option is unavailable to you go back and follow the guide on creating a GitHub account. Note that it is not easy to change the name of your repository, so choose wisely!

    **(B)** Fork a repository and make it private. Scroll down to the bottom and choose ``Import code``. Fill in the the url of the repository (in this case https://github.com/tudat/tudatAssignments.git). Click ``Begin import``.

    **(C)** Fork a repository and keep it public. Visit any repository on GitHub and click ``Fork`` at the top of the page.

**Step 2: Clone the repository in SmartGit**
    Open SmartGit, go to ``Repository`` and choose ``Clone..`` Copy in the repository URL and choose ``Next``. If your repository is private, you will need to enter your GitHub credentials! Select a local path.

    .. tip:: It is highly recommended to create a new directory myApplication (in this guide, substitute for the repository name) under tudatBundle/tudatApplications.

**Step 3: Copy and modify the template (A only)**

    .. note:: Space Project: Only go through this step if installing a new empty repository. do not go through this step for the space project!

    Open your file browser and navigate to ``tudatBundle/tudatExampleApplications/templateApplication``. Copy the contents to ``tudatBundle/tudatApplications/myApplication``. Rename ``TemplateApplication`` to ``MyApplication``. Open ``MyApplication/CMakeLists.txt`` and replace the CMake project name ``TemplateApplication`` with your application.

**Step 4: Add application to the Tudat Bundle**
    Open the top-level ``tudatBundle/CMakeLists.txt`` and add the directory to the bottom of the file using the folder names as just set-up::

        add_subdirectory( "${PROJECTROOT}/tudatApplications/myApplication/MyApplication")

    In case you are following the space project, add the following instead::

        add_subdirectory( "${PROJECTROOT}/tudatApplications/myApplication/spaceProjectAssignment1")

**Step 5: Commit changes to the application repository**

    .. note:: Space Project: This step is only applicable when you have made changes to your code. you can skip this step for now.

    When the application is modified or new code is added, it is recommended to commit these changes to the repository and push these changes to your GitHub account. In the case below the added code that was copied from the template is committed. This process, however, is general for any change to your application. Select the changes and right-click and choose ``Commit..`` Write a short commit message describing the changes just made. You can now choose to ``Commit`` (just record) or ``Commit & Push`` (record and update GitHub).

**Step 6: Invite a collaborators to (private) repository**
    Go to the ``Settings/Collaborators`` of the repository. Add other users that you want to invite to the repository.

You have now reached the end of the installation documentation and are ready to set-up your own applications and add those of others. It is a good idea to have a look at all the applications that already came with Tudat and follow the Tutorials and Documentation. See you around!
