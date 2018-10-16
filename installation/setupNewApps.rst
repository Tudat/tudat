.. _setupNewApps:

   .. warning:: This page is due for removal, and to be replaced with the new :ref:`createNewApps` page.

Setup New Applications
======================

This guide will show how to start adding new applications to the Bundle environment. The guide assumes that you you have already set-up the build environment, and downloaded, configured and build the Tudat libraries.

New applications are typically added to the ``tudatBundle/tudatApplications`` folder, which is empty by default. Moreover, these new applications are expected to have their own github repository, since the tudatBundle git repository is set to not track any changes inside the tudatApplications folder.

You can either start a new (empty) repository or continue developing from an existing repository (for instance by retrieving it from Github). You also have the options to make the repository private or public. This leads to the following three slightly different approaches:

- **A:** Create a new empty repository (private or public).
- **B:** Continue an existing public repository and make it private (for AE4867 students).
- **C:** Continue an existing private/public repository and keep it as is.

Note this same **A**, **B**, **C** notation is used further on in the guide as well. You can skip sections not applicable to your use case.

**Step 1: Setting up a repository for your app**
    **(A + B)** Go to https://github.com/new and create a new repository for your application. It is recommended (and for AE4867 **required**) to set the scope to private. If the private option is unavailable to you, go back and follow the guide on creating a GitHub account. Note that it is not easy to change the name of your repository, so choose wisely!

    **(B)** Create your own private repository from the code. Scroll down to the bottom and choose ``Import code``. Fill in the the url of the repository (in this case https://github.com/tudat/tudatAssignments.git). Click ``Begin import``.

    **(C)** Fork a repository and keep it public. Visit any repository on GitHub and click ``Fork`` at the top of the page.

**Step 2: Clone the repository in SmartGit**
    Open SmartGit, go to ``Repository`` and choose ``Clone...`` Copy in the repository URL and choose ``Next``. If your repository is private, you will need to enter your GitHub credentials! Select the local path ``tudatBundle/tudatApplications/myApplication``. In general, select the name ``myApplication`` according to your preferences. For AE4867, use ``tudatBundle/tudatApplications/numericalAstrodynamics`` (create this directory yourself).

**Step 3: Copy and modify the template (A only)**

   .. note:: AE4867: Only go through this step if installing a new empty repository. Do not go through this step for AE4867!

   Open your file browser and navigate to ``tudatBundle/tudatExampleApplications/templateApplication``. Copy the contents to ``tudatBundle/tudatApplications/myApplication``. Rename the  ``TemplateApplication`` folder to ``MyApplication``. Open ``MyApplication/CMakeLists.txt`` and replace the CMake project name ``TemplateApplication`` with your application:

   .. code-block:: cmake

      # Specify minimum CMake version required.
      cmake_minimum_required(VERSION 2.6)

      # Specify project name.
      project(myApplication)

   Furthermore, one can change the name of your code from ``helloWorld`` to ``myApplication``. Don't forget to change the names at the end of the ``CMakelists.txt`` accordingly:

   .. code-block:: cmake

      # Add helloWorld application.
      add_executable(application_myApplication "${SRCROOT}/myApplication.cpp")
      setup_executable_target(application_myApplication "${SRCROOT}")
      target_link_libraries(application_myApplication tudat_gravitation tudat_basic_astrodynamics ${Boost_LIBRARIES} )

   A detailed explanation of the ``CMakeLists.txt`` file can be found in :ref:`applicationCMakeLists`. 


**Step 4: Add application to the Tudat Bundle**
    Open the top-level ``tudatBundle/CMakeLists.txt`` and add the directory to the bottom of the file using the folder names as just set-up::

        add_subdirectory( "${PROJECTROOT}/tudatApplications/myApplication/MyApplication")

    In case you are following AE4867, add the following instead::

        add_subdirectory( "${PROJECTROOT}/tudatApplications/numericalAstrodynamics/assignment1")
        add_subdirectory( "${PROJECTROOT}/tudatApplications/numericalAstrodynamics/assignment2")

**Step 5: Commit changes to the application repository**

    .. note:: AE4867: This step is only applicable when you have made changes to your code. You can skip this step for now.

    When the application is modified or new code is added, it is recommended to commit these changes to the repository and push these changes to your GitHub account. In SmartGit, select the changes and right-click and choose ``Commit...``. Write a short commit message describing the changes just made. You can now choose to ``Commit`` (just record) or ``Commit & Push`` (record and update GitHub). 

**Step 6: Invite a collaborators to (private) repository**
    Go to the ``Settings/Collaborators`` of the repository. Add other users that you want to invite to the repository.

You have now reached the end of the installation documentation and are ready to set-up your own applications and add those of others. It is a good idea to have a look at all the applications that already came with Tudat and follow the Tutorials and Documentation. See you around!
