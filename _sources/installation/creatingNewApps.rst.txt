.. _createNewApps:

Creating Tudat Applications
===========================

Now that you've installed Tudat, you can run all of the example applications and unit tests at your leisure. On this page, we will describe how to create other Tudat applications, specifically tuned to your research/education project. There are two broad options for creating applications: pulling an existing application from Github (or elsewhere), or creating an application from scratch. Below, both are discussed. New applications are typically added to the ``tudatBundle/tudatApplications`` folder, which is empty by default.  

   .. note:: It is not the goal of this page to give a full-fledged introduction to repository management, git or Github. Much more details can be found online, for instance in the `Resources to learn Git <https://try.github.io/>`_ or `Github Help <https://help.github.com/>`_. An extensive online book on git is given `here <https://git-scm.com/book/en/v2>`_ This page is meant to serve as a starting point for using git with Tudat and using examples that are very close to a typical user's first Tudat experience.

.. _gettingExistingApp:

Getting an Existing Tudat Application
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Below, we discuss how to retrieve an application from an existing Github repository, how to retrieve any updates that may come available after you've downloaded it, and how to send the code to your own Github account. As an example, we will use the ``tudatAssignments`` repository. Depending on whether your repository is to be public or private, follow Step 1a, or Step 1b.

      .. warning:: IMPORTANT! When setting up a repository for use in the AE4866 or AE4868 courses, it is **required** to set the repository to private (and follow step 1b, not 1a).

**Step 1a (for a public repository): Fork the code on Github**

   The first step in the workflow of retrieving a public application is to "fork" it. Forking means that you will create a copy on your own Github page, which is identical to the original. The reason why we do this is that, typically, you will *not* have permission to modify the original repository, and have to save any modifications you make under your own account. To do so, first make sure you are logged into Github, and then go the repository page for the application you want to retrieve, and click the "fork" button at the top (see screenshot below).

   .. figure:: images/forkingExample.png

   After doing so, you will have created your own local version of the ``tudatAssignments`` repository.

**Step 1b (for a private repository): Import the code on Github**

   The first step in the workflow of retrieving a private application is create an empty new Repository, by going to ``https://github.com/YourName`` (fill in your own name), and clicking ``New`` under ``Repositories`` (see screenshot below):

   .. figure:: images/newRepository.png

   Type a name for your repository (make sure it is descriptive, so not ``NumAstroAssignment``, but, for instance, ``NumericalAstrodynamics2018Assignment``), and set the repository to **private** (see screenshot below). 

   .. figure:: images/newRepositoryPrivate.png

   You may, but need not, write a repository description. Click ``Create repository``. 

   On the bottom of page you are redirected to, click ``Import Code``. Fill in the name of the repository from which you want to retrieve code (for instance ``https://github.com/Tudat/tudatAssignments``), and click ``Begin Import``. Once the import is successfully completed, you will receive an e-mail, with a link to your new, private, repository.

   For a private repository, you can control who can view/modify your repository. To add other users to your repository, go to the ``Settings/Collaborators`` of the repository, and add other users that you want to invite.

      .. warning:: IMPORTANT! When setting up a repository for use in the AE4866 or AE4868 courses, invite :literal:`@dominicdirkx` and :literal:`@transferorbit`. For these courses, it is prohibited to invite your fellow students to a private repository for an individual assignment.

**Step 2: Clone the application to your system**

   Now, the next step is to create a version of the code on your computer, which is linked to your own Github page. Using the terminal (or the ``tudat_shell.bat`` program when using Windows), navigate to the ``tudatBundle/tudatApplications`` directory. Then, use the following command in the terminal::

      git clone https://github.com/YourName/tudatAssignments.git

   where ``YourName`` should be replaced with your Github account name (so that it corresponds to the URL to where you've forked the repository).


**Step 3: Adding the application to your Tudat compilation**

   Open the top-level ``tudatBundle/CMakeLists.txt`` and add the directory to the bottom of the file using the folder where you've just cloned the repository to::

      add_subdirectory( "${PROJECTROOT}/tudatApplications/tudatAssignments")

   When opening the tudatBundle project in Qt, you'll see your new application(s) added to the project tree, and the corresponding executable added to the list of applications you can run.

Pulling and Pushing Application Code
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Now, we are ready to discuss two distinct but related aspects of developing your code: retrieving modifications from the original repository, called pulling (here: ``https://github.com/tudat/tudatAssignments.git``) and uploading your modifications to your own repository, called pushing (here: ``https://github.com/YourName/tudatAssignments.git``).

**Step 1: Setting up your remotes**

   In git, a ``remote`` is an address of an external repository (in this case on Github). You can set any number of remotes you like for a given repository. You can view all your remotes for a given Git repository by using the ::

      git remote -v

   command in your terminal. Running this command will likely result in the output (for now, forget about the distinction between pull and fetch)::

      origin	https://github.com/YourName/tudatAssignments (fetch)
      origin	https://github.com/YourName/tudatAssignments (push)
 
   Typically, you will have two: an ``origin`` and an ``upstream``, which is also the convention we'll stick to here. The ``origin`` remote is the one from where you've cloned the repository, in this case your own Github version of the ``tudatAssignments`` repository. This remote will have been set automatically when cloning the code. With the way the repository is set up, you are ready to update your own Github version of ``tudatAssignments``. To also retrieve new code from the original tudat repository, we need to add an additional remote, the ``upstream``. To do so, use the following terminal command::

      git remote add upstream https://github.com/tudat/tudatAssignments.git

   Rerunning the ``git remote -v`` command should now result in::

      origin	https://github.com/YourName/tudatAssignments (fetch)
      origin	https://github.com/YourName/tudatAssignments (push)
      upstream	https://github.com/tudat/tudatAssignments (fetch)
      upstream	https://github.com/tudat/tudatAssignments (push)

**Step 2: Making local commits**

   Before you can push changes to your Github account, you must first ``commit`` (save) these changes locally on your computer. A commit provides a snapshot of the current version of the code, to which you can return at later points in time.

   .. tip:: When finishing a part of the code to your own satisfaction, or making clear progress in your work, commit your code. This does not override any old or later commits, but will provide a way to go back to your current version of the code.

   The first step in making commits is usually to check what has changed w.r.t. the previous commit. To check this, type the command::

      git status
   
   This should given an output similar to that given below:

   .. figure:: images/gitStatusExample.png

   In this example, we have modified two existing files, and created a new file. 

   Now, before committing, you must ``stage`` changes for commit. To stage all changes shown by the ``git status`` command for commit, use::

      git add . 

   For the example given above, this will result in (after running ``git status`` again):

   .. figure:: images/gitAddExample.png
   
   If you only want to stage a single file, or folder, use::

      git add FolderName/
      git add FileName.ext
   
   where ``.ext`` is just an arbitrarily chosen extension. You can also use::

      git add FolderName/FileName.ext

   to stage a single file in a folder. You can combine as many ``git add`` commands as you like to stage all your files for commit.

   Now, committing your code is done by::

      git commit -m "Your commit text here"

   The text between the quotes will show up in your commit log, and should ideally describe the current state of your code: which changes have you made since the last commit?
  
**Step 3: Pushing your commits to Github**

   After committing the code, you will have made a snapshot of the current version of the code, on your local system only. If you want to share it with others, the best way is to push it to an online repository (typically Github). Assuming you've set up your remotes as defined above, you can use::

      git push origin master

   This will push your code to the ``origin`` remote. The ``master`` term denotes the current branch you are working on. Branch management is beyond the scope of this tutorial, and you may safely ignore these issues for now.

   If you have set up your repository to be private, you will be prompted to enter your Github username and password. After a push is succesfull, you should see your changes on the Github page for your applications, available for everyone (in case it is public) or a selected few (in case it is private).

Updating Your Local Repository from Github
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

After following the above guide, you'll have created a git repository on your computer, based on some remote from Github. Typically, commits are regularly done to Github repositories, and it may happen that you want to update your local code with the modifications of the remote. For this example, we'll assume that there has been some modification of the ``upstream`` remote (i.e., at ``https://github.com/tudat/tudatAssignments``), and that you want to update your local version of the code (both on your computer, and on your own Github page). 

**Step 1: Committing any local changes**

   Before pulling the latest code from Github, commit any modifcations you have made, using the ``git add`` and ``git commit`` commands described above. Not doing so will cause the following error message when pulling::

      error: Your local changes to the following files would be overwritten by merge:
         file_name
      Please, commit your changes or stash them before you can merge.
      Aborting

   This error is given as a safety measure, since pulling the latest version of the code may inadvertently, and irreversibly, overwrite your own changes. 

   .. note:: In case it is your intention to overwrite the changes you have made locally, you can use the command ``git reset --hard``. Note however, that this step is **irreversible**!

**Step 2: Fetching and Pulling the Remote**

   The next step in updating the code is to type::

      git fetch upstream

   The ``fetch`` command does not update the code on your computer, but makes your local git repository aware of any changes make to a remote (the ``upstream`` in the above example). Next, you will ``pull`` the code from the ``upstream`` with the following command:

      git pull upstream master

   Note that we are still assuming that only the ``master`` branch is relevant for our current application. The ``pull`` command will have one of two possible outputs (assuming you have correctly performed step 1). Either no error is given, and the pull has been succesful, or there are conflicts with changes you have made, which will give the following error message::

      Pull is not possible because you have unmerged files.
      Please, fix them up in the work tree, and then use 'git add/rm <file>'
      as appropriate to mark resolution, or use 'git commit -a'.
 
   In case you get this message, go to step 3.

**Step 3: Solving Conflicts (if needed)**

   As is often the case, changes you have committed on your own computer will not be compatible with changes that have been made to the remote you are pulling. The list of files with merge conflicts will be shown when using the ``git status`` command. The resulting merge conflicts are typically corrected manually, where the user decides how to update the code after a pull. A merge conflict will show up in your code as::

      <<<<<<< HEAD
      Remote modifications
      =======
      Your modifications
      >>>>>>> master

   Clearly, this code will not compile anymore. You can change this block to either::

      Remote modifications

   or::

      Your modifications

   or something else entirely, as you see fit for the case at hand. 

   After correcting all conflicts, use the ``git add`` and ``git commit`` commands to commit your merged code.

Creating a New Tudat Application
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For some projects, you will want to start your own application repository from scratch. Here, we briefly explain how to set this up, while details of the code itself (e.g., CMake settings) are discussed in the following sections.

**Step 1: Initializing the repository**

   To create a new git repository, use the terminal to navigate to the directory of this new repository and type::

      git init

   This will create a new, empty, repository in the current directory. Using the same ``git add`` and ``git commit`` commands as above, you can add files to the repository as you see fit. 

   Before (or after) doing so, you can add a ``.gitignore`` file to your repository (see Tudat repository for a typical example). This file can contain a list of files, directories, file extensions, etc., that git will normally *ignore* when using the ``git status`` or ``git add`` commands. For example, you may want to keep ``.dat`` files, or a ``bin/`` directory outside of your repository. As an example, the ``tudatApplications`` directory is in the ``.gitignore`` list of ``tudatBundle``, as application commits are not added to the bundle repository.

**Step 2: Creating a Github Repository**

   Now that you've created a local repository on your system, you need to create a new Github project, to which you can push your code. On the `Github main page <https://github.com/>`_, click ``Start a Project`` (make sure you are logged in first). You will be prompted to provide some baisc information on your new repository (and declare it public or private). After clicking ``Create Repository``, your new (empty) Github repository will be created.

   Now, we need to tell your local repository where this new Github project is located. Using the same tools as above, use the ``git remote add`` command to add your new repository as the ``origin``. For instance::

      git remote add origin https://github.com/UserName/MyNewTudatApplication.git

   You are now free to push your code to this repository.

.. _writingCMakeLists:

Writing Your CMakeLists, and Starting Your Code
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The above guides show you how to push, pull, commit, etc., using the git version control system. In this last part of the guide on how to set up new applications, we will discuss the basic aspects that the ``CMakeLists.txt`` for your application, and your C++ code, must adhere to. Note that this part of the guide is primarilly relevant if you want to create your own application code from scratch. However, it will also provide insight into how/why to modify the ``CMakeLists.txt`` file for a project you've pulled from Github.

To make your life easier, we have created a ``TemplateApplication`` directory in the example applications. You can copy and past the ``CMakeLists.txt`` file (see `file on Github <https://github.com/Tudat/tudatExampleApplications/blob/master/templateApplication/TemplateApplication/CMakeLists.txt>`_) in this directory to your application. To adapt it to your specific needs, you will typically only need to make minimal modifications. Add the bottom of the ``CMakeLists.txt`` file, you'll see::

   # Add helloWorld application.
   add_executable(application_HelloWorld "${SRCROOT}/helloWorld.cpp")
   setup_executable_target(application_HelloWorld "${SRCROOT}")
   target_link_libraries(application_HelloWorld tudat_gravitation tudat_basic_astrodynamics ${Boost_LIBRARIES} )

These lines of CMake code add an application to your project, which can then be compiled and run (note that lines starting with ``#`` are treated as comments in CMake). In this case, it is the ``helloWorld.cpp`` file, located in the ``${SRCROOT}`` directory (which denotes the directory of the ``CMakeLists.txt`` file).

.. note:: Any C++ application must contain one, and only one, ``int main`` function. The file containing this function may contain any number of additional function definitions.

Often, it is only in these lines where you will modify the ``CMakeLists.txt`` file. As an example, say you want to add the code in the ``myNewTudatApplication.cpp`` file (located in the same directory as your ``CMakeLists.txt``), and compile into an executable named ``application_MyNewApplication``::

   # Add myNewApplication application.
   add_executable(application_MyNewApplication "${SRCROOT}/myNewTudatApplication.cpp")
   setup_executable_target(application_MyNewApplication "${SRCROOT}")
   target_link_libraries(application_MyNewApplication ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES} )

assuming you have already added your application to your Tudat bundle CMakeLists file (see Step 4 of :ref:`gettingExistingApp`). Depending on the details of your application, the final line ``target_link_libraries`` may look slightly different. There are many options to change it, but for most Tudat applications it will be sufficient to use the above line, or::

   target_link_libraries(application_MyNewApplication ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES} )

which is required if any of the state-estimation-related functionalities are needed (variational equations propagation, observtion models, acceleration partials, orbit determination, etc.). If you don't need this functionality, using ``${TUDAT_PROPAGATION_LIBRARIES}`` can save some compilation time.

.. tip:: When you receive a compilation error with the words ``undefined reference to ...``, or ``Undefined symbols for architecture x86_64``, this is typically indicative of the ``target_link_libraries`` being set incorrectly.

Adding More Files to Your Application
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. note:: This section assumes that you have already gained some basic knowledge and experience of (Tudat) code development.

At some point in the development of your application, you may reach a point where you don't want to cram all your application code into a single file (like the ``myNewTudatApplication.cpp`` file in the application above), and you'll want to spread your code over multiple files, as we do in Tudat. Here, we'll show an example in which we want to add the code in the following files::

   newFile1.cpp
   newFile1.h
   NewFolder/newFile2.cpp
   NewFolder/newFile2.h
   NewFolder/newFile2.cpp
   NewFolder/newFile2.h

Unfortunately, it is not sufficient to add the correct ``#include`` statements in your code. Doing so, and not modifying the CMakeLists file, will result in the undefined reference error mentioned above. First, the following code must be added, just before your ``add_executable`` commands::

   # Set the source files.
   set(MY_NEW_APPLICATION_SOURCES
     "${SRCROOT}$/newFile1.cpp"
     "${SRCROOT}$/NewFolder/newFile2.cpp"
     "${SRCROOT}$/NewFolder/newFile3.cpp"
   )
   
   # Set the header files.
   set(MY_NEW_APPLICATION_HEADERS
     "${SRCROOT}$/newFile1.h"
     "${SRCROOT}$/NewFolder/newFile2.h"
     "${SRCROOT}$/NewFolder/newFile3.h"
   )

   # Add static libraries.
   add_library(my_new_application_libration STATIC ${MY_NEW_APPLICATION_SOURCES} ${MY_NEW_APPLICATION_HEADERS})
   setup_library_target(my_new_application_libration "${SRCROOT}{AERODYNAMICSDIR}")

This code will add the library ``my_new_application_libration`` to your project. This library will contain the compiled code of the ``MY_NEW_APPLICATION_SOURCES`` source and ``MY_NEW_APPLICATION_HEADERS`` header files. Now, the final step, to allow these six new files to be used in your application, is to update the ``target_link_libraries`` to::

   target_link_libraries(application_MyNewApplication my_new_application_libration ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES} )

.. note:: The names chosen here for the new source/header files, library, etc., are for illustrative purposes only. Feel free to modify them as you see fit.

As a final point, be aware that you may add any number of applications to your CMakeLists file (each with exactly one ``int main`` function. As an example, below is a selection of the ``SatellitePropagatorExamples`` CMakeLists.txt (slightly edited for readibility)::

   # Add Galileo constellation application.
   add_executable(application_GalileoConstellationSimulator "${SRCROOT}/galileoConstellationSimulator.cpp")
   setup_executable_target(application_GalileoConstellationSimulator "${SRCROOT}")
   target_link_libraries(application_GalileoConstellationSimulator ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES} )

   # Add JSON-based Apollo propagation    
   add_executable(application_ApolloEntryJSON "${SRCROOT}/apolloCapsuleEntryJSON.cpp")
   setup_executable_target(application_ApolloEntryJSON "${SRCROOT}")
   target_link_libraries(application_ApolloEntryJSON json_interface_library ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES} )

   # Add simulated Earth orbiter simulated POD example
   add_executable(application_EarthOrbiterStateEstimation "${SRCROOT}/earthOrbiterStateEstimation.cpp")
   setup_executable_target(application_EarthOrbiterStateEstimation "${SRCROOT}")
   target_link_libraries(application_EarthOrbiterStateEstimation ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES} )

In the same folder, three separate applications have been added.
