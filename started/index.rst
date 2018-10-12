.. TUDAT_Documentation documentation master file, created by
   sphinx-quickstart on Mon Jul  3 06:08:42 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Getting Started
===============
This page is aimed at everyone who is using Tudat for the first time, or for users who want to refresh their knowledge on Tudat. The getting started guide will show new users where to start when they want to use Tudat, but it can also help users who have specific problems by pointing them to a specific part of the wiki where there might be more information on their problem. 

First Steps in Tudat
~~~~~~~~~~~~~~~~~~~~
For everyone who has already chosen Tudat for their next project: welcome! If more information is needed on the features and/or the possible applications of Tudat, please see the `introduction page <http://tudat.tudelft.nl/>`_, which lists several research projects that have made use of Tudat, and furthermore, it gives a small list of features that are currently present in Tudat.

Installation
************
The first step to get started with Tudat is to install all the software that is necessary to write, compile, and run your Tudat code. For a list of these programs and a clear guide on how to install them, users are referred to the `installation guide <http://tudat.tudelft.nl/installation/setupDevEnvironment.html>`_. Which specific programs need to be installed depends on your operating system, windows users can go `to this page <http://tudat.tudelft.nl/installation/setupDevWindows.html>`_, OSx users `here <http://tudat.tudelft.nl/installation/setupDevMacOs.html>`_, and users working on a Unix machine can go to the `Linux page <http://tudat.tudelft.nl/installation/setupDevLinux.html>`_. 

Once you have installed all the software corresponding to your operating system, the next step is to download the Tudat bundle. This bundle contains all the necessary external libraries and code that allows a user to work with Tudat. There are several options for downloading the Tudat bundle, all of them can be found `here <http://tudat.tudelft.nl/installation/downloadTudatBundle/index.html>`_. Which one you should use depends on what you plan to do with Tudat, and your skill with git. Any basic user of Tudat is referred to the `user download guide <http://tudat.tudelft.nl/installation/downloadTudatBundle/smartgitUser.html>`_. 

Once this is done, it is time to configure the Tudat libraries. This part of the installation process usually takes the most time as it encompasses the building of the code and running the unit tests. The steps that need to be taken are the same for every user and can be found on the `configure Tudat libraries page <http://tudat.tudelft.nl/installation/configureTudatLibraries.html>`_. Please do not forget to run the unit tests as these will show you if the installation was performed correctly. 

If the unit tests output shows that all the test have passed, you are ready to get started with Tudat! If there are some unit tests that have failed, or something else went wrong during the installation, users can go to the `troubleshooting page <http://tudat.tudelft.nl/installation/troubleshooting/index.html>`_ to get solutions to problems that have been encountered before. If the error is not on the troubleshooting page, users can go to the issues page on `github <https://github.com/Tudat/tudat/issues>`_. Please look through the previous issues on github first for the error before something new is posted, as there is a large chance that someone else has already encountered the error. If the error has not been encountered before, a new issue can be started and one of the developers will help as soon as possible. 

Once the installation process is done, it is time to start learning how to actually use Tudat and set up a new application. 

Learning Tudat
**************
Tudat is written in the C++ language. Thus, it is important that the user understands how this language works and how they can write their own code in C++. The :ref:`preknowledgeIndex` page has explanations of concepts used in Tudat and links to external tutorials on C++ that a user can follow to get to know this language. This page also contains explanations of several external libraries that Tudat uses. These external libraries are not written by the Tudat developers, but they are needed within Tudat to be able to run most of the applications. Before Tudat can be used, some basic concepts of these libraries need to be understood by the user, which are all explained there.

Once these basic concepts are all understood, it is time to get into the Tudat code. The best way to start learning Tudat is to go through the examples on the :ref:`walkthroughsIndex` page. These example applications are used to teach the user all the concepts needed to be able to write their own application. It contains examples on how to set up an environment, how to set up all accelerations, how to change settings for bodies and vehicles, how to integrate and propagate the state of the vehicle, and much more. For new Tudat users, it is recommended that they go through all of these examples before they proceed to the next step.

The last step before a user can start with Tudat is to set up a new application in which the user can build their own simulation. How this is done is explained here: :ref:`setupNewApps`. 

Now that a new application has been set up, it is time to start building your own Tudat application. As for most user applications the examples given before are not complex enough, users are referred to the rest of the wiki on this website to get a more in depth look at all the options that are included in Tudat. The next section on this page will give a quick overview of how the wiki is structured and where a user needs to go if they have a specific problem.

The Wiki
~~~~~~~~
As Tudat is a large and complex software distribution, it needs a way of explaining all of its different applications and options. This is where the wiki comes in. The rest of this website contains most of the information needed to understand Tudat, and tutorials on how to get started with Tudat (see `First Steps in Tudat`_). This section will give a quick overview of the wiki and give links to the pages that contain certain information.

Most of the information needed can be found on the :ref:`tutorialsIndex` page. The first two sections have already been discussed in `First Steps in Tudat`_ and are only needed if a user is new to Tudat. The page containing most of the information needed on setting up a simulation is the :ref:`tudatFeaturesIndex` page. This page discusses all the libraries in Tudat that can be used to set up a simulation, and it discusses the various options there are when using these libraries. The page contains the following subpages:

- The :ref:`tudatFeaturesAstroIndex` page. This page contains basic information on various astrodynamics tools that can be used in Tudat, e.g. orbital elements, frame translations, frame rotations, etc.
- The :ref:`tudatFeaturesMathIndex` page. All the basic mathematical tools that can be used in Tudat are listed here. It contains information on: interpolators, integrators, and probability distributions.
- The :ref:`tudatFeaturesEnvironmentIndex` page. One of the most important aspects of Tudat is the environment set-up. This page goes over how to set up different environments in the solar system (e.g. gravity fields, atmospheres, rotational models, etc.) and which options the user has for these environments.
- The :ref:`tudatFeaturesAccelerationIndex` page. Another important concept of Tudat is the acceleration set-up. Acceleration(s) on the vehicle/body that is propagated in the simulation need to be defined inside the simulation. This concept is explained on this page, together with all the options available to the user.
- The :ref:`tudatFeaturesSimulatorIndex` page. This is one of the most important aspects of Tudat, as on this page the simulator framework is discussed. All the options the user has, plus the necessary steps the user needs to take to set-up the simulator are discussed on this page.
- The `Estimation Set-Up <http://tudat.tudelft.nl/tutorials/tudatFeatures/estimationSetup/index.html>`_ page. This page discusses the options the user has to simulate observations and estimations of the state of the vehicle, and other variables, during the simulation.
- The :ref:`tudatFeaturesOtherIndex` page. This page lists all other libraries used in Tudat. It contains information on how to input/output certain data, information on SPICE (the library that contains information on certain solar system bodies), and information on JSONCPP. 

After the :ref:`tudatFeaturesIndex` page, there are also pages on the :ref:`matlabInterface` and the :ref:`jsonInterface` for users interested in using these options for their applications. 

Not all the different classes and methods are explained in this wiki. If users can't find a certain feature, or they come across something in Tudat that isn't explained in the wiki, they can go to the `Doxygen <http://doxygen.tudat.tudelft.nl/>`_ page. This page contains all of the classes and namespaces contained in Tudat, including which input parameters are needed and what the output is of certain methods. 

For developers and users who need some more information on advanced subjects, there is the :ref:`devGuideIndex`. In here, there is some more information on: :ref:`githubBasics`, :ref:`updatingTudat`, :ref:`qtBasics`, :ref:`applicationCMakeLists`, :ref:`extendingJSON`, :ref:`extendingMATLAB`, and :ref:`howToWriteTheWiki`.

A request to put something on the wiki can always be made by posting an issue on `github <https://github.com/Tudat/tudat/issues>`_.







