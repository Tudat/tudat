.. _tutorialsIndex:

Tutorials and Documentation
===========================

Now that you have downloaded, compiled, and tested Tudat, you can start setting up your own applications. At present, we have three interfaces for Tudat:

* C++
* JSON files
* Matlab

The C++ interface allows you access to the full functionality of Tudat, but requires all input to be defined as C++ code. Using the JSON/Matlab interface, you are able to use majority of numerical state propagation features of Tudat, although the ability to define thrust or aerodynamic guidance is limited. Note that the JSON/Matlab interface are new as of the time of writing (finalized early October 2017), and the majority of documentation and examples relate to the C++ interface.

To use the C++ interface we recommend to start going through the :ref:`preknowledgeIndex` pages, since these contain essential information on C++. Once you feel confortable, we advise you to look at the :ref:`walkthroughsIndex`, which systematically explain some of the example applications included in the Tudat Bundle. These are very helpful to understand the interaction between the Tudat libraries and how a simulation is properly set up. Next, the :ref:`tudatFeaturesIndex` pages discuss the details behind a number of Tudat libraries common in most simulations.

The final two pages, :ref:`jsonInterface` pages, and :ref:`matlabInterface` pages, discuss in detail the manner in which you can perform numerical orbit propagation _without_ directly being exposed to the C++ interface.

.. toctree::
   :hidden:

   gettingStarted/index
   applicationWalkthroughs/index
   tudatFeatures/index
   jsonInterface/index
   matlabInterface/index
   troubleShooting/index

