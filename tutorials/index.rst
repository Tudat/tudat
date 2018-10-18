.. _tutorialsIndex:

Tutorials and Documentation
===========================

Now that you have downloaded, compiled, and tested Tudat, you can start setting up your own applications. At present, we have three interfaces for Tudat:

* C++
* JSON files
* Matlab (unsupported)

The C++ interface allows you access to the full functionality of Tudat, but requires all input to be defined as C++ code. Using the JSON interface, you are able to use majority of numerical state propagation features of Tudat, although the ability to, for instance, define thrust or aerodynamic guidance is limited. Note that the JSON interface is relatively new as of the time of writing (finalized early October 2017), and the majority of documentation and examples relate to the C++ interface. The Matlab interface uses the JSON interface, but is currently not supported by the Tudat development team.

Depending on whether you are planning to use the C++ or JSON interface, we recommend that you start reading the documentation at:

   * :ref:`walkthroughsIndex` (for the C++ interface). Here, all the example applications included in the Tudat Bundle are systematically explained.
   * :ref:`jsonInterface` (for the JSON interface). 

The :ref:`tudatFeaturesIndex` pages discuss the details behind the various interfaces and options that Tudat provides, with a focus on numerical state propagation and estimation. In these pages, you will find a comprehensive list of all your options and limitations when using Tudat.

Several key C++ features: mainly :class:`shared_ptr`, :class:`make_shared`, :class:`bind`, :class:`function` and various Eigen types (vectors, matrices, etc.) are discussed in the :ref:`preknowledgeIndex` pages. These are very helpful to understand the code structure of the Tudat libraries. 

If encountering any problems please visit the :ref:`troubleshootingTutorial` pages. 

.. toctree::
   :hidden:

   applicationWalkthroughs/index
   tudatFeatures/index
   basics/index
   jsonInterface/index
   matlabInterface/index
   Doxygen <http://doxygen.tudat.tudelft.nl>
   troubleshooting/index




