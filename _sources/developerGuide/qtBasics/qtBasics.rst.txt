.. _qtBasics:

Qt Basics
=============

Qt Creator is a cross-platform integrated development environment (IDE), and is suggested to be used by Tudat users new to coding and compiling C++ code. This page describes some of the basics of Qt Creator. The layout of the IDE is shown below:

.. figure:: images/QtCreatorLayout.png


Project Tree
~~~~~~~~~~~~
The project tree is indicated by the red area. It contains the top-level ``CMakeLists.txt`` which contains the settings for the Tudat project, this file is discussed in :ref:`CMakeListTutorial`. It contains the external libraries used by Tudat: boost, cspice, eigen, nrlmsise-00 and sofa. The tudat folder contains all the code of the Tudat project as documented in :ref:`tudatFeaturesIndex`. Furthermore, the example applications as discussed in :ref:`walkthroughsIndex` can be found in the last folder in the project tree. Your application can be added to this project tree list. How this is done can be found :ref:`setupNewApps`.

Output console
~~~~~~~~~~~~~~
The output console (indicated with green) is used to displays any issues encountered during compilation which help pinpoint any compilation problems. Furthermore, it is used to display the output of search results, applications and compilation of application. The last tab `General Messages` shows the output of CMake. 
 
Warnings can be muted. Go to the issues tab of the output console and click icon with the yellow triangle with exclemation mark. Now only the errors are shown. If at any time you also want to see the warnings, just click the icon again.
