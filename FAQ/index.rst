.. _faqIndex:

FAQ
===============

This page will give some answers to frequently asked questions about Tudat. If your question is not anwsered on this page, please open an issue on the Tudat github `page <https://github.com/Tudat/tudat/issues/new>`_.

	- *Q: I would like to use Tudat, where should I start?*
	   A: the :ref:`gettingStarted` page contains a detailed explanation for new users of how to start using Tudat.

	- *Q: I have problems with the installation of the Tudat bundle, how can I solve these problems?* 
	   A: the :ref:`troubleshootingInstallation` page contains several solutions to common errors that can occur during installation. If this does not help, a list of issues raised by users can be found `here <https://github.com/Tudat/tudat/issues>`_. If this still does not contain a solution, a new issue can be made on `github <https://github.com/Tudat/tudat/issues/new>`_ so that the developers can help you.

	- *Q: I would like to use Tudat, but I don't know how to code in C++. Can I still use it?*
	   A: To use Tudat, a little knowledge of C++ is necessary. You can learn :literal:`C++` from scratch using `this <http://www.cplusplus.com/doc/tutorial/>`_ page. If you have already used Matlab before, `this <http://runge.math.smu.edu/Courses/Math6370_Spring13/Lec2.pdf>`_ page lists the differences between the two languages. The same page exists for Python users, and can be found `here <https://pdfs.semanticscholar.org/9ad1/030685050e949d1a3d6d92bababcbe075e07.pdf>`_.
	
	- *Q: I am on Windows and I would like to use the command line for git, how can I do this?*
	   A: A command line tool was developed for Tudat, you can find it in your tudat bundle under external/tools/tudat_shell. Within this tool you can use all the git commands.

	- *Q: Some of the unit test on Windows fail due to high RAM usage, how can I solve this?*
	   A: Turn the COMPILE_HIGH_ACCURACY_ESTIMATION_TESTS off in the Projects tab.

	- *Q: I want to add some features to Tudat, can I do this, and if yes, how?*
	   A: You can always add features to Tudat. If you think that it can also be used by others, you can make a pull request to the tudat repository. The :ref:`devGuideIndex` page contains guides on how to develop extra features for Tudat.

	- *Q: How can I find the definition of a specific class/function/variable in Tudat?*
	   A: In Qt, you can right click on what you would like to inspect and select the "Follow Symbol Under Cursor" option. 

	- *Q: How can I find the usages of a specific class/function/variable in the whole Tudat bundle?*
	   A: In Qt, you can right click on what you would like to find the usages of and select the "Find Usages" option. 

	- *Q: How can I run CMake again from Qt?*
	   A: In the Project Tree tab in Qt (on the left side of your screen) you can right-click on the top-level TudatBundle icon and select the "Run CMake" option. The output will be shown in the General Messages tab.

	- *Q: When running a Tudat propagation application, my application stops running but I don't know why. How can I find out?*
	   A: The dynamics Simulator class has a function: :literal:`getPropagationTerminationReason`, which allows you to retrieve a pointer to the :literal:`PropagationTerminationDetails`. This object has a function called :literal:`getPropagationTerminationReason`, which returns the reason for termination of the propagation.

	- *Q: I have my own custom aerodynamic coefficient file/function that does not work with any of the given options in Tudat, how can I implement this?*
	   A: Tudat has a :literal:`customAerodynamicCoefficientInterface` which can be used for this purpose (see the class definition for more information on how to use this). You do need to make a :literal:`std::function` for your coefficients to be able to use it.

	- *Q: When making a Pagmo problem with a propagation step in it, the optimization stops after a while due to high RAM usage, how can I fix this?*
	   A: Put some variables (especially the creation of the :literal:`bodyMap`) in the constructor of the problem. This will make sure that some vectors will not grow unnecesarily large after several generations.

	- *Q: I get a CSpice error telling me it cannot load any more kernels when running a optimization application, how can I solve this?*
	   A: Make sure that you are not loading the Spice kernels in the fitness function, but either in the constructor or somewhere else that will not be called upon by each time the fitness function is called.

	- *Q: How can I repress the output of the :literal:`dependentVariableSettings`?*
	   A: Set the second argument of the class to false/0.

	- *Q: When making a Pagmo problem with a propagation step in it, the optimization stops after a while due to high RAM usage, how can I fix this?*
	   A: Put some variables (especially the creation of the :literal:`bodyMap`) in the constructor of the problem. This will make sure that some vectors will not grow unnecesarily large after several generations.

	

	
