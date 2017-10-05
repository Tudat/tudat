.. _matlabInterface_installation:

Installation
============

1. Clone the `tudatBundle:json <https://github.com/aleixpinardell/tudatBundle/tree/json/>`_ repository and its submodules:

  .. code-block:: txt
    
    git clone -b json --recursive https://github.com/aleixpinardell/tudatBundle

  If you don't have git installed on your machine, download it from `here <https://git-scm.com/downloads>`_.

2. Install `CMake <https://cmake.org>`_.

3. Install MATLAB R2016b or later.

4. After choosing a definitive location for your tudatBundle directory, run the script :class:`tudatBundle/matlabInterface/setup.m` in MATLAB.

5. Run the script :class:`tudatBundle/matlabInterface/build.m` in MATLAB to compile all the necessary targets for the MATLAB Interface to work.

6. If you get an error or prefer to compile the targets manually, you can use Qt Creator (see how to :ref:`setupDevelopmentEnvironment`). After the project tudatBundle has been built, run the MATLAB Interface's unit tests by writing in MATLAB's Command Window:

  .. code-block:: matlab
    
    tudat.test()

If you want to move, rename or delete your tudatBundle/matlabInterface directory, remove it first from MATLAB's path by using :literal:`pathtool`.

When building with Qt Creator and running the tests in MATLAB, if all of them fail, the source of the problem may have to do with MATLAB being unable to locate some of the dynamic libraries. Pick one of the tests and run in from Qt creator. It it passes, then in Qt creator go to Projects > Build > Build Environment and check the value of the variable PATH. If some of the paths indicated there refer to Qt, you will have to add them to the MATLAB Interface manually. For instance, in MATLAB's Command Window:

.. code-block:: matlab

  tudat.PATH('C:\Qt\Tools\mingw492_32\bin')
  
This setting is set permanently, so you won't need to run this command again in future MATLAB sessions. If you need to specify several paths, use :literal:`;` as separator (on UNIX systems, using :literal:`:` as separator is also valid).

After performing this step, try running the tests again. If all of them fail, open an issue on GitHub. If only some of them fail, open an individual issue on GitHub for each. Depending on the severity and percentage of features affected by the problem, you may be able to run the examples in the :class:`tudatBundle/matlabInterface/Examples` folder.
