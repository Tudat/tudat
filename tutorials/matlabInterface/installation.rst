.. _matlabInterface_installation:

Installation
============

1. Clone the `tudatBundle:json <https://github.com/aleixpinardell/tudatBundle/tree/json/>`_ repository and its submodules:

  .. code-block:: txt
    
    git clone -b json --recursive https://github.com/aleixpinardell/tudatBundle

2. Install `CMake <https://cmake.org>`_.

3. After choosing a definitive location for your tudatBundle directory, run the script :class:`matlabInterface/setup.m` in MATLAB.

4. Run the script :class:`matlabInterface/build.m` in MATLAB to compile all the necessary targets for the MATLAB Interface to work.

5. If you get an error or prefer to compile the targets manually, you can use Qt Creator. After the project tudatBundle has been built, run the MATLAB Interface's unit tests by writing in MATLAB's Command Window:

  .. code-block:: matlab
    
    tudat.test()

If you want to move, rename or delete your tudatBundle/matlabInterface directory, remove it first from MATLAB's path by using :literal:`pathtool`.
