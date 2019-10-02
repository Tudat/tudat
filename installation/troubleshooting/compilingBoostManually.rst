.. _troubleShootingInstallationCompilingBoostManually:

Compiling boost manually
========================
   
The top-level CMakeLists.txt of Tudat Bundle, downloads, extracts, configures, builds and installs Boost for you. Although this process is completely automated it can happen that it fails somewhere doing the former.

.. note:: If the automated process fails, it is necessary to take note of where it fails please copy the output of CMake for specifics.

In case it fails perform the following steps:

1. Find the Boost version TudatBundle is trying to build.
    1. Open tudatBundle/CMakeLists.txt
    2. Look for the uncommented (without a # in front) instance of set(BoostVersion 1.XX.0).
2. https://sourceforge.net/projects/boost/files/boost/
    1. Pick the version corresponding to your version. Do not select beta.
    2. It doesn't matter which archive type you select, generally pick .tar.bz2 for Linux and macOS and .zip for Windows.
3. Unpack the folder somewhere, for instance /home/user/boost or c:\boost.
4. Open terminal emulator and go to the Boost folder.
5. Run bootstrap:
    1. Linux and macOS: ./bootstrap.sh --with-toolset=gcc
    2. Windows: .\bootstrap.bat gcc
6. If successful, run bjam2:
    1. Linux and macOS: ./b2 toolset=gcc link=static threading=multi --build-dir=Build stage variant=release --layout=tagged cxxflags=-std=c++11 --with-filesystem --with-system --with-thread --with-regex --with-date_time --with-test
    2. Windows: .\b2.exe toolset=gcc link=static threading=multi --build-dir=Build stage variant=release --layout=tagged cxxflags=-std=c++11 --with-filesystem --with-system --with-thread --with-regex --with-date_time --with-test
7. In case of errors try, to identify if bjam fails for each module or only for select modules.
8. Re-run the b2 command several times, each time with only one thread and a different --with-[module] argument.

.. note::Make logs of bootstrap and b2 command
       
.. note:: 

   If you have made changes to your configuration it is important to clean your cache.
    - Navigate to your build folder.
    - Delete ``CMakeCache.txt``.
