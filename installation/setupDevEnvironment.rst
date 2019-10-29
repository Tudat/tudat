.. _setupDevelopmentEnvironment:

Setup Development Environment
=============================

The first step is to set up the Integrated Development Environment (IDE) that will be used to write, compile and run your code. We recommend to use `Qt Creator <https://www.qt.io/ide/>`_ as your IDE due to its multi-platform capability and its popularity among students and staff. Note that the compiler and the precise installation steps will vary depending on your OS, so please refer to the section that is relevant to you. 

.. tip::
   
   As the build system (CMake and GCC) used for Tudat is native in Linux, it is recommended to work with Tudat in a Linux environment, such as Ubuntu. This can speed up CMake parse times up to a factor of 10 and may also speed up the compilation process itself. Moreover, installation of the development environment is more straightforward on Linux as compared to Windows. If you do not want to give up on your Windows installation, you can still benefit from Linux by exploiting a dual boot configuration; instructions on how to do this can be found **here**.

.. toctree::
   :maxdepth: 1

   setupDevLinux
   setupDevWindows
   setupDevMacOs
