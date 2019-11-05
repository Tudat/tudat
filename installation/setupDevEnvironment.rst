.. _setupDevelopmentEnvironment:

Setup Development Environment
=============================

The first step is to set up the Integrated Development Environment (IDE) that will be used to write, compile and run your code. We recommend to use `Qt Creator <https://www.qt.io/ide/>`_ as your IDE due to its multi-platform capability and its popularity among students and staff. Note that the compiler and the precise installation steps will vary depending on your OS, so please refer to the section that is relevant to you. 

.. tip::
   
It is recommended to work with Tudat in a Unix-based environment, such as macOS or Ubuntu. This can speed up CMake parse times up to a factor of 10 and may also speed up the compilation process itself. Moreover, installation of the development environment is more straightforward compared to Windows. Tudat will work and install perfectly fine on Windows, but working with it, and especially developing it, will be less user-friendly. If you do not want to give up on your Windows installation, you can still benefit from Linux by exploiting a dual boot configuration; instructions on how to do this can be found `here <https://itsfoss.com/install-ubuntu-1404-dual-boot-mode-windows-8-81-uefi/>`_ for example.

.. toctree::
   :maxdepth: 1

   setupDevLinux
   setupDevWindows
   setupDevMacOs
