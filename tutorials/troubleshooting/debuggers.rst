.. _debuggers:

Debuggers
=========
Debuggers are essential in finding the source of your problems. They are able to provide you with information about where in your code a problem occurs. Several debuggers are available two of which are discussed here. 

Debugging in Qt
~~~~~~~~~~~~~~~
This build in debugger is available on all platforms. It can be started from within Qt by selecting ``Debug -> Start Debugging -> Start Debugging``. 

In the debugger console the output of the debugger can be seen. This often gives information on where your problem arises. The following is an example of the output of the Qt built in debugger for an error occuring in the FlightConditions update function. 

.. figure:: images/QtDebugger.png

Debugging in Terminal
~~~~~~~~~~~~~~~~~~~~~
On linux the debugger can be called from Terminal. It is installed by using (for debian-based distributions)::

   $ sudo apt-get update
   $ sudo apt-get install gdb

Then navigate to the folder with your application usually something like::

   $ cd tudatBundle/tudatApplication/myApplicationsFolder/bin/applications/

Then run gdb with the following command::

   $ gdb ./myApplicationName

Followed by::

   (gdb) run

This will run the application and stops when it stumbles upon debugging symbols. If your application crashes the error can be traced back at this point by using::

   (gdb) backtrace

When applied to the same example as for the Qt build-in debugger, gdb provide the following output in the terminal:

.. literalinclude:: DebuggingOutput.txt



