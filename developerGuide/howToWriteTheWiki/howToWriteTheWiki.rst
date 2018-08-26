.. _howToWriteTheWiki:

How to Write the Wiki
=====================
The goal of this page is to document how this wiki is written, with the goal of helping future users document their code and maintain the current content.

The Documentation Framework
~~~~~~~~~~~~~~~~~~~~~~~~~~~
The backbone of this wiki is reStructuredText, which is simply plain text with a special markdown syntax. ReStructuredText is written to :literal:`.rst` source files, which are then compiled to produce :literal:`.html` files.

The "compiler" used is Sphinx (http://www.sphinx-doc.org/en/stable/), a documentation generator widely used by the Python community but also useful for other software projects. Sphinx allows the user separate style work from documentation work, in a similar way that LaTeX documentation works. This wiki uses a style provided by Read the Docs (https://readthedocs.org/), again widely used by a number of software projects due to its beauty, flexibility and cross-compatibility with a number of vieweing devices.

Ultimately, the produced :literal:`.html` files can be hosted online and viewed under a web browser. This wiki is hosted under the Github Pages framework, which allows easy hosting of the :literal:`.html` files, an integrated backup solution and collaboration with git version control.

A top-level view of the documentation framework can be seen in the following figure:

.. figure:: images/wikiProcess.jpg   
   :align: center

Installing the Documentation Tools
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The documentation using reStructuredText and Sphinx requires relatively few tools and these are very easy to install. At the moment, documentation support is only given to users with a **Debian Linux** distribution or with **macOS**. If you have any other OS, you can make use of a Virtual Machine which will allow you to use a different OS without partitioning your drive. Please refer to Vmware Workstation Player (https://www.vmware.com/go/downloadplayer) or Oracle VirtualBox (https://www.virtualbox.org/wiki/Downloads) for virtual machine solutions.

**1. Install Python 2.7**
   Python 2.7 is necessary to run Sphinx. Please note that this may be already included in your system. You can install it using the following command::

      // Linux
      sudo apt-get install python

      // Mac
      sudo easy_install python // although Macs usually come with python built-in

**2. Install Sphinx**
   You can install the Sphinx documentation generator using the following commands::

      // Linux
      sudo apt-get install python pip
      pip install sphinx sphinx-autobuild

      // Mac
      sudo easy_install pip
      python -mpip install matplotlib --user
      pip install sphinx sphinx-autobuild --user

**3. Install Graphviz**
   The Graphviz software is used for creating the diagrams in the documentation, it is installed using::

      // Linux
      sudo apt-get install graphviz

      // Mac
      sudo easy_install graphviz

**4. Clone the** :literal:`source` **repository**
   This repository contains the :literal:`.rst` source files that make up the wiki. These are stored under the :literal:`source` branch within the :literal:`tudat` repository of https://github.com/Tudat/tudat. Place the contents of this repository inside a :literal:`source` folder. This can be done with the command::

      cd ... // path where you want the folder to be created
      git clone -b source https://github.com/Tudat/tudat.git source

**5. Clone the** :literal:`gh-pages` **repository**
   This repository contains the :literal:`.html` files that build the wiki website. These are stored under the :literal:`gh-pages` branch within the :literal:`tudat` repository of https://github.com/Tudat/tudat. Place the contents of this repository inside a :literal:`build/html` folder, next to your :literal:`source` folder. This can be done with the command::

      cd ... // path where you want the folder to be created
      git clone -b gh-pages https://github.com/Tudat/tudat.git build/html

**6. Create the** :literal:`make.bat` **and the** :literal:`Makefile` **files**
   These two files should be placed next to the folders. They allow to build the wiki from the terminal with few commands.

   Contents of the :literal:`make.bat`:

   .. code-block:: bat

      @ECHO OFF

      pushd %~dp0

      REM Command file for Sphinx documentation

      if "%SPHINXBUILD%" == "" (
      	set SPHINXBUILD=python -msphinx
      )
      set SOURCEDIR=source
      set BUILDDIR=build
      set SPHINXPROJ=TUDelftAstrodynamicsToolbox

      if "%1" == "" goto help

      %SPHINXBUILD% >NUL 2>NUL
      if errorlevel 9009 (
   	echo.
	echo.The Sphinx module was not found. Make sure you have Sphinx installed,
	echo.then set the SPHINXBUILD environment variable to point to the full
	echo.path of the 'sphinx-build' executable. Alternatively you may add the
	echo.Sphinx directory to PATH.
	echo.
	echo.If you don't have Sphinx installed, grab it from
	echo.http://sphinx-doc.org/
	exit /b 1
      )

      %SPHINXBUILD% -M %1 %SOURCEDIR% %BUILDDIR% %SPHINXOPTS%
      goto end

      :help
      %SPHINXBUILD% -M help %SOURCEDIR% %BUILDDIR% %SPHINXOPTS%

      :end
      popd

   .. note ::
      Note that for macOS, the :literal:`make.bat` file is not necessary.

   Contents of the :literal:`Makefile`:

   .. literalinclude:: Makefile.txt
      :language: makefile

.. tip ::
   Make sure that :literal:`Makefile` does not have an extension. On macOS, for instance, if you use TextEdit to create the file, it is possible that the file is automatically saved as :literal:`.txt` (with hidden extension).

**7. Compiling the wiki**
   To compile the wiki, navigate using the terminal to the folder containing the :literal:`make.bat` file and launch the following command::

      make html

  This will compile the :literal:`.html` files which you can view using your favorite web browser.

   The overall files structure should look as follows::

      root
      | 
      | build
      |     |
      |     | html
      |          |
      |          | index.html
      |          | ...
      |   
      | source
      |      |
      |      | conf.py
      |      | index.rst
      |      | ...
      |
      | make.bat
      |
      | makefile

Updating Documentation
~~~~~~~~~~~~~~~~~~~~~~
If changes are made to the documentation (wiki and/or doxygen) . These changes need to be pushed to the correct repository. 

**Wiki**

For the wiki the contents of ``root/build/html/`` need to be pushed to https://github.com/Tudat/tudat/tree/gh-pages. 

**Doxygen**

For the Doxygen the contents of ``tudatBundle/tudat/Tudat/Documentation/`` need to be pushed to https://github.com/Tudat/tudatDoxygen. 

Wiki Style Guide
~~~~~~~~~~~~~~~~
As mentioned above, the wiki is built using :literal:`.rst` files that contain a special mark-up. This files are linked together using :literal:`toctree` commands, that link the source files together.

Linking :literal:`.rst` Files
*****************************
The wiki source files are organized in different folders, where each folder contains an :literal:`index.rst` file that serves are the header of such folder. The source files should be organized in folder with meaningful names::

      source
      | 
      | index.rst
      |
      | installation
      |            |
      |            | index.rst
      |            | downloadTudatBundle.rst
      |            | downloadTudatBundle
      |            |                   |
      |            |                   | commandLine.rst
      |            |                   | smartgitDeveloper.rst
      |            |                   | smartgitUser.rst	
      |            | ...
      |
      | tutorials
      |         |
      |         | index.rst
      |         | prerequiredKnowledge
      |         |                    | 
      |         |                    | index.rst
      |         |                    | ...
      |         | ...
      |        
      | developerGuide
      |              | index.rst
      |              | ...
      | ...


The files are then linked together in a top-down fashion from the top :literal:`index.rst` file. Such linking is done by means of :literal:`toctree` commands, which are placed within the source files:

   .. code-block:: rst
      
      .. toctree::
         :maxdepth: 2
         :hidden:
         :caption: Contents

         self
         installation/index
         tutorials/index
         developerGuide/index

   where:

   - :literal:`.. toctree::`
      Launches the :literal:`toctree` command. Here we see a special syntax widely used in reStructuredText which are used to launch special commands.
   - :literal:`:maxdepth: 2`
      Defines the maximum number of levels to show in the :literal:`toctree`. Note that the corresponding source pages will still be linked, even if not shown. This option becomes irrelevant if :literal:`hidden` is used.
   - :literal:`:hidden:`
      Hides the :literal:`toctree` entirely from the page.
   - :literal:`:caption:`
      Defines the name of the :literal:`toctree`. In the Tudat wiki, this is only used in the top-level :literal:`toctree`.
   - :literal:`self`, :literal:`installation/index`, ...
      Define the location of the corresponding source files to be linked in the :literal:`toctree`. 

.. note:: Having an :literal:`index.rst` file in each folder is done simply due to convention. You could use a different file name as long as you refer to it properly in the :literal:`toctree`.

Indentation in reStructuredText
*******************************
Once the wiki source files have been linked, it is time to start writing your text. While doing so, it is highly emphasized to nest all similar content among different categories and captions. To do so, you need to use the indentation system of reStructuredText, which consist on sequentially adding blocks of 3 whitespaces depending on the indentation level. For instance:

.. code-block:: rst

   First level
   ===========
      Bla bla bla starts after 3 whitespaces.

      Second level
      ~~~~~~~~~~~~
         Bla bla bla starts after 6 whitespaces

         Third level
         ***********
         Bla bla bla starts after 9 whitespaces

Where all text with an indentation of at least 3 whitespace belongs to the caption "First level", all text with an indentation of at least 6 whitespace belongs to the caption "Second level" and all text with an indentation of at least 9 whitespace belongs to the caption "Third level".

.. warning:: It is important to realize that the "3-whitespaces" rule is applicable to all commands within reStructuredText. In fact, many commands require to add all the relative content to such command with an indentation of 3 whitespaces for the command to actually work. For instance:

   .. code-block:: rst

      .. code-block:: cpp

         This text WILL be recognized within code-block.

      .. code-block:: cpp

      This text WILL NOT be recognized within code-block.

Use of Special Commands
***********************
One of the most powerful features of reStructuredText is the ability to use complex mark-up. It is required that you make use of it to ensure that the wiki stays pretty and clear. The best advice we can give you is to go through existing wiki pages and cross-check against the source code that generates them. You can do so by clicking on "View page source" on the top-right corner of the page.

The following commands are available within reStructuredText

- :literal:`.. warning::`, :literal:`.. note::`, :literal:`.. tip::`

   These commands place coloured boxes which highlight particular content relevant to the current page. :literal:`.. warning::` should be used for critical information will likely lead to issues if ignored. :literal:`.. note::` should be used to extend the information beyond the regular paragraphs for further clarification, but is not critical for the end-user. :literal:`.. tip::` should be used for non-critical information that can help the user in solve further problems than the one treated.

- :literal:`.. code-block:: cpp`

   Sets a box with a monospaced box that contains C++ (cpp) marked-up code. You can use alternative designators such as reStructuredText (rst), MakeFile (makefile) or Windows Batch files (bat), among others.

- :literal:`**this text goes in bold**`
   Makes the text enclosed by double multiplication signs **bold**.

- :literal:`*this text goes in italics*`
   Makes the text enclosed by single multiplication signs *italics*.

- :literal:`:literal:```
   Introduces text in a red monospaced font. Should be used to refer to "software objects", such as clickable buttons in GUIs, files, folders, variables, etc.

- :literal:`:class:```
   Introduces text in a black monospaced font. Should be used to refer to C++ classes. In Tudat, these are generally identified by variable names where each word is capitalized, such as in :class:`ExampleClass`.

- :literal:`.. _exampleReferenceLink:`
   Introduces a referable link, typically used at the start of sections. Similarly to cross-links used in LaTeX.

- :literal:`:ref:`exampleReferenceLink``
   Introduces an inline link to the page designated with :literal:`.. _exampleReferenceLink:`.

- :literal:`.. class:: ExampleClass`
   Creates a special environment for class objects. If a :literal:`:class:`ExampleClass`` exists within the wiki, it will exist as a clickable object that redirects the user to the environment defined here.

- :literal:`use of semi-colon after a paragraph::`
   Introduces a box of monospaced text without mark-up. Generally used to define file locations.


Upcoming Improvements
~~~~~~~~~~~~~~~~~~~~~

The following issues on the `github project page <https://github.com/Tudat/tudat/labels/documentation>`_ need to be resolved. Other improvements to the documentation are:

- High priority.
  
- Medium priority.

   - Update Qt project tree page to new layout occuring in latest Qt release.
   - Move gh-pages and source branch of tudat repo to seperate repo.

- Low Priority.
 
   - Explain important functions in Basics/utilities.h
   - Classes which need to be defined:
      GroundStationState, PPNParameterSet, BasicSolidBodyTideGravityFieldVariations, ConstantObservationBias, ConstantRelativeObservationBias, JsonSimulationManager, MAssPropagationSettings, CustomStatePropagationSettings, MultiTypePropagationSettings, MultiArcPropagationSettings
   - Fix large literals which go beyond page margins such as in `Defining estimation input <http://tudat.tudelft.nl/tutorials/tudatFeatures/estimationSetup/estimationExecution.html#defining-estimation-input>`_. 
   - Extensively refer to the Tudat libraries in the JSON interface :ref:`jsonInterface_keys`.
   - Update the definition of the observation bias/bias estimation.
   - Make the differences between Windows and *NIX instructions clear in <http://tudat.tudelft.nl/installation/troubleshooting/compilingBoostManually.html>_.

