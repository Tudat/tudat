.. _basicJsonInterfaces:

.. role:: jsonkey

Basic JSON Interfaces
=====================

In this section you will be introductd to the JSON file format and to why it is used. Then, a brief description of how to use it for the purpose of Tudat applications will be given.

What Is JSON
~~~~~~~~~~~~

The term JSON stands for JavaScript Object Notation, and is a file format which, according to its Wikipedia page, uses human-readable text to transmit data objects consisting of attribute-value pairs and array data types. This is in fact, where one of JSON's main advantages can be found, which is explained in the next section. 

Writing a JSON file is a very simple process. The text you write will firstly have to be included in a series of curly brackets. These brackets usually enclose an object, which can include some member elements, which further describe the object. Some of these elements can be objects themselves. For instance, to describe a person, the Wikipedia page provides the following neat example:

.. code-block:: json

	{
	  "firstName": "John",
	  "lastName": "Smith",
	  "isAlive": true,
	  "age": 27,
	  "address": {
	    "streetAddress": "21 2nd Street",
	    "city": "New York",
	  },
	  "phoneNumbers": [
	    {
	      "type": "home",
	      "number": "212 555-1234"
	    },
	    {
	      "type": "office",
	      "number": "646 555-4567"
	    }
	  ],
	  "children": [],
	  "spouse": null
	}

Below, some of the main characteristics of how a JSON file is written will be outlined:

	- Note how the whole text is enclosed in two curly brackets. These denote the beginning and end of the file.
	- The object "person" is descibed by a series of characteristics (such as :jsonkey:`firstName`, :jsonkey:`isAlive`, :jsonkey:`phoneNumbers`, etc.) which are called *keys*. To each key corresponds a *value*, which can be of the following types:

		- **String**: text sorrounded by two quotation marks :literal:`"..."` (as in :jsonkey:`firstName`)
		- **Boolean**: either :literal:`true` or :literal:`false` (as in :jsonkey:`isAlive`)
		- **Number**: one or more digits (with no quotation marks) (as in :jsonkey:`age`)

			.. note:: In JSON there is no difference between integer and floating point representation. Every number is stored as a :literal:`double`.

		- **Object**: a series of key-value pairs, enclosed in curly brackets (as in :jsonkey:`address`)
		- **Array**: a series of values (of any type), enclosed in square brackets (as in :jsonkey:`phoneNumbers`)
		- :literal:`null`: an empty value (as in :jsonkey:`spouse`)

Why Use JSON
~~~~~~~~~~~~

Using JSON comes with a few practical advantages:

	- **Clear and concise structure**

		You will quickly realize while looking at the example applications in :ref:`jsonInterface_tutorials`, that writing an application in JSON is much clearer and more concise than with C++. In terms of clarity, the list of key and values gives a very neat overview of the simulation, and allows for quick additions/modifications to the code. In terms of conciseness, take the single unperturbed spacecraft example. The JSON file contains less than 70 lines, whereas the C++ version has slightly more than 150.

	- **No need to speak C++** (in most cases)

		You can use a JSON file as the main way of modeling your Tudat application. The file, which should contain the settings and conditions for all the objects/processes you want to simulate, is parsed by Tudat and the corresponding settings are automatically loaded. There may be, however, some cases where you still need to write a little code in C++. This is required, for instance, if you need to write an aerodynamic guidance algorithm, or to modify the default behavior of an object.

	- **No need to recompile after making modifications**

		In case you were to change your mind about some details, for example, the initial conditions of the spacecraft, there is no need to recompile the application in Qt Creator. Thus, instead of pressing the large "play" button in the bottom-left corner, you can simply press the smaller "play" button in the "Application Output" tab.

How to Use JSON
~~~~~~~~~~~~~~~

There are two ways of using the JSON interface with Tudat. In the first case, you can write a full numerical simulation in a JSON file which provides propagation output, whereas in the second option is to use a JSON file to load simulation settings and extend it with C++ code yourself. The latter is usefult if, for example, you need to write a custom guidance algorithm for your application. The second interface is explained in :ref:`jsonInterface_tutorials_lifetimeMaximisationCPP`.  The first option can be called in Terminal and Qt by following the procedure below.

The :literal:`json_interface` application can be used to run Tudat propagations by providing the path to a JSON input file as command-line argument in Terminal::

	cd // insert the path to the tudatBundle directory
	tudatBundle/tudat/bin/json_interface main.json

.. note:: The extension can be omitted, in which case :literal:`.json` will be assumed. The argument can also be the path to a directory containing a :literal:`main.json` file, and it can be omitted if a :literal:`main.json` file exists in the current directory.

It is also possible to use the :literal:`json_interface` application from Qt Creator. In the Projects > Run tab, choose :literal:`json_interface` in "Run configuration" and write the (absolute) path to your input file in "Command line arguments", as shown in the image below.

.. image:: qt.png

Note that details about the contents of the JSON file is explained in the upcoming :ref:`jsonInterface_tutorials` page.