.. _extendingJSON_enhancedDeserialization:

.. role:: jsontype
.. role:: jsonkey

Enhanced deserialisation
========================

This is how a basic JSON input file for Tudat looks like:

.. literalinclude:: main.json
  :linenos:
  :language: json
  :caption: :class:`main.json`
  :name: main-json

The :class:`nlohmann::json` object created by parsing the contents of this file, referred to as :literal:`mainJson`, is used throughout the examples in this tutorial to explain the features that have been added to the :ref:`jsonInterface` regarding :ref:`extendingJSON_enhancedDeserialization`, :ref:`extendingJSON_enhancedValueAccess` and :ref:`extendingJSON_enhancedValueSet`.

Modular files
~~~~~~~~~~~~~

The JSON Interface library introduces a set of functions that can be used to parse modular JSON files. A modular JSON file is a JSON file in which special syntax is used to include (parts of) other JSON files (see :ref:`jsonInterface_modularFiles`). The following definitions are introduced:

  - **Declaration file**: file in which a JSON key and corresponding value are defined.
  - **Parent file**: file from which the declaration file is directly referenced.
  - **Root file**: file provided as input argument to the :literal:`json_interface` application.

The function :literal:`void parseModularJSON( nlohman::json& jsonObject, const boost::filesystem::path& filePath, ... )` declared in :class:`Tudat/JsonInterface/Support/deserialization.h` can be used to parse a modular JSON file. This function combines recursively the contents of all the files referenced from :literal:`jsonObject` obtained by parsing :literal:`filePath`, and updates the :literal:`jsonObject` (passed by reference).

Each of the referenced JSON files is parsed individually using the function :literal:`nlohmann::json readJSON( const boost::filesystem::path& filePath, ... )` declared in :class:`Tudat/JsonInterface/Support/deserialization.h`. In addition to parsing the contents of the file at :literal:`filePath` using the :literal:`nlohmann::json::parse` function, this function adds the following features:

  - If the file to be parsed contains a syntax error (i.e. invalid JSON syntax), the :literal:`nlohmann::json::parse` throws an error indicating the byte in which the syntax error is found. The :literal:`readJSON` catches this error and uses this information to throw an error indicating the line and column of the syntax error.
    
  - The file paths in the input file(s) are made relative to the root file's directory. To do so, we have to indicate that a given string represents a relative path by using the syntax :literal:`"@path(pathRelativeToDeclarationFile)"` which will be replaced to :literal:`"pathRelativeToRootFile"`.
    
  - The following substrings found inside **any** JSON string are replaced during run-time by:
    
    - :literal:`${FILE_DIR}`: absolute path of the directory where the declaration file is located.
    - :literal:`${FILE_STEM}`: filename (without extension) of the declaration file.
    - :literal:`${FILE_NAME}`: filename (with extension) of the declaration file.
    - :literal:`${PARENT_FILE_DIR}`: absolute path of the directory where the parent file is located.
    - :literal:`${PARENT_FILE_STEM}`: filename (without extension) of the parent file.
    - :literal:`${PARENT_FILE_NAME}`: filename (with extension) of the parent file.
    - :literal:`${ROOT_FILE_DIR}`: absolute path of the directory where the root file is located.
    - :literal:`${ROOT_FILE_STEM}`: filename (without extension) of the root file.
    - :literal:`${ROOT_FILE_NAME}`: filename (with extension) of the root file.


Mergeable files
~~~~~~~~~~~~~~~

In addition to modular JSON files, in which the content of (part of) other JSON files is included, the :literal:`json_interace` also introduces support for mergeable JSON files (see :ref:`jsonInterface_multicase`). Since the :literal:`json_interface` application expects a root JSON file containing an object, when providing a root JSON file containing an array, the contents of this array can be merged to generate a single JSON object to be used to set up the simulation. For instance, the following input file:

.. code-block:: json

  [
    "$(main.json)",
    {
      "bodies.asterix.initialState.eccentricity": 0
    }
  ]

can be merged leading to a JSON object identical to the one that would have been obtained by parsing :literal:`main.json`, but with the initial eccentricty of the body :literal:`asterix` set to :literal:`0`. The file must be first de-modularized (by replacing :literal:`"$(main.json)"` by an object retrieved from that file) and then merged, by (re-)defining the keys indicated in the second element of the array with the corresponding values. Thus, the returned :literal:`nlohmann::json` has value type :jsontype:`object`. Note that the following file would not result in the same merged :literal:`nlohmann::json` object:

.. code-block:: json

  [
    "$(main.json)",
    {
      "bodies": {
        "asterix": {
          "initialState": {
            "eccentricity": 0
          }
        }
      }
    }
  ]
  
since this would re-define the key :jsonkey:`bodies` of :literal:`main.json` to be an object containing only one element (the body :literal:`asterix`) whose only property would be an :literal:`initialState` with an :literal:`eccentricity` set to :literal:`0`.

In order to merge a :class:`nlohmann::json` of value type :jsontype:`array` into a one of value type :jsontype:`object`, the function :literal:`void mergeJSON( nlohmann::json& jsonObject, const boost::filesystem::path& filePath )` declared in :class:`Tudat/JsonInterface/Support/deserialization.h` is used. If the passed :literal:`jsonObject` is not of value type :jsontype:`array`, this function does nothing.


Full deserialisation sequence
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

All the features described previously are combined into the function :literal:`nlohmann::json getDeserializedJSON( const boost::filesystem::path& filePath )` declared in :class:`Tudat/JsonInterface/Support/deserialization.h`. This function replaces relative paths, combines modular files and merges objects when possible. This is the function that should be called when creating a :class:`nlohmann::json` object to be used to set up a simulation. In general, this function is only called once during the life-cycle of the application.

When testing individual parts of the JSON Interface library, the input :class:`nlohmann::json` object is not necessarily of value type :jsontype:`object`, and thus this function cannot be used, as the expected object may be of value type :jsontype:`array`. Thus, in :ref:`extendingJSON_unitTesting`, modular and mergeable JSON files are not used and the function :literal:`parseJSON` is used instead. In practice, for non-modular files, the only thing this function does is replacing strings such as :literal:`"@path(text)"` by :literal:`"text"`.

.. note:: Immediately before or after parsing a JSON input file (using either :literal:`getDeserializedJSON` or :literal:`parseJSON`), the current working directory must be changed to the root input file's directory using :literal:`boost::filesystem::current_path( rootFile.parent_path( ) )`, so that the relative paths defined in the obtained :class:`nlohmann::json` are valid.

