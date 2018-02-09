.. _extendingMATLAB_basics:

Basics
======

No external packages or functions are used in the MATLAB interface. All the necessary code is available in the directory :class:`tudatBundle/matlabInterface/MatlabInterface` or through the use of built-in MATLAB functions. The MATLAB interface uses the built-in functions `jsondecode <https://nl.mathworks.com/help/matlab/ref/jsondecode.html>`_, `jsonencode <https://nl.mathworks.com/help/matlab/ref/jsonencode.html>`_ and `split <https://nl.mathworks.com/help/matlab/ref/split.html>`_, which are only available since version R2016b.


Code structure
~~~~~~~~~~~~~~

The contents of the MATLAB interface are structured in three directories:

  - **Examples**: contains example scripts.
  - **MatlabInterface**: contains the source code. Includes class definitions for creating settings objects, as well as some packages (directories starting with "+").
  - **UnitTests**: contains scripts used to generate input files for the :ref:`jsonInterface`.
  
Additionally, these three files can be found in the root directory:

  - **setup.m**: adds the directory :class:`tudatBundle/matlabInterface` to the user's MATLAB path definition permanently.
  - **build.m**: builds all the required targets (using CMake and the C++ code in the tudatBundle) for the MATLAB interface to work properly.
  - **tudat.m**: defines the class :class:`tudat` that is used to include all the source directories in the current session's path, amongst other things.


The :class:`tudat` class
************************

After running :literal:`setup.m` (to be done only once), one can write :literal:`tudat.load()` (in the current or future MATLAB sessions). In this way, all the classes and functions inside :class:`tudatBundle/matlabInterface/MatlabInterface` become available directly through their name, except for the code contained in package directories, which has to be accessed through :literal:`packageName.functionOrClassName(...)`.

It is also possible to run all the unit tests associated to the MATLAB/JSON interfaces by calling :literal:`tudat.test()`. This class also stores some paths, such as the absolute path to the MATLAB interface, its source directory, its unit tests directory or the (optional) file :literal:`settings.mat`, created when the user defines custom paths for the tudatBundle and associated targets. By default, the :class:`tudatBundle` directory is the parent directory of the :class:`matlabInterface` directory, but this can be changed by calling :literal:`tudat.find(absolutePathToTudatBundle)`.

The :class:`tudat` class can also be used to set/get user settings, stored in the file :literal:`settings.mat`:

.. code-block:: matlab
  :caption: :class:`MatlabInterface/tudat.m`

  properties (Constant, Hidden)
    ...
    bundlePathKey = 'bundlePath'
    defaultBundlePath = fileparts(tudat.rootdir);
    ...
  end
  ...
  methods (Static, Hidden)
    function path = bundle(newPath)
        if nargin == 0  % get
            try
                path = tudat.settings.(tudat.bundlePathKey);
            catch
                path = tudat.defaultBundlePath;
            end
        else  % set
            updateSetting(tudat.bundlePathKey,newPath);
        end
    end
  ...
  end
  
If the settings named :literal:`bundlePath` is not stored in the file :literal:`settings.mat`, or if this file does not exist, the default bundle path will be returned when calling :literal:`tudat.bundle`. The setting can be modified by calling :literal:`tudat.bundle(customPath)`. In the same way, support for additional user settings can be added in the future. Note that these properties and methods are usually hidden, as the user will rarely use them. The settings that can be currently customised are: :literal:`bundlePath`, :literal:`binaryPath` (the path to the binary :literal:`json_interface`), :literal:`testsSourcesDirectoryPath` (the path to the directory containing the JSON Interface unit tests C++ sources) and :literal:`testsBinariesDirectoryPath` (the path to the directory containing the C++ unit tests binaries).


Deserialisation (decoding)
~~~~~~~~~~~~~~~~~~~~~~~~~~

The built-in function :literal:`jsondecode` fails when the provided text to be decoded (or deserialised) contains line breaks. A custom :literal:`json.decode` function, located in :class:`MatlabInterface/+json/decode.m` has been written to bypass this limitation.

However, currently the MATLAB interface does not support creating settings objects from JSON text. This means that, if one decodes the string :literal:`{"mass": 5000, "referenceArea": 10}`, a struct will be obtained, but not an object of :class:`Body` class. To add support for creating settings objects from JSON text, it would be necessary to write custom constructor implementations or static functions for all the classes in the :class:`MatlabInterface` directory or to provide this for the :class:`jsonable` class, from which all the settings objects classes derive, in a similar way as done for serialisation.


Serialisation (encoding)
~~~~~~~~~~~~~~~~~~~~~~~~

Serialisation consists in converting MATLAB variables to JSON-formatted text. This is done using the built-in function :literal:`jsonencode`. This function generates a single-line text, which is difficult to read for long files. Thus, a custom :literal:`json.encode` function, located in :class:`MatlabInterface/+json/encode.m`, has been written. This function has a second input argument, :literal:`tabsize`, used to determine the number of spaces for each indentation level. For instance, a value of :literal:`2` adds line breaks and uses two space for each indentation level. When not specified, the value of :literal:`2` is used. When set to :literal:`0`, no spaces or line breaks are added, resulting in a single-line text.

Before encoding an object into JSON-formatted text, it is necessary to convert that object to something that is serialisable, basically to a number, string, boolean or array/map containing these objects. This is done by using the function :literal:`json.jsonize`. For numbers, strings and booleans, this function does nothing, as these types are directly serialisable. For structs and cell arrays, this function *jsonizes* each of their elements. The built-in encoding function only supports :class:`containers.Map` when its keys are strings, so :literal:`json.jsonize` extends this functionality by converting the keys of non-char :class:`containers.Map` to string by using the format specification :literal:`%g`. The most important addition of the :literal:`json.jsonize` is, however, support for custom classes and enumerations.

Enumerations
************

For enumerations, the name of the enumeration value is used. This means that nothing has to be added to the enumeration definition for it to be serialisable. For instance, consider the following enumeration:

.. code-block:: matlab
  :caption: :class:`MatlabInterface/MassRateModel/MassRateModels.m`

  classdef RotationModels
    enumeration
        simple
        spice
    end
  end

If we use the built-in function and write :literal:`jsonencode(RotationModels.simple)` we get:

.. code-block:: txt

  Error using jsonencode
  Unable to encode objects of enumeration class RotationModels as JSON-formatted text.

However, if we use the custom function and write :literal:`json.encode(RotationModels.simple)` we get the text :literal:`"simple"`. Writing :literal:`json.jsonize(RotationModels.simple)` returns the char variable :literal:`simple`, i.e. :literal:`json.jsonize` converts a variable to something that is serialisable, while :literal:`json.encode` actually serialises it.


Custom classes
**************

For classes, the serialisation process is a bit more complex. Although the built-in function supports encoding objects of custom classes, by converting them first to a struct in which each of its field names corresponds to a class' property, this fails when any of the properties is set to be an object that is not supported by :literal:`jsonencode`, such as an enumeration value. Additionally, in some cases it may not be necessary to use all the properties of a class in its JSON representation. Thus, for the classes used in the MATLAB interface, we define a :literal:`jsonize` method, which is then called by the :literal:`json.jsonize` function. This means that this function, as well as :literal:`json.encode`, do not support custom classes that do not implement the method :literal:`jsonize`.

Rather than implementing a custom :literal:`jsonize` method for all the custom classes, all the classes that should be convertible to JSON were made to derive from the class :class:`jsonable`, defined in :class:`MatlabInterface/jsonable.m`, which derives from the built-in class :class:`handle` (similar to a shared pointer, i.e. if a :class:`handle` is modified, all the objects containing it will also be updated automatically).

The :class:`jsonable` class has the following (hidden) methods:

- :literal:`getProperties`: returns the list of properties to be used for the JSON representation. This list includes all the properties of the (derived) class on which it is called **that are not marked as Transient**. For instance, for the :class:`Body` class:
  
  .. code-block:: matlab
    :caption: :class:`MatlabInterface/Body/Body.m`
  
    classdef Body < jsonable
      properties
          useDefaultSettings
          initialState
          mass
          rotationalState
          referenceArea
          aerodynamics
          atmosphere
          ephemeris
          gravityField
          gravityFieldVariation
          radiationPressure
          rotationModel
          shapeModel
      end
      properties (Transient)
          name
      end
    ...
    end
  
  when it is converted to JSON, the property :literal:`name` is ignored, because it is used internally by the MATLAB interface but it is not needed by the JSON interface or Tudat.


- :literal:`isMandatory(property)`: returns whether the property named :literal:`property` is mandatory. If a property is marked as mandatory, but its value is empty, the :literal:`json.encode` function will fail. In order to know which properties are mandatory, the derived classes have to implement the (hidden) method :literal:`getMandatoryProperties`. For instance, for :class:`Body`:
  
  .. code-block:: matlab
    :caption: :class:`MatlabInterface/Body/Body.m`
  
    classdef Body < jsonable
      ...
      methods (Hidden)
          function mp = getMandatoryProperties(obj)
              mp = {};
          end
      end
    end
  
  In this case, none of the properties of :class:`Body` are mandatory, because for a non-perturbed problem, we do not need to define any of its properties, so an empty body (i.e. :literal:`[]` or :literal:`null`) is valid. In other cases, some of the properties *are* mandatory:
  
  .. code-block:: matlab
    :caption: :class:`MatlabInterface/Propagator/Propagator.m`
  
    classdef Propagator < jsonable
      ...
      methods (Hidden)
          function mp = getMandatoryProperties(obj)
              mp = {'integratedStateType','bodiesToPropagate'};
          end
      end
    end
    
  .. code-block:: matlab
    :caption: :class:`MatlabInterface/Propagator/MassPropagator.m`
  
    classdef MassPropagator < Propagator
      ...
      methods (Hidden)
          function mp = getMandatoryProperties(obj)
              mp = getMandatoryProperties@Propagator(obj);
              mp = horzcat(mp,{'massRateModels'});
          end
      end
    end

  Note that, for :class:`MassPropagator`, the mandatory properties are those of :class:`Propagator`, from which it derives, i.e. :literal:`integratedStateType` and :literal:`bodiesToPropagate`, and also the property :literal:`massRateModels`. Thus, we get first the mandatory properties of the base class by calling :literal:`getMandatoryProperties@Propagator(obj)`, and then add additional mandatory properties for the derived class.


- :literal:`isPath(property)`: returns whether the property named :literal:`property` represents a path. Paths are automatically converted to JSON using the format specification :literal:`@path(%s)`. By default, this returns :literal:`false`. Derived classes can override this method:

  .. code-block:: matlab
    :caption: :class:`MatlabInterface/Body/Atmosphere/TabulatedAtmosphere.m`

    classdef TabulatedAtmosphere < Atmosphere
        properties
            file
        end
        ...
        methods (Hidden)
            function p = isPath(obj,property)
                p = strcmp(property,'file');
            end
        end
    end


- :literal:`isempty`: overrides the implementation of the :class:`handle` class, so that an instance of :class:`jsonable` will be empty if (an only if) all of its non-transient properties (i.e. the properties returned by :literal:`getProperties`) are empty. There is no need to override this method in derived classes of :class:`jsonable`.

- :literal:`jsonize`: returns an object that can be serialised as JSON-formatted text. Basically, it returns a struct containing the values for the non-transient properties (i.e. the properties returned by :literal:`getProperties`) that are non-empty. If an optional property has an empty value (e.g. :literal:`[]`), this field will not be defined in the returned struct. There is no need to provide a custom implementation of this method in derived classes; we only override this method if we want the JSON representation to include keys other than those associated to non-transient properties. This is done only once in the MATLAB interface in :class:`MatlabInterface/Termination/Termination.m`


Other types
***********

Some variables other than custom classes or enumerations are not serialisable either. For instance, complex numbers. These types have to be identified in the function :literal:`json.jsonize`, a custom implementation can be provided in this file:

.. code-block:: matlab
  :caption: :class:`MatlabInterface/+json/jsonize.m`
  
  function j = jsonize(object)
    ...
    try
      j = object.jsonize();
    catch ME
        ...
        elseif strcmp(ME.identifier,'MATLAB:structRefFromNonStruct')
            if isnumeric(object)
                if any(imag(object(:)))
                  ...
                

First, we try calling the :literal:`jsonize` method. If it is not defined, a :literal:`MATLAB:structRefFromNonStruct` error will be thrown. By catching this error, we can e.g. check whether the object is an imaginary number and provide a custom implementation for converting it to a value that is serialisable. For instance, the :ref:`jsonInterface` uses :literal:`boost::literal_cast` to cast strings to :class:`std::complex`, so we provide the required code so that :literal:`json.jsonize(0.5-4i)` returns the char :literal:`(0.5,-4)`, which is serialised as the string :literal:`"(0.5,-4)"`.

