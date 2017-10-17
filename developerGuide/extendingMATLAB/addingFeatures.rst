.. _extendingMATLAB_addingFeatures:

Adding features
===============

For each settings-containing class, a directory is added to :class:`tudatBundle/matlabInterface/MatlabInterface`. The name of this directory is irrelevant for the code to work, but it is preferable to the same name for the directory as for the base class. The name of the files do have to coincide with the name of the classes/enumerations. For instance, for the :class:`Integrator` class, we will create the file :class:`MatlabInterface/Integrator/Integrator.m`. In that same directory, we will add all the files for the derived classes (in this case only :class:`VariableStepSizeIntegrator.m`) as well as the files for the enumerations needed by these classes (in this case, :class:`Integrators.m` and :class:`RungeKuttaCoefficientSets.m`). Note that the name of the enumerations is in plural, and does not include the word :class:`Type` at the end. In the same way, the settings-containing classes do not add the word :class:`Settings` at the name: this is done to avoid repeating :class:`Settings` and :class:`Type` everywhere in the code unnecessarily, as in the MATLAB interface all the custom classes are used to define settings and all the enumerations to defined types, so it is not necessary to explicitly add this information to their names.

In some cases, when a class can only be contained by another class, it is recommended to put that class (and its derived classes and associated enumerations) in a directory inside the parent's class directory. For instance, the file :class:`MatlabInterface/Body/GravityField/GravityField.m` is in a directory inside the :class:`MatlabInterface/Body` directory because only :class:`Body` objects can contain :class:`GravityField` objects when setting up a valid Tudat simulation.


Enumerations
~~~~~~~~~~~~

An enumeration is defined by writing its possible values names. These names must match the string representations defined in the C++ code. Note that the unsupported enum values are also included, but commented out:

.. code-block:: matlab
  :caption: :class:`MatlabInterface/Acceleration/Accelerations.m`

  classdef Accelerations
      enumeration
          undefined
          pointMassGravity
          aerodynamic
          cannonBallRadiationPressure
          sphericalHarmonicGravity
          mutualSphericalHarmonicGravity
          % thirdBodyPointMassGravity
          % thirdBodySphericalHarmonicGravity
          % thirdBodyMutualSphericalHarmonicGravity
          thrust
          relativisticCorrection
          empirical
      end
  end



Settings classes
~~~~~~~~~~~~~~~~

Properties and methods
**********************

A settings-containing class defines a set of properties, whose names must match exactly those of the JSON keys. The properties that must not be used during conversion to JSON are marked as transient:

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
      properties (Transient, Dependent)
          dragCoefficient
          radiationPressureCoefficient
      end
      ...
  end

Dependent properties can be used as shortcuts to other properties. For instance, :literal:`Body.dragCoefficient` sets/gets the value of :literal:`Body.aerodynamics.dragCoefficient`, and :literal:`Body.radiationPressureCoefficient` sets/gets the value of :literal:`Body.radiationPressure.Sun.radiationPressureCoefficient`. Thus, we need to provide set/get methods. For instance:

.. code-block:: matlab
  :caption: :class:`MatlabInterface/Body/Body.m`

  classdef Body < jsonable
  ...
  methods
      ...
      function value = get.dragCoefficient(obj)
          value = obj.aerodynamics.dragCoefficient;
      end

      function set.dragCoefficient(obj,value)
          obj.aerodynamics.dragCoefficient = value;
      end
      ...
  end

The classes provided in the MATLAB Interface do not provide functionality, they are just classes containing settings that can be converted to JSON. Thus, they do not contain many methods (normally, just those overriding the methods of :class:`jsonable` to determine which properties are mandatory, paths, etc.). One exception is the :class:`Simulation` class, which does include some methods that add functionality, such as those for running the simulation, adding bodies, defining results to be exported to JSON files or to be imported from a Tudat propagation into MATLAB, etc.

For other classes, in addition to getters and setters for dependent properties and overriding the function needed to determine their JSON representation, we need to write the following methods:

- **Setters for properties storing enumeration values**, in order to enable using a char (in addition to an enumeration value) when defining it. For instance:
  
  .. code-block:: matlab
    :caption: :class:`MatlabInterface/Integrator/VariableStepSizeIntegrator.m`
    
    classdef VariableStepSizeIntegrator < Integrator
    properties
        ...
        rungeKuttaCoefficientSet
        ...
    end
    ...
    methods
        ...
        function set.rungeKuttaCoefficientSet(obj,value)
            if ~isa(value,'RungeKuttaCoefficientSets')
                value = RungeKuttaCoefficientSets(value);
            end
            obj.rungeKuttaCoefficientSet = value;
        end
    ...
    end
  
  In this way, these two assignments are equivalent and both result in a valid integrator:
  
  .. code-block:: matlab
    
    integrator.rungeKuttaCoefficientSet = RungeKuttaCoefficientSets.rungeKuttaFehlberg78;
    integrator.rungeKuttaCoefficientSet = 'rungeKuttaFehlberg78';
    
  Without writing the method :literal:`set.rungeKuttaCoefficientSet`, if the user writes a value that is not in the enumeration list, the second assignment does not generate any warnings or error, an invalid integrator will be generated and an error will be thrown when running Tudat. By writing the :literal:`set.rungeKuttaCoefficientSet` method, the second assignment will throw a MATLAB error if the provided enum value is not in the enumeration list, preventing the generation of invalid JSON files.
  
- **Getters for properties storing body names**, in order to be able to provide either body objects or body names. For instance:
  
  .. code-block:: matlab
    :caption: :class:`MatlabInterface/Acceleration/Thrust.m`
    
    classdef Thrust < Acceleration
        properties
            ...
            centralBody
        end
        
        methods
            ...
            function bodyName = get.centralBody(obj)
                if isa(obj.centralBody,'Body')
                    bodyName = obj.centralBody.name;
                else
                    bodyName = obj.centralBody;
                end
            end
            ...
        end
    end
    
  When getting the :class:`centralBody` property of a :class:`Thrust` object, a char will be returned even if the user assigned a :class:`Body` object to this property; in that case, the transient property :literal:`name` of the body will be used. Thus, the following two assignments are equivalent and valid:
  
  .. code-block:: matlab
    
    thrust.centralBody = 'Earth';
    thrust.centralBody = Earth;
    
  Note that :class:`Earth` is just a pre-defined :class:`Body` with the property :literal:`name` set to :literal:`'Earth'` and the property :literal:`useDefaultSettings` set to :literal:`true` (so that its properties are obtained automatically from Spice if not specified manually). Pre-defined objects for the other celestial bodies also exist.


Constructors
************

We can provide constructors for the settings-containing classes. Usually, in the classes provided in the MATLAB interface, a type containing an enumeration value is provided as first (and only) argument and assigned to a :literal:`type` property:

.. code-block:: matlab
  :caption: :class:`MatlabInterface/Acceleration.m`

  classdef Acceleration < jsonable
      properties
          type
      end
      
      methods
          function obj = Acceleration(type)
              if nargin >= 1
                  obj.type = type;
              end
          end
          ...
      end
  end

Providing constructors with additional arguments is discouraged, as this usually results in code that is not readily understandable. For instance, the following must be avoided:

.. code-block:: matlab
  
  integrator = VariableStepSizeIntegrator(10,'rungeKuttaFehlberg78',1,1000,1E-9,1E-9);

If the user writes this in their MATLAB script, it is difficult to know what each parameter represents. If one checks the constructor definition (in MATLAB, only one constructor can be provided per class), it will look like :literal:`VariableStepSizeIntegrator(varargin)` in most cases, which does not provide useful information. Moreover, if more properties are added in the future, these could only be added at the end of the list of input arguments; otherwise existing code would be broken (or worse, would still run but generating in wrong results). Thus, it is recommended to create instances of settings-containing objects using an empty constructor whenever possible:

.. code-block:: matlab
  
  integrator = VariableStepSizeIntegrator();
  integrator.initialStepSize = 10;
  integrator.rungeKuttaCoefficientSet = 'rungeKuttaFehlberg78';
  integrator.minimumStepSize = 1;
  integrator.maximumStepSize = 1000;
  integrator.relativeErrorTolerance = 1E-9;
  integrator.absoluteErrorTolerance = 1E-9;

In any case, **an empty constructor must always be provided** (i.e. the length of :literal:`varargin` must be allowed to be zero). In this way, we can write:

.. code-block:: matlab
  
  simulation = Simulation();
  ...
  simulation.integrator.stepSize = 20;
  json.export(simulation,'main.json');

because in the constructor of :class:`Simulation`, the property :literal:`integrator` is set to :literal:`Integrator()`. If we skip this step, the code above would still be valid, but the property :literal:`simulation.integrator` would be a :class:`struct` instead of an :class:`Integrator` object, so all the functionality defined in the :class:`Integrator` class for converting to JSON would be lost. This means that, all the properties storing settings-containing objects must be initialised in the constructor with their corresponding empty constructors to prevent this problems. For instance:

.. code-block:: matlab
  :caption: :class:`MatlabInterface/Acceleration/Thrust.m`

  classdef Thrust < Acceleration
      properties
          direction
          magnitude
          dataInterpolation
          specificImpulse
          frame
          centralBody
      end
      
      methods
          function obj = Thrust()
              obj@Acceleration(Accelerations.thrust);
              obj.direction = ThrustDirection();
              obj.magnitude = ThrustMagnitude();
              obj.dataInterpolation = DataInterpolation();
          end
          ...
      end
  end

Note that before we call the base class' constructor first by writing :literal:`obj@BaseClass(...)`.

For optional properties, we can also set them during construction, but in that case the used empty constructor **must** return an empty object. For instance, consider the constructor of the :literal:`Body` class:

.. code-block:: matlab

  function obj = Body(name,useDefaultSettings)
      obj.initialState = State();
      obj.aerodynamics = ConstantAerodynamics();
      obj.atmosphere = Atmosphere();
      obj.ephemeris = Ephemeris();
      obj.gravityField = GravityField();
      obj.gravityFieldVariation = GravityFieldVariation();
      obj.rotationModel = RotationModel();
      obj.shapeModel = ShapeModel();
      if nargin >= 1
          obj.name = name;
          if nargin >= 2
              obj.useDefaultSettings = useDefaultSettings;
          end
      end
  end

Thanks to this constructor, one can write:

.. code-block:: matlab

  body = Body('Vehicle');
  % body.aerdoynamics = ConstantAerodynamics();
  body.aerodynamics.forceCoefficients = [2 0 0];
  body.aerodynamics.referenceArea = 10;

without the need to write the line that has been commented out. However, :literal:`ConstantAerodynamics()` must provide an empty object, otherwise the property :literal:`body.aerodynamics` will be defined in the JSON file even when the user does not provide any aerodynamics settings. For :class:`Atmosphere`, :class:`Ephemeris`, :class:`GravityField`, etc. this is trivial, but for :class:`ConstantAerodynamics` this is more difficult because it derives from :class:`Aerodynamics` and its property :literal:`type` should be set to :literal:`'constant'`, resulting in an object that is not empty. Thus, this is only possible because the JSON Interface uses the default value :literal:`"constant"` for the type of the aerodynamics settings when no provided in the input file.
