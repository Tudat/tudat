.. _tudatTabulatedAtmosphere:

Tabulated Atmosphere Settings
=============================
The tabulated atmosphere can be customized to your needs. In particular, the generated atmosphere can vary from the very simple, such as an atmosphere that only depends on altitude and gives values of pressure, density and temperature to a fully specified atmosphere, dependent on altitude, longitude, latitude and time, and providing information on many atmospheric properties. 

In this chapter the main distinction is found between the altitude-only tabulated atmosphere and the more flexible multi-dimensional tabulated atmosphere. At the end of the chapter you will also find a list of built-in tabulated atmospheres and a brief description on how to use them.

.. tip::
   
   In :ref:`walkthroughsTabulatedAtmosphere` you can find an example where the tabulated atmosphere settings explained below are put into practice. Once you have read this chapter, we recommend trying out for yourself how this part of the environment works, by using the provided example.

.. _generalTabulatedAtmosphere:

General Tabulated Atmosphere
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In case a sophisticated tabulated atmosphere model is needed, one can use these descriptive settings. The full default constructor looks as follows:

.. code-block:: cpp

   std::map< int, std::string > atmosphereTableFile = ...
   std::vector< AtmosphereIndependentVariables > independentVariablesNames = ...
   std::vector< AtmosphereDependentVariables >& dependentVariablesNames = ...
   double specificGasConstant = ...
   double ratioOfSpecificHeats = ...
   std::vector< interpolators::BoundaryInterpolationType > boundaryHandling = ...
   std::vector< std::vector< std::pair< double, double > > > defaultExtrapolationValue = ...

   bodySettings[ "Earth" ]->atmosphereSettings = std::make_shared< TabulatedAtmosphereSettings >( atmosphereTableFile, 
                                                                                                     independentVariablesNames, 
                                                                                                     dependentVariablesNames, 
                                                                                                     specificGasConstant, 
                                                                                                     ratioOfSpecificHeats, 
                                                                                                     boundaryHandling, 
                                                                                                     defaultExtrapolationValue )

where :literal:`atmosphereTableFile` is a map of files containing information on the atmosphere, where the key of the map is an integer and the mapped value is the file name. Each file should correspond to a different dependent variable. The order of both independent and dependent parameters needs to be specified in the :literal:`independentVariablesNames` and :literal:`dependentVariablesNames` vectors, respectively. In :ref:`tudatWriteMultiDimensionalFiles` you can find a description of the format in which the files in :literal:`atmosphereTableFile` need to be.

.. note::

  The order of independent and dependent variables in the file is not important, however, it needs to match the order given by the :literal:`independentVariablesNames` and :literal:`dependentVariablesNames` variables, respectively. 

.. warning:: 

  The tabulated atmosphere needs to contain information on at least density, pressure and temperature.

The next two parameters, :literal:`independentVariablesNames` and :literal:`dependentVariablesNames`, specify, as mentioned above, the independent and dependent variables that are used in the atmosphere. For the independent variables, the following are supported:

  - :literal:`altitude_dependent_atmosphere`
  - :literal:`longitude_dependent_atmosphere`
  - :literal:`latitude_dependent_atmosphere`
  - :literal:`time_dependent_atmosphere`

whereas, for dependent variables, the supported variables are:

  - :literal:`density_dependent_atmosphere`
  - :literal:`pressure_dependent_atmosphere`
  - :literal:`temperature_dependent_atmosphere`
  - :literal:`gas_constant_dependent_atmosphere`
  - :literal:`specific_heat_ratio_dependent_atmosphere`
  - :literal:`molar_mass_dependent_atmosphere`

If the atmosphere you are modeling has constant gas constant and ratio of specific heats, then you can use the following two inputs, i.e., :literal:`specificGasConstant` and :literal:`ratioOfSpecificHeats`, to specify their values. However, keep in mind that these values will be disregarded if you have already specified :literal:`gas_constant_dependent_atmosphere` or :literal:`specific_heat_ratio_dependent_atmosphere` as dependent variables. In fact, another constructor exists where you can simply omit these two inputs, if you have already specified them.

Then, it's time to specify how the atmosphere will behave in case one or more of the independent variables goes out of range. It can happen that the altitude of the spacecraft is temporarely above the defined altitude of the tabulated atmosphere you provided. In this case, the default behavior is extrapolation. This means that Tudat will use the data provided to try and predict what the atmosphere outside its defined range looks like. As you can imagine, this can go very wrong sometimes. To avoid strange results, more options exist:

  - :literal:`throw_exception_at_boundary`: throw a runtime error

  - :literal:`use_boundary_value`: use the closest defined value for these conditions

  - :literal:`extrapolate_at_boundary`: extrapolate by using cubic spline in case of only one independent variable, or multi-linear interpolation in case more than one independent variable is defined. This is the **default** behavior.

  - :literal:`use_default_value`: use the value provided by the next input parameter, i.e., :literal:`defaultExtrapolationValue`.

.. tip::
  Each of the methods above, with the exception of the first one (i.e., the one that gives the runtime error), comes with the option of also giving a warning in case the variable goes out-of-range. You just need to add :literal:`_with_warning` at the end of the enumeration. For example, you can use :literal:`extrapolate_at_boundary_with_warning` to extrapolate and give a warning. 

You will need to specify the preferred behavior in the :literal:`boundaryHandling` parameter, which is the next input. Note that here, you can input the value in a few different formats:

  .. method:: As a single BoundaryInterpolationType

     The method you define will be used for each independent variable. 

  .. method:: As an std::vector< BoundaryInterpolationType >

     In this case, on the other hand, you can specify one out-of-range interpolation method for each independent variable. For instance, you may want to give an error in case the altitude value is out-of-range, but you may want to extrapolate in case time is requested outside its defined domain.

Unfortunately, it is not yet possible to define a :literal:`boundaryHandling` method for each extremity of the independent variable range (e.g., give an error in case altitude goes below lowest value and extrapolate if it goes beyond the highest value). 

If you decided to set :literal:`use_default_value` or :literal:`use_default_value_with_warning` as methods for :literal:`boundaryHandling`, then you will also need to specify the extrapolation value, called :literal:`defaultExtrapolationValue`. In the default constructor, the input format is as follows. For each dependent variable (e.g., temperature, gas constant, etc.), you need to define a vector, where each element in the vector corresponds to a pair of default values. A pair is defined for each independent variable (e.g., altitude, time, etc.). Moreover, each element in the pair (defined with the command :literal:`std::make_pair( X, Y )`), corresponds to the default value to be output in case the independent variable goes out-of-range at the lower or at the higher bound. For instance, if the altitude is defined between 100 and 1000 kilometers and the value at 50 km is requested, then the output value will be :literal:`X`, whereas if the value at 1500 km is requested, then the output will be :literal:`Y`. 

From the above description of the :literal:`defaultExtrapolationValue` it is clear that its use is very flexible. However, it can also be confusing and overly detailed. For this reason, a few extra constructors exist, where the :literal:`defaultExtrapolationValue` can be defined with a simpler format:

  .. method:: As a single DependentVariableType

     The same default value will be used for each dependent and independent variables and regardless of whether the lower or higher boundary has been overrun. 

  .. method:: As an std::vector< DependentVariableType >

     With this format, you can define a default value for each dependent variable, with will be used regardless of which independent variable goes out-of-range and whether the lower or higher boundary has been overrun. 

Altitude-dependent Tabulated Atmosphere
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Old Constructor
***************

Atmosphere model with properties (pressure, density, temperature) read in from a file, which is independent of time, latitude and longitude (i.e., only altitude is used to retrieve the atmospheric properties). 

.. code-block:: cpp

  std::string atmosphereFile = ...
  bodySettings[ "Earth" ]->atmosphereSettings = std::make_shared< TabulatedAtmosphereSettings >( atmosphereFile ); 

which will read the atmospheric properties from the file ``atmosphereFile`` (with as many columns as dependent variables).

.. note::

  These settings are compatible with the previous definition of the class :literal:`TabulatedAtmosphereSettings`, thus no change in your previous files is needed. 

New Constructor
***************

Similarly to the multi-dimensional tabulated atmosphere, you can also specify more parameters for the tabulated atmosphere. 

   .. code-block:: cpp

      // Define inputs
      std::string atmosphereFile = "...";

      std::vector< aerodynamics::AtmosphereDependentVariables > dependentVariables;
      dependentVariables.push_back( aerodynamics::pressure_dependent_atmosphere );
      dependentVariables.push_back( aerodynamics::temperature_dependent_atmosphere );
      dependentVariables.push_back( aerodynamics::specific_heat_ratio_dependent_atmosphere );
      dependentVariables.push_back( aerodynamics::density_dependent_atmosphere );

      specificGasConstant = 197.0;

      ratioOfSpecificHeats = 1.3;

      interpolators::BoundaryInterpolationType boundaryHandling = interpolators::use_default_value_with_warning;

      double defaultExtrapolationValue = TUDAT_NAN;

      // Set tabulated atmosphere
      bodySettings[ "Earth" ]->atmosphereSettings = std::make_shared< TabulatedAtmosphereSettings >( atmosphereFile, dependentVariables, specificGasConstant, ratioOfSpecificHeats, boundaryHandling, defaultExtrapolationValue );

For starters, you can choose from any of the dependent variables listed in the section above (note that even here, density, pressure and temperature need to be present). Moreover, you can specify custom values for the gas constant and ratio of specific heats (in the example above, values for Mars have been used). Then, the settings for interpolation can be defined. Their meaning is the same as for the multi-dimensional atmosphere, but only a single value can be input for both variables.

.. tip::

   If you want to create a one-dimensional tabulated atmosphere where the independent variable is not altitude (which is the default for the constructors above), you can always use one of the constructors in :ref:`generalTabulatedAtmosphere` and only input one file in the :literal:`std::map`.

Included Atmosphere Tables
~~~~~~~~~~~~~~~~~~~~~~~~~~

For Earth and Mars, one will find some built-in files that can be used to simulate the atmosphere at varying degrees of detail. These can be found in :literal:`Tudat/External/AtmosphereTables/`.

Earth
*****

Two atmosphere tables are available for Earth:

  - :literal:`USSA1976Until86kmPer100m.dat`: provides values of density, pressure and temperature for a range of altitudes from 0 to 86 kilometers with 100 meter spacing

  - :literal:`USSA1976Until100kmPer100mUntil1000kmPer1000m.dat`: same as above, but the range of altitudes is extended until 1000 kilometers, where the spacing is 100 meters until 100 kilometers, then becomes 1000 meters.

As the name suggests, both of these tables are based on the US Standard Atmosphere (USSA) 1976 model. 

Mars
****

Also Mars comes with two tabulated atmospheres:

  - :literal:`MCDMeanAtmosphere.dat`
     Provides values of density, pressure, temperature, gas constant, specific heat ratio and molar mass for a range of altitudes from 50 to 10000 kilometers with logarithmic spacing (i.e., the spacing goes from about 250 meters to 50 kilometers). This table was generated by averaging the atmosphere over longitude, latitude and time.

  - :literal:`MCDMeanAtmosphereTimeAverage`
     Provides the same depenent values as above (without molar mass), but as a function of altitude, longitude and latitude. In fact, the name above indicates a folder, not a file, where one can find 5 files, each corresponding to one dependent variable. The altitude range is the same, but the spacing, although still logarithmic, is larger (500 versus 1000 altitude elements). Spacing for longitude and latitude is constant at roughly 6 and 4 degrees, respectively. This table was generated by averaging the atmosphere over time.

These tables are based on the Mars Climate Database (MCD), and were generated via the web interface provided by ESA and others, which is available `here <http://www-mars.lmd.jussieu.fr/mcd_python/>`_.

.. _tudatWriteMultiDimensionalFiles:

Writing Multi-dimensional Files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The general format of multi-dimensional files is as follows::
	
	<number of independent variables> # block 1

	<independent variables> # block 2

	<dependent variables for each independent variable combination> # block 3

In the first block, one can find the number of independent variables. This will be used by the file reader to know how many lines of independent variables to look for, which are located in the second block. Then, in the final block, the dependent variables can be seen. Depending on how many independent variables there are, the storage format is adapted:

   - **Two-dimensional Files** (N - first dimension, M - second dimension)

      In the two-dimensional case, the data is stored as a matrix, where along the rows the values are stored according to the first independent variables, and along the columns they are stored according to the second. ::
   
         (1,1) ... (1,M)
           .         .
         (N,1) ... (N,M)

   - **Three-dimensional Files** (L - third dimension)

      The third dimension is added by stacking the two-dimensional slices vertically. This results to something similar to what is show below. ::
   
         (1,1,1) ... (1,M,1)
            .           .
         (N,1,1) ... (N,M,1)
         
            .           .

         (1,1,L) ... (1,M,L)
            .           .
         (N,1,L) ... (N,M,L)

   - **Four-dimensional Files** (P - forth dimension)

      The forth dimenison can be added in two methods. For the first one, which is shown below, the forth-dimension is added by stacking the three-dimensional slices horizontally: ::
   
         (1,1,1,1) ... (1,M,1,1)   .   (1,1,1,P) ... (1,M,1,P)
             .             .       .       .             .
         (N,1,1,1) ... (N,M,1,1)   .   (N,1,1,P) ... (N,M,1,P)
                       
             .             .       .       .             .

         (1,1,L,1) ... (1,M,L,1)   .   (1,1,L,P) ... (1,M,L,P)
             .             .       .       .             .
         (N,1,L,1) ... (N,M,L,1)   .   (N,1,L,P) ... (N,M,L,P)

      Alternatively, you can also decide to store the data in a similar manner as the thrid-dimension is stored, i.e., by stacking the forth-dimensional data points vertically: ::
   
         (1,1,1,1) ... (1,M,1,1)
             .             .
         (N,1,1,1) ... (N,M,1,1)
         
             .             .

         (1,1,L,1) ... (1,M,L,1)
             .             .
         (N,1,L,1) ... (N,M,L,1)
         
             .             .

         (1,1,1,P) ... (1,M,1,P)
             .             .
         (N,1,1,P) ... (N,M,1,P)
         
             .             .

         (1,1,L,P) ... (1,M,L,P)
             .             .
         (N,1,L,P) ... (N,M,L,P)

Future Additions
~~~~~~~~~~~~~~~~

   - Define a :literal:`boundaryHandling` method for each boundary (for each independent variable).
   - Input interpolator settings.