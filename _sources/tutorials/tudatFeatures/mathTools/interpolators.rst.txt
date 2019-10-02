.. _tudatFeaturesInterpolators:

Interpolators
=============
For many applications, an interpolation routine will be an important part of your pre- and/or post-processing. Tudat comes with a number of interpolation routines, each of which is handled through the same class architecture.

As is the case for most other model that you create in Tudat, you provide settings to a create function to set up an interpolator. Before going into the details of this process, below is a list of the interpolation routines in Tudat:

	- One-dimensional:
	   - Linear interpolation
	   - Piece-wise constant interpolator
	   - Cubic spline interpolation
	   - Lagrange interpolation
	   - Hermite interpolation

	- Multi-dimensional:
	   - Multi-linear interpolation

Clearly, the first five options are interpolation routines with a single independent variable (i.e., they are one-dimensional), and the final one can be used for interpolation with an arbitrary number of independent variables. In this chapter, you will first introduced to how the one-dimensional interpolators are used, and then to how the multi-dimensional ones are implemented. These will be followed by a description of the various methods that can be used in case the requested independent variable for interpolation is outside the defined domain. 

One-dimensional Interpolator Creation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Creating interpolators of a single independent variable is most easily done in Tudat by creating an :class:`InterpolatorSettings` object, and passing this to the ``createOneDimensionalInterpolator`` function.

.. class:: InterpolatorSettings

   Class containing settings to setup an interpolator.

   .. code-block:: cpp

      // Load data.
      std::map< double, Eigen::Vector6d > stateMap;
      stateMap = ....
      
      // Create interpolator
      std::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings =
         std::make_shared< interpolators::InterpolatorSettings >( linear_interpolator ); 
      std::shared_ptr< OneDimensionalInterpolator< double, Eigen::Vector6d > > interpolator =
      	   interpolators::createOneDimensionalInterpolator(
      		   stateMap, interpolatorSettings );

      // Interpolate
      Eigen::Vector6d interpolatedResult = interpolator->interpolate( 0.0 );

   Note that if the data you interpolate is a different type, this should be reflected in how the interpolator is created. For instance, for ``std::map< double, Eigen::VectorXd >`` input, use

   .. code-block:: cpp

      // Load data.
      std::map< double, Eigen::VectorXd > stateMap;
      stateMap = ....

      // Create interpolator
      std::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings =
         std::make_shared< interpolators::InterpolatorSettings >( linear_interpolator ); 
      std::shared_ptr< OneDimensionalInterpolator< double, Eigen::VectorXd > > interpolator =
      	   interpolators::createOneDimensionalInterpolator(
      		   stateMap, interpolatorSettings );

      // Interpolate
      Eigen::VectorXd interpolatedResult = interpolator->interpolate( 0.0 );

   In this example, the ``stateMap`` contains the data that is interpolated, using the ``double`` key (time) as independent variable and the ``Eigen::Vector6d`` value (state) as dependent variable. The interpolation type is linear, and the ``interpolatorSettings`` object is created by passing only the argument ``linear_interpolator``.

The interpolator itself is of the type ``std::shared_ptr< OneDimensionalInterpolator< double, Eigen::Vector6d > >``, where the ``double`` and ``Eigen::Vector6d`` denote the (in)dependent variables of the interpolation. For you application they may be different (for instance :literal:`double` and :literal:`double`, or :literal:`long double` and :class:`Time`). The interpolator is created by calling the ``createOneDimensionalInterpolator`` function with the map containing the data and the interpolator settings as input.

The different interpolator types are handled in a similar manner:

- **Linear:** Requires no additional information to the ``createOneDimensionalInterpolator`` and is defined by using the :class:`InterpolatorSettings` base class (with ``linear_interpolator`` as argument).

- **Piece-wise constant:** Requires no additional information to the ``createOneDimensionalInterpolator`` and is defined by using the :class:`InterpolatorSettings` base class (with ``piecewise_constant_interpolator`` as argument).

- **Cubic spline:** Requires no additional information to the ``createOneDimensionalInterpolator`` and is defined by using the :class:`InterpolatorSettings` base class (with ``cubic_spline_interpolator`` as argument). 

   .. note:: Note that for the cubic spline implementation, natural boundary conditions are imposed (2nd derivatives at the boundaries equal to zero) and the first derivatives are continuous throughout the curve.

- **Hermite spline:** This interpolation requires the definition of the values, as well as the derivatives of the curve at each of the nodes. The interpolator is defined by using the ``InterpolatorSettings`` base class (with ``hermite_spline_interpolator``) as argument.
   
   .. code-block:: cpp

		// Load data for values at node points
		std::map< double, Eigen::Vector6d > stateMap;
		stateMap = ....

		// Load data for first derivatives at node points
		std::map< double, Eigen::Vector6d > stateDerivativeMap;
		stateDerivativeMap = ....

		// Load pair of default values for extrapolation (empty by default)
		std::pair< DependentVariableType, DependentVariableType > defaultExtrapolationValues;
		defaultExtrapolationValues = ...

		// Create interpolator
		std::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings =
		 std::make_shared< interpolators::InterpolatorSettings >( hermite_spline_interpolator ) 
		std::shared_ptr< OneDimensionalInterpolator< double, Eigen::Vector6d > > interpolator =
			 interpolators::createOneDimensionalInterpolator(
				 stateMap, interpolatorSettings, defaultExtrapolationValues, stateDerivativeMap );

		// Interpolate
		Eigen::Vector6d interpolatedResult = interpolator->interpolate( 0.0 );

Lagrange Interpolator
*********************
This interpolation routine uses an :math:`n`-th degree polynomial to approximate a function from :math:`n+1` data points. In our implementation, you can use a large data set of :math:`m` data points (with :math:`m > n`) to generate a set of interpolating polynomials. When interpolating a data point, the interpolation routine will automatically select the polynomial where the requested data point lies between the two middle points, to prevent wild oscillations (which occur at the edge of the polynomial). At the boundaries of the full interval, a cubic spline interpolator is used.

To create a Lagrange interpolator, the number of data points used for each interpolating polynomial should be defined, using the dedicated derived class :class:`LagrangeInterpolatorSettings`. The input argument for this class is the amount of points per polynomial. An example, when using 8 data points per polynomial, is described below:

.. class:: LagrangeInterpolatorSettings

   Derived class used for the settings of a lagrange interpolator.

   .. code-block:: cpp

	   // Load data for values at node points
	   std::map< double, Eigen::Vector6d > stateMap;
	   stateMap = ....

	   // Create interpolator
	   std::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings =
		   std::make_shared< interpolators::LagrangeInterpolatorSettings >( 8 ) 
	   std::shared_ptr< OneDimensionalInterpolator< double, Eigen::Vector6d > > interpolator =
			   interpolators::createOneDimensionalInterpolator(
				   stateMap, interpolatorSettings );

	   // Interpolate
	   Eigen::Vector6d interpolatedResult = interpolator->interpolate( 0.0 );

Multi-dimensional Interpolator Creation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Currently, only one multi-dimensional interpolator is supported in Tudat. This is the multi-linear interpolator, which is created by also using the :class:`InterpolatorSettings` class, and by specifying :literal:`multi_linear_interpolator` as first input. To create the interpolator object, a different function needs to be used. This is the :literal:`createMultiDimensionalInterpolator` function, whose inputs are as follows:

   - :literal:`independentValues`
      This is a vector of vectors, where each external vector corresponds to one independent variable, and the internal vector simply lists the independent variable points. 

   - :literal:`dependentData`
      This variable is of type :literal:`boost::multi_array`, which is very similar to a multi-dimensional matrix in MATLAB and is generally created by reading data from a file by using the :class:`MultiArrayFileReader` class.

   - :literal:`interpolatorSettings`
      This is simply a pointer to the :class:`InterpolatorSettings` class described above.

Interpolation When Independent Variables are Out-of-range
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In some instances, it can happen that the independent variable at which interpolation has to take place, falls outside the list of independent variables defined in the interpolator. In this case, the default behavior is extrapolation, but as you can imagine, this may lead to very wrong and strange behaviors, especially if a method like cubic spline is selected. For this reason, the user has the option of changing this so called, boundary handling method, or :literal:`boundaryHandling`, by specifying an extra input to the :class:`InterpolatorSettings` object. The allowed methods for independent variable out-of-range are the following:

   - :literal:`throw_exception_at_boundary`: throw a runtime error

   - :literal:`use_boundary_value`: use the closest defined value for these conditions

   - :literal:`extrapolate_at_boundary`: apply extrapolation (default behavior)

   - :literal:`use_default_value`: use the value provided by :literal:`defaultExtrapolationValue` (defined as an input to the function :literal:`createOneDimensionalInterpolator` or :literal:`createMultiDimensionalInterpolator`)

.. tip::
   Each of the methods above, with the exception of the first one (i.e., the one that gives the runtime error), comes with the option of also giving a warning in case the variable goes out-of-range. You just need to add :literal:`_with_warning` at the end of the enumeration. For example, you can use :literal:`extrapolate_at_boundary_with_warning` to extrapolate and give a warning. 

For a one-dimensional interpolator, the input should be a single :literal:`DependentVariableType` value, whereas for a multi-dimensional one, you can also specify one method per each independent variable (i.e., for each dimension). Unfortunately, it is not yet possible to define a :literal:`boundaryHandling` method for each extremity of the independent variable range (e.g., give an error in case an independent variable is requested below lowest value and extrapolate if beyond the highest value). 

In case you have specified :literal:`use_default_value`: or :literal:`use_default_value_with_warning`: as the boundary handling method, you should also specify the default value that the extrapolator has to output. As previously mentioned, this is done by adding an extra input to the :literal:`createOneDimensionalInterpolator` or :literal:`createMultiDimensionalInterpolator` functions. The format of this input, which is named :literal:`defaultExtrapolationValue`, depends on the dimensionality, and is as follows: 

   .. method:: std::pair< DependentVariableType, DependentVariableType >

      For **one-dimenisonal** interpolators, this variable is defined as a pair. The first value of the pair (accessed with the method :literal:`.first`) denotes the default extrapolation value in case the independent variable is requested below the lower domain boundary. The second value (:literal:`.second`), on the other hand, is output in case the independent variable is asked above the higher domain boundary.

   .. method:: std::vector< std::pair< DependentVariableType, DependentVariableType > >

      For **multi-dimensional** interpolators, it is defined as a vector of pairs. In this case, the pairs are used in the same manner as above, but a pair is defined for each independent variable. 
