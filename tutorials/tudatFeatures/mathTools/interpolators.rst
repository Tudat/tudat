.. _tudatFeaturesInterpolators:

Interpolators
=============
For many applications, an interpolation routine will be an important part of your pre- and/or post-processing. Tudat comes with a number of interpolation routines, each of which is handled through the same class architecture.

As is the case for most other model that you create in Tudat, you provide settings to a create function to set up an interpolator. Before going into the details of this process, below is a list of the interpolation routines in Tudat:

    - Linear single-variable interpolation
    - Cubic spline interpolation
    - Lagrange interpolation
    - Hermite interpolation
    - N-dimensional linear interpolation

The first four options are interpolation routines with a single independent variable. The final one can be used for interpolation from an arbitrary number of independent variables.

Single-variable interpolator creation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Creating interpolators of a single independent variable is most easily done in Tudat by creating an ``InterpolatorSettings`` object, and passing this to the ``createOneDimensionalInterpolator`` function.

.. code-block:: cpp

    // Load data.
    std::map< double, Eigen::Vector6d > stateMap;
    stateMap = ....

    // Create interpolator
    boost::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings =
        boost::make_shared< interpolators::InterpolatorSettings >( linear_interpolator ) 
    boost::shared_ptr< OneDimensionalInterpolator< double, Eigen::Vector6d > > interpolator =
            interpolators::createOneDimensionalInterpolator(
                stateMap, interpolatorSettings );

    // Interpolate
    Eigen::Vector6d interpolatedResult = interpolator.interpolate(0.0);

In this example, the ``stateMap`` contains the data that is interpolated, using the ``double`` key (time) as independent variable and the ``Eigen::Vector6d`` value (state) as dependent variable. The interpolation type is linear, and the ``interpolatorSettings`` object is created by passing only the argument ``linear_interpolator``.

The interpolator itself is of the type ``boost::shared_ptr< OneDimensionalInterpolator< double, Eigen::Vector6d > >``, where the ``double`` and ``Eigen::Vector6d`` denote the (in)dependent variables of the interpolation. For you application they may be different (for instance double and double, or long double and Time). The interpolator is created by calling the ``createOneDimensionalInterpolator`` function with the map containing the data and the interpolator settings as input.

The different interpolator types are handled in a similar manner:

    - **Linear:** Requires no additional information to the ``createOneDimensionalInterpolator`` and is defined by using the ``InterpolatorSettings`` base class (with ``linear_interpolator`` as argument).
    - **Cubic spline:** Requires no additional information to the ``createOneDimensionalInterpolator`` and is defined by using the ``InterpolatorSettings`` base class (with ``cubic_spline_interpolator`` as argument).

        .. note:: Note that for the cubic spline implementation, natural boundary conditions are imposed (2nd derivatives at the boundaries equal to zero) and the first derivatives are continuous throughout the curve.

    - **Hermite spline:** This interpolation requires the definition of the values, as well as the derivatives of the curve at each of the nodes. The interpolator is defined by using the ``InterpolatorSettings`` base class (with ``hermite_spline_interpolator`` as argument.

        .. code-block:: cpp
 
            // Load data for values at node points
            std::map< double, Eigen::Vector6d > stateMap;
            stateMap = ....

            // Load data for first derivatives at node points
            std::map< double, Eigen::Vector6d > stateDerivativeMap;
            stateDerivativeMap = ....

            // Create interpolator
            boost::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings =
                boost::make_shared< interpolators::InterpolatorSettings >( hermite_spline_interpolator ) 
            boost::shared_ptr< OneDimensionalInterpolator< double, Eigen::Vector6d > > interpolator =
                    interpolators::createOneDimensionalInterpolator(
                        stateMap, interpolatorSettings, stateDerivativeMap );

            // Interpolate
            Eigen::Vector6d interpolatedResult = interpolator.interpolate(0.0);

Lagrange interpolator
~~~~~~~~~~~~~~~~~~~~~
This interpolation routine uses an nth degree polynomial to approximate a function from (n+1) data points. In our implementation, you can use a large data set of m data points (with m > n) to generate a set of interpolating polynomials. When interpolating a data point, the interpolation routine will automatically select the polynomial where the requested data point lies between the two middle points, to prevent wild oscillations (which occur at the edge of the polynomial). At the boundaries of the full interval, a cubic spline interpolator is used.

To create a Lagrange interpolator, the number of data points used for each interpolating polynomial should be defined, using the dedicated derived class ``LagrangeInterpolatorSettings``. For example, when using 8 data points per polynomial, the following code should be used:

.. code-block:: cpp

    // Load data for values at node points
    std::map< double, Eigen::Vector6d > stateMap;
    stateMap = ....

    // Create interpolator
    boost::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings =
        boost::make_shared< interpolators::LagrangeInterpolatorSettings >( 8 ) 
    boost::shared_ptr< OneDimensionalInterpolator< double, Eigen::Vector6d > > interpolator =
            interpolators::createOneDimensionalInterpolator(
                stateMap, interpolatorSettings );

    // Interpolate
    Eigen::Vector6d interpolatedResult = interpolator.interpolate(0.0);

