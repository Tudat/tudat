.. _tudatFeaturesHodographicShaping:

Hodographic shaping
===================

The hodographic shape-based method shapes the velocity of the spacecraft using cylindrical coordinates. Two different implementations of this shaping method exists, one using time as independent variable, and the other using polar angle. The first one has been implemented here.
This shaping method relies on the combination of integrable and differentiable base functions to shape the spacecraft velocity. The contribution of each of those base functions is weighted by a coefficient. Three base functions per velocity component at least are required so that their associated coefficients can be chosen to ensure the boundary conditions are fulfilled. The addition of any other base function, and thus of so-called free coefficients, introduces an extra degree of freedom in the problem, making it possible to optimise the shape-based trajectory to minimise the required deltaV.

Setting up a hodographically shaped trajectory
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A low-thrust trajectory designed by hodographic shaping can be created with the :literal:`HodographicShaping` class, which inherits from the base class :literal:`ShapeBasedMethodLeg`. A settings class for hodographic shaping has also been implemented (:literal:`HodographicShapingLegSettings`, see :ref:`tudatFeaturesSetUpLowThrustTrajectory` for more details) to construct :literal:`HodographicShaping` objects more easily.

.. class:: HodographicShaping

The :literal:`HodographicShaping` class is defined as follows:

.. code-block:: cpp
   
      HodographicShaping(
            const Eigen::Vector6d initialState,
            const Eigen::Vector6d finalState,
            const double timeOfFlight,
            const int numberOfRevolutions,
            simulation_setup::NamedBodyMap& bodyMap,
            const std::string bodyToPropagate,
            const std::string centralBody,
            std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > >& radialVelocityFunctionComponents,
            std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > >& normalVelocityFunctionComponents,
            std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > >& axialVelocityFunctionComponents,
            const Eigen::VectorXd freeCoefficientsRadialVelocityFunction,
            const Eigen::VectorXd freeCoefficientsNormalVelocityFunction,
            const Eigen::VectorXd freeCoefficientsAxialVelocityFunction,
            std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings =
                std::shared_ptr< numerical_integrators::IntegratorSettings< > >( ) )


where the inputs are:

	- :literal:`initialState`
		State of the spacecraft at departure.

	- :literal:`finalState`
		State of the spacecraft at arrival.

	- :literal:`timeOfFlight`
		Time of flight required for the shape-based trajectory (boundary condition automatically fulfilled from the way the hodographic shaping is implemented).

	- :literal:`numberOfRevolutions`
		Required number of revolutions before the spacecraft reaches its final state.

	- :literal:`bodyMap`
		Map of pointers to :literal:`Body` objects involved in the low-thrust trajectory.

	- :literal:`bodyToPropagate`
		Name of the spacecraft to be propagated.

	- :literal:`centralBody`
		Name of the central body of the low-thrust trajectory.

	- :literal:`radialVelocityFunctionComponents`
		Vector of :literal:`BaseFunctionHodographicShaping` objects, containing the definition of the base functions used to shape the radial velocity of the spacecraft during the low-thrust trajectory.

	- :literal:`normalVelocityFunctionComponents`
		Vector of :literal:`BaseFunctionHodographicShaping` objects, used to shape the normal velocity of the spacecraft during the low-thrust trajectory.

	- :literal:`axialVelocityFunctionComponents`
		Vector of :literal:`BaseFunctionHodographicShaping` objects, used to shape the axial velocity of the spacecraft during the low-thrust trajectory.

	- :literal:`freeCoefficientsRadialVelocityFunction`
		Vector of free coefficients for the radial velocity function (weighting coefficients for the additional base functions (any base function starting from the fourth one, if more than 3 base functions are provided in :literal:`radialVelocityFunctionComponents`)).

	- :literal:`freeCoefficientsNormalVelocityFunction`
		Vector of free coefficients for the normal velocity function (weighting coefficients for the additional base functions (any base function starting from the fourth one, if more than 3 base functions are provided in :literal:`normalVelocityFunctionComponents`)).

	- :literal:`freeCoefficientsAxialVelocityFunction`
		Vector of free coefficients for the axial velocity function (weighting coefficients for the additional base functions (any base function starting from the fourth one, if more than 3 base functions are provided in :literal:`axialVelocityFunctionComponents`)).

	- :literal:`integratorSettings`
		Integrator settings (empty by default), used to propagate the spacecraft mass or thrust profiles, or to numerically propagate the fully perturbed trajectory (as a means to assess the quality of the analytical shaped-based preliminary design).


Hodographic shaping base functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The constructor of the :literal:`HodographicShaping` class requires as inputs one vector of :literal:`BaseFunctionHodographicShaping` objects for each of the three velocity components. All the base functions to hodographically shape a trajectory are derived from the base class :literal:`BaseFunctionHodographicShaping`. The shaping functions for the three velocity components are not fixed in hodographic shaping, but rather are left to be selected by the user (to be chosen among a set of base functions). All these hodographic shaping base functions are defined in the classes presented below.  

.. class:: BaseFunctionHodographicShaping

This is the base class to derive the base functions used for hodographic shaping. Each class derived from this base class retains the following methods:
	
	- :literal:`evaluateFunction`
		Returns the value of the base function for a specific independent variable value.
	- :literal:`evaluateDerivative`
		Returns the derivative function of the base function at a specific independent variable value.
	- :literal:`evaluateIntegral`
		Returns the integral function of the base function at the specific independent variable.

The following classes inherit from the base class :literal:`BaseFunctionHodographicShaping`:

.. class:: ConstantFunctionHodographicShaping

Class for the constant function. The class is defined as follows and the function always returns the value 1.0:

.. code-block:: cpp
	
	ConstantFunctionHodographicShaping( ) 

.. class:: SineFunctionHodographicShaping

Class for the sine function. The class is defined as follows:

.. code-block:: cpp	
	
	SineFunctionHodographicShaping( double frequency )
	
where the input is the frequency of the sine function.
 
.. class:: CosineFunctionHodographicShaping

Class for the cosine function. The class is defined as follows:

.. code-block:: cpp	

	CosineFunctionHodographicShaping( double frequency )
	
where the input is the frequency of the cosine function.

.. class:: ExponentialFunctionHodographicShaping

Class for the exponential function. The class is defined as follows:

.. code-block:: cpp
	
    ExponentialFunctionHodographicShaping( double exponent )
	
where the input is the multiplying coefficient of the independent variable within the exponential function.

.. class:: ScaledExponentialFunctionHodographicShaping

Class for the scaled exponential function. The class is defined as follows:

.. code-block:: cpp	
	
	ScaledExponentialFunctionHodographicShaping( double exponent, double scaleFactor )
	
where the inputs are the multiplying coefficient of the independent variable within the exponential function, and a scaling factor to be applied in front of the multiplying coefficient.

.. class:: ExponentialSineFunctionHodographicShaping

Class for the exponential sine function. The class is defined as follows:
	
.. code-block:: cpp

	ExponentialSineFunctionHodographicShaping( double exponentExponentialFunction, double frequencySineFunction )

where the inputs are, in order, the exponent coefficient of the exponential function, and the frequency of the sine function. The exponential sine function returns the product of the sine and exponential functions.

.. class:: ScaledExponentialSineFunctionHodographicShaping

Class for the scaled exponential sine function. The class is defined as follows:

.. code-block:: cpp	

	ScaledExponentialSineFunctionHodographicShaping( double exponentExponentialFunction, double frequencySineFunction, double scaleFactor )
	
where the inputs are the exponent of the exponential function, the frequency of the sine function, and a scaling factor (to be applied to the exponential function, as implemented in :literal:`ScaledExponentialFunctionHodographicShaping`). This function returns the product of the scaled exponential function with the sine function.

.. class:: ExponentialCosineFunctionHodographicShaping

Class for the exponential cosine function. The class is defined as follows:

.. code-block:: cpp	

	ExponentialCosineFunctionHodographicShaping( double exponentExponentialFunction, double frequencyCosineFunction )
	
where the inputs are the exponent of the exponential function, and the frequency of the cosine function. This function returns the product of the cosine and exponential functions.

.. class:: ScaledExponentialCosineFunctionHodographicShaping

Class for the scaled exponential cosine function. The class is defined as follows:

.. code-block:: cpp
	
	ScaledExponentialCosineFunctionHodographicShaping( double exponentExponentialFunction, double frequencyCosineFunction, double scaleFactor )
	
where the inputs are the exponent of the exponential function, the frequency of the cosine function, and a scaling factor (to be applied to the exponential function, as implemented in :literal:`ScaledExponentialFunctionHodographicShaping`). This function returns the product of the scaled exponential function with the cosine function.

.. class:: PowerFunctionHodographicShaping

Class for the power function. The class is defined as follows:

.. code-block:: cpp
	
	PowerFunctionHodographicShaping( double exponent )
	
where the input parameter is the exponent of the power function. This function returns the value of the input parameter to the power indicated by the input variable.

.. class:: ScaledPowerFunctionHodographicShaping

Class for the scaled power function. The class is defined as follows:

.. code-block:: cpp
	
	ScaledPowerFunctionHodographicShaping( double exponent, double scaleFactor )
	
where the inputs are the exponent of the power function, and a scaling factor to be applied in front of this power function. This function returns the product between a scaling factor and the power function, as defined by the class :literal:`PowerFunctionHodographicShaping`.

.. class:: PowerSineFunctionHodographicShaping

Class for the power sine function. The class is defined as follows:

.. code-block:: cpp
	
	PowerSineFunctionHodographicShaping( double exponentPowerFunction, double frequencySineFunction )
	
where the inputs are the exponent of the power function and the frequency of the sine function. This function returns the product of the sine and the power functions (defined by the classes :literal:`SineFunctionHodographicShaping` and :literal:`PowerFunctionHodographicShaping`, respectively).

.. class:: ScaledPowerSineFunctionHodographicShaping

Class for the scaled power sine function. The class is defined as follows:

.. code-block:: cpp
	
	ScaledPowerSineFunctionHodographicShaping( double exponentPowerFunction, double frequencySineFunction, double scaleFactor )
	
where the inputs are the exponent of the power function, the frequency of the sine function, and a scaling factor (to be applied to the power function, as implemented in :literal:`ScaledPowerFunctionHodographicShaping`). This function returns the product of the scaled power function (:literal:`ScaledPowerFunctionHodographicShaping`) with the sine function (:literal:`SineFunctionHodographicShaping`).

.. class:: PowerCosineFunctionHodographicShaping

Class for the power cosine function. The class is defined as follows:

.. code-block::	cpp

	PowerCosineFunctionHodographicShaping( double exponentPowerFunction, double frequencyCosineFunction )
	
where the inputs are the exponent of the power function and the frequency of the cosine function. This function returns the product of the cosine and the power functions (:literal:`CosineFunctionHodographicShaping` and :literal:`PowerFunctionHodographicShaping`, respectively).

.. class:: ScaledPowerCosineFunctionHodographicShaping

Class for the scaled power cosine function. The class is defined as follows:

.. code-block:: cpp
	
	ScaledPowerCosineFunctionHodographicShaping( double exponentPowerFunction, double frequencyCosineFunction, double scaleFactor )
	
where the inputs are the exponent of the power function, the frequency of the cosine function, and a scaling factor (to be applied to the power function, as implemented in :literal:`ScaledPowerFunctionHodographicShaping`). This function returns the product of the cosine function (:literal:`CosineFunctionHodographicShaping`) with the scaled power function (:literal:`ScaledPowerFunctionHodographicShaping`).


Setting up the hodographic base functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To facilitate the creation of the base functions described above, settings classes have been implemented for all of them. A :literal:`BaseFunctionHodographicShaping` object can then be created from the corresponding settings with the function :literal:`createBaseFunctionHodographicShaping`. It takes the type of the base function, as well as a :literal:`BaseFunctionHodographicShapingSettings` object as inputs.

.. class:: BaseFunctionHodographicShapingSettings

This is the base class for base function settings for hodographic shaping. Several base function settings classes inherit from this base class:

.. class:: TrigonometricFunctionHodographicShapingSettings

Class of settings for cosine and sine base functions. It is defined as follows, taking the frequency of the trigonometric function as input:

.. code-block:: cpp

	TrigonometricFunctionHodographicShapingSettings( const double frequency )

.. class:: ExponentialFunctionHodographicShapingSettings

Class of settings for exponential functions (either simply exponential or scaled exponential functions). It is defined as follows (the scaled factor is by default set to 1):

.. code-block:: cpp

	ExponentialFunctionHodographicShapingSettings( const double exponent,
                                                   const double scaleFactor )

.. class:: ExponentialTimesTrigonometricFunctionHodographicShapingSettings

Class of settings for base functions obtained by multiplying exponential functions with trigonometric ones (This includes the following hodographic shaping base functions:  :literal:`ExponentialSineFunctionHodographicShaping`, :literal:`ExponentialCosineFunctionHodographicShaping`,
:literal:`ScaledExponentialSineFunctionHodographicShaping`, :literal:`ScaledExponentialCosineFunctionHodographicShaping`). The class is defined as follows (the scaled factor is by default set to 1.0):
	
.. code-block:: cpp

	ExponentialTimesTrigonometricFunctionHodographicShapingSettings( const double exponent,
                                                                     const double frequency,
                                                                     const double scaleFactor )

.. class:: PowerFunctionHodographicShapingSettings

Class of settings for power functions (either simply power or scaled power base functions). The class is defined as follows (the scaled factor is by default set to 1.0):

.. code-block:: cpp
	
	PowerFunctionHodographicShapingSettings( const double exponent,
                                             const double scaleFactor )

.. class:: PowerTimesTrigonometricFunctionHodographicShapingSettings

Class of settings for base functions defined by multiplying power functions with trigonometric ones (this includes the following hodographic shaping base functions: :literal:`PowerSineFunctionHodographicShaping`, :literal:`PowerCosineFunctionHodographicShaping`, :literal:`ScaledPowerSineFunctionHodographicShaping`, :literal:`ScaledPowerCosineFunctionHodographicShaping`). This settings class is defined as follows (the scaled factor is by default set to 1.0):

.. code-block:: cpp
	
	PowerTimesTrigonometricFunctionHodographicShapingSettings( const double exponent,
                                                               const double frequency,
                                                               const double scaleFactor )

.. _tudatFeaturesHodographicShapingOptimisation:

Optimising the hodographically shaped trajectory
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
	 
As explained before, hodographic shaping allows for a flexible number of base functions for each velocity component. Three base functions at least must be provided so that the boundary conditions can be satisfied, but a higher number of base functions can also be used, which turns into the creation of n-9 degrees of freedom (n being the total number of base functions provided, over the three velocity components). The weighting coefficients for these additional base functions are free parameters and can be tuned to minimise the deltaV required by the shaped trajectory. This thus transforms the hodographically shaping problem into an optimisation problem where the best set of free parameters, leading to the minimum deltaV, is to be found.

A pre-defined optimisation problem compatible with the PAGMO library has been implemented to this end.

.. class:: HodographicShapingOptimisationProblem

This class sets up the optimisation of the hodographically shaped trajectory, and its constructor is defined as follows:

.. code-block:: cpp

	HodographicShapingOptimisationProblem(
            Eigen::Vector6d initialState,
            Eigen::Vector6d finalState,
            const double timeOfFlight,
            const int numberOfRevolutions,
            simulation_setup::NamedBodyMap& bodyMap,
            const std::string bodyToPropagate,
            const std::string centralBody,
            std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > >& radialVelocityFunctionComponents,
            std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > >& normalVelocityFunctionComponents,
            std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > >& axialVelocityFunctionComponents,
            std::vector< std::vector< double > >& freeCoefficientsBounds )

where the input parameters are:

	- :literal:`initialState`
		State vector at departure

	- :literal:`finalState`
		State vector at arrival.

	- :literal:`timeOfFlight`
		Time-of-flight of the shaped trajectory.

	- :literal:`numberOfRevolutions`
		Expected number of revolutions of the shaped trajectory.

	- :literal:`bodyMap`
		Map of pointers to :literal:`Body` objects defining the trajectory environment.

	- :literal:`bodyToPropagate`
		Name of the spacecraft to be propagated.

	- :literal:`centralBody`
		Name of the central body of the trajectory.

	- :literal:`radialVelocityFunctionComponents`
		Vector of :literal:`BaseFunctionHodographicShaping` objects, containing the definition of the base functions used to shape the radial velocity.

	- :literal:`normalVelocityFunctionComponents`
		Vector of :literal:`BaseFunctionHodographicShaping` objects, containing the definition of the base functions used to shape the normal velocity.
		
	- :literal:`axialVelocityFunctionComponents`
		Vector of :literal:`BaseFunctionHodographicShaping` objects, containing the definition of the base functions used to shape the axial velocity.
	
	- :literal:`freeCoefficientsBounds`
		Vector containing the lower and upper bounds for the free coefficients of the hodographic shaping method.


The :literal:`fitness` function creates the hodographically shaped trajectory corresponding to the base functions provided as inputs, and to a given set of free parameters. It then returns the deltaV associated with this shaped trajectory.

The :literal:`get_bounds` function simply returns the bounds for the free coefficients to be optimised, which are already provided as inputs of the :literal:`HodographicShapingOptimisationProblem` constructor.	

 

	
