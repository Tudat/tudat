.. _tudatFeaturesIntegrators:

Integrators
===========
.. warning:: The pages in this guide describe the software implementation of numerical integrators within Tudat, thus the reader is expected to be already familiar with the mathematics behind such tools.

Many applications require solving sets of first-order differential equations, where in the field of astrodynamics this is regularly done when using propagators, as described in :ref:`tudatFeaturesSimulatorIndex`. The purpose of a numerical integrator is to find a solution to a system of differential equations that satisfies a set of initial conditions, a problem commonly referred to in literature as an Initial-Value Problem (IVP). Tudat includes a framework for using numerical integrators, where several types are offered depending on how the step-size is discretized:

    - **Single-step methods:** These integrators use derivative information from a single step. Among these methods, Tudat includes the Euler method, the Runge-Kutta 4 method and several Runge-Kutta variable step-size methods.
    - **Multi-step methods:** These integrators use derivative information from multiple steps. Note that at the moment multi-step methods have not been implemented in Tudat.

Setting up a Numerical Integrator
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Creating an Integrator Object
*****************************
All numerical integrators presented in this page follow the same framework, where the first step is to create the :literal:`integrator` object, which has the :literal:`NumericalIntegratorXdPointer` type:

.. code-block:: cpp

    NumericalIntegratorXdPointer integrator
                    = boost::make_shared< EulerIntegratorXd >(
                        stateDerivativeFunction,
                        intervalStart,
                        initialState );

Although in this case :literal:`NumericalIntegratorXdPointer` is a *typedef* for a boost :literal:`shared_ptr`, the constructor can be created explicitly. To initialize an Euler integrator, the following arguments need to be specified:

- :literal:`EulerIntegratorXd` is the input argument for boost's :literal:`make_shared`, which defines the derived integrator class that we want to use.
- :literal:`stateDerivativeFunction` passes the state derivative function that defines the differential equation to be integrated. See below for more details.
- :literal:`intervalStart` provides the start of the integration interval and must have the type of the independent variable.
- :literal:`initialState` provides the initial state and must have the type of the dependent variable.
 
It is important to emphasize that the provided state derivative function must have a particular structure for the numerical integration to work. A sample state derivative function is defined here:

.. code-block:: cpp

    static inline Eigen::VectorXd stateDerivativeFunction( const double time,
                                                           const Eigen::VectorXd& state )
    {
        ...
    }

Any user-defined state derivative function must satisfy the following requirements:

- The first input argument provides the independent variable, which must be a :literal:`const double` and in this case is named :literal:`time`.
- The second input argument provides the dependent variable, which type is **free** to be defined. In this case, the type to be integrated is a :literal:`Eigen::VectorXd` which is passed by reference.
- The output argument type must coincide with the dependent variable type.

The contents of the state derivative function can be freely defined according to the user's needs, but it is important to respect the function's input-output structure.

Using the Integrator Object
***************************
Once the :literal:`integrator` object has been created, the multiple functions available within the :class:`NumericalIntegrator` class can be accessed. The following functions are available and they are exemplified using the :literal:`integrator` and the :literal:`stateDerivativeFunction` declared above:

- :literal:`getNextStepSize( )`

    Returns the step-size at the next integration step. The return type is defined by the type of the independent variable:

    .. code-block:: cpp

            double nextStepSize = integrator->getNextStepSize( );

- :literal:`getCurrentState( )`

    Returns the state at the current integration step. The return type is defined by the type of the dependent variable.

    .. code-block:: cpp

        Eigen::VectorXd currentState = integrator->getCurrentState( );

- :literal:`getCurrentIndependentVariable( )`

    Returns the value of the independent variable at the current integration step. The return type is defined by the type of the dependent variable.

    .. code-block:: cpp

        double currentIndependentVariable = integrator->getCurrentIndependentVariable( );

- :literal:`performIntegrationStep( stepSize )`

    Perform an integration step with the step-size fed to the first argument and return the integration result.

    .. code-block:: cpp

        Eigen::VectorXd stateAfterIntegrationStep = integrator->performIntegrationStep( nextStepSize );

- :literal:`integrateTo( intervalEnd , initialTimeStep )`

    Performs an integration until the specified :literal:`intervalEnd`, given the provided :literal:`initialTimeStep`. This function returns the state at the end of the interval.

    .. code-block:: cpp

        Eigen::VectorXd stateAtIntervalEnd = integrator->integrateTo( intervalEnd , initialTimeStep );

.. note:: The functions described above are virtual functions and thus redefined for each integrator method described in this page. Selection of the integrator method is made at the stage of creating the integrator object, where selection of the functions is taken care of by the implementation framework.

.. tip:: It is possible to integrate backwards in time by choosing an initial time step that is smaller then zero. If time is the chosen termination condition, the time at which the integration starts should also be larger then the final time.

Selecting a Numerical Integrator
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Euler Integrator
****************
The Euler integrator is the simplest integrator available but is also a first-order method, meaning that the global error is proportional to the step-size. Thus, its use for high-accuracy applications is not encouraged. Creating the Euler integrator object is done as follows:

.. code-block:: cpp

    NumericalIntegratorXdPointer integrator
                    = boost::make_shared< EulerIntegratorXd >(
                        stateDerivativeFunction,
                        intervalStart,
                        initialState );

Runge-Kutta 4 Integrator
************************
The Runge-Kutta 4 (RK4) integrator is a fourth-order fixed step-size integrator, thus performing better than the Euler integrator. The RK4 integrator is created as follows:

.. code-block:: cpp

    NumericalIntegratorXdPointer integrator
                    = boost::make_shared< RungeKutta4IntegratorXd >(
                        stateDerivativeFunction,
                        intervalStart,
                        initialState );

Runge-Kutta Variable Step-size Integrator
*****************************************
The Runge-Kutta variable step-size integrator involves a number of methods where the step-size is adjusted throughout the integration interval to bound the numerical error. Creating such integrators differs from the Euler integrator and the Runge-Kutta 4 fixed step-size methods:

.. code-block:: cpp

    RungeKuttaVariableStepSizeIntegratorXd integrator(
                rungeKuttaCoefficients,
                stateDerivativeFunction,
                initialTime,
                initialState,
                minimumStepSize,
                maximumStepSize,
                relativeErrorTolerance,
                absoluteErrorTolerance )

where the following arguments are necessary:

- :literal:`RungeKuttaVariableStepSizeIntegratorXd` is the input argument for boost's :literal:`make_shared`, which defines the derived integrator class that we want to use.
- :literal:`rungeKuttaCoefficients` provides the set of coefficients that define the particular variable step-size method being used. A number of Runge-Kutta coefficient sets are available in Tudat:

    .. code-block:: cpp
    
        // Runge-Kutta-Fehlberg 4(5)
        RungeKuttaCoefficients rungeKuttaCoefficients =
                RungeKuttaCoefficients::get( RungeKuttaCoefficients::rungeKuttaFehlberg45 );

        // Runge-Kutta-Fehlberg 5(6)
        RungeKuttaCoefficients rungeKuttaCoefficients =
                RungeKuttaCoefficients::get( RungeKuttaCoefficients::rungeKuttaFehlberg56 );

        // Runge-Kutta-Fehlberg 7(8)
        RungeKuttaCoefficients rungeKuttaCoefficients =
                RungeKuttaCoefficients::get( RungeKuttaCoefficients::rungeKuttaFehlberg78 );

        // Runge-Kutta-Fehlberg 8(7) Dormand-Prince
        RungeKuttaCoefficients rungeKuttaCoefficients =
                RungeKuttaCoefficients::get( RungeKuttaCoefficients::rungeKuttaFehlberg87 );

- :literal:`stateDerivativeFunction` passes the state derivative function that defines the differential equation to be integrated.
- :literal:`initialTime` provides the initial value of the independent variable.
- :literal:`initialState` provides the initial state and must have the type of the dependent variable.
- :literal:`minimumStepSize` defines the minimum step-size that the variable step-size integrator can take.
- :literal:`maximumStepSize` defines the maximum step-size that the variable step-size integrator can take.
- :literal:`relativeErrorTolerance` defines the relative error tolerance.
- :literal:`absoluteErrorTolerance` defines the absolute error tolerance.

Using a Numerical Integrator to Propagate an Orbit
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The numerical integrators described in this page are commonly used to propagate the orbit of spacecraft and celestial bodies. The reader is referred to :ref:`tudatFeaturesIntegratorSettings`, which discusses the how the numerical integrator fit in the simulator framework of Tudat.







