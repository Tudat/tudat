.. _tudatFeaturesFilters:

Filters
=======
Filters, in astrodynamics, are generally used to estimate the state of a system, based on a model description and measurement data. Tudat provides a few filters, based on the Kalman filter concept, that you can use for these or any other purpose you see fit. The filtering techinques currently offered in Tudat are:

   - **Linear Kalman Filter** (LKF): an extremely effective and versatile procedure for combining noisy sensor outputs to estimate the state of a linear system with uncertain dynamics
   - **Extended Kalman Filter** (EKF): extends the LKF to non-linear applications, by using a first-order Taylor series approximation of the system around the current estimate
   - **Unscented Kalman Filter** (UKF): similar to the EKF, but addresses the non-linearity by using a deterministic sampling approach (instead of using system and measurement Jacobians)

You can find more information on the mathematical background of these filters in virtually any guidance, navigation and control book. However, the references used to define the filters in Tudat are:

   - **LKF**: E. Mooij, AE4870B - Re-entry Systems, Lecture Notes, Delft University of Technology, 2016.
   - **EKF**: Ogata, K., Discrete-Time Control Systems, 2nd ed. Pearson Education Asia, 2002.
   - **UKF**: Wan, E. and Van Der Merwe, R., “The Unscented Kalman Filter for Nonlinear Estimation,” in Adaptive Systems for Signal Processing, Communications, and Control Symposium. Institute of Electrical and Electronics Engineers, 2000, pp. 153–158.

On this page, you will find a description of how to create a Kalman filter object, either based on the linear, extended, or unscented principle. In general, you will notice that the creation of a filter, requires the input of two template arguments. These are called :literal:`IndependentVariableType` and :literal:`DependentVariableType`, and they respectively denote the type (i.e., :literal:`double`, :literal:`long double`, etc.) for the independent variable (usually time) and the dependent variables. Note that the dependent variable is internally defined as an :literal:`Eigen::Matrix< DependentVariableType, Eigen::Dynamic, 1 >` object (which corresponds to a vector of unspecified length). 

Linear Kalman Filter
~~~~~~~~~~~~~~~~~~~~

.. class:: LinearKalmanFilter

   .. code-block:: cpp

      boost::shared_pointer< LinearKalmanFilter< IndependentVariableType, DependentVariableType > > linearFilter = 
         boost::make_shared< LinearKalmanFilter< IndependentVariableType, DependentVariableType > >(
                     stateTransitionMatrix,
                     controlMatrix,
                     measurementMatrix,
                     systemUncertainty, 
                     measurementUncertainty,
                     filteringTimeStep,
                     initialTime, 
                     initialEstimatedStateVector, 
                     initialEstimatedStateCovarianceMatrix );

   In the code above, the inputs are given by:

      - :literal:`stateTransitionMatrixFunction`: a function returning the state transition matrix, as a function of time, state and control input
      - :literal:`controlMatrixFunction`: a function returning the control matrix as a function of time, state and control input
      - :literal:`measurementMatrixFunction`: a function returning the measurement matrix as a function of time, state and control input
      - :literal:`systemUncertainty`: a matrix defining the uncertainty in modeling of the system
      - :literal:`measurementUncertainty`: a matrix defining the uncertainty in modeling of the measurements
      - :literal:`filteringTimeStep`: a scalar representing the value of the constant filtering time step
      - :literal:`initialTime`: a scalar representing the value of the initial time
      - :literal:`initialStateVector`: a vector representing the initial (estimated) state of the system; it is used as first a-priori estimate of the state vector
      - :literal:`initialCovarianceMatrix`: a matrix representing the initial (estimated) covariance of the system; it is used as first a-priori estimate of the covariance matrix
      - :literal:`integrator`: a pointer to integrator to be used to propagate state

   .. note:: Currently, the creation of the LKF is only supported in the method above. A class such as :class:`LinearKalmanFilterSettings` has not yet been implemented.

Extended Kalman Filter
~~~~~~~~~~~~~~~~~~~~~~

.. class:: ExtendedKalmanFilter

   .. code-block:: cpp

      boost::shared_ptr< FilterBase< IndependentVariableType, DependentVariableType > > filterObject;
      filterObject = createFilter( filterSettings,
                                   systemFunction,
                                   measurementFunction,
                                   stateJacobianFunction,
                                   stateNoiseJacobianFunction,
                                   measurementJacobianFunction,
                                   measurementNoiseJacobianFunction );

   Here, the inputs are as follows:

      - :literal:`filterSettings`: a pointer to a :class:`ExtendedKalmanFilterSettings` object, described at the end of this section
      - :literal:`systemFunction`: a function returning the state as a function of time and state vector; this can be a differential equation if the :literal:`integratorSettings` is set (i.e., if it is not a :literal:`nullptr`)
      - :literal:`measurementFunction`: a function returning the measurement as a function of time and state
      - :literal:`stateJacobianFunction`: a function returning the Jacobian of the system w.r.t. the state; its input values can be time and state vector
      - :literal:`stateNoiseJacobianFunction`: a function returning the Jacobian of the system function w.r.t. the system noise; its input values can be time and state vector
      - :literal:`measurementJacobianFunction`: a function returning the Jacobian of the measurement function w.r.t. the state; its input values can be time and state vector
      - :literal:`measurementNoiseJacobianFunction`: a function returning the Jacobian of the measurement function w.r.t. the measurement noise. The input values can be time and state vector

   The filter settings above, can be defined for an extended Kalman filter as:

   .. class:: ExtendedKalmanFilterSettings

      .. code-block:: cpp

         boost::shared_ptr< FilterSettings< IndependentVariableType, DependentVariableType > > filterSettings = 
            boost::make_shared< ExtendedKalmanFilterSettings< IndependentVariableType, DependentVariableType > >( 
                        systemUncertainty,
                        measurementUncertainty,
                        filteringTimeStep,
                        initialTime,
                        initialStateVector,
                        initialCovarianceMatrix,
                        integratorSettings );

      where each element is defined as:

         - :literal:`systemUncertainty`: a matrix defining the uncertainty in modeling of the system
         - :literal:`measurementUncertainty`: a matrix defining the uncertainty in modeling of the measurements
         - :literal:`filteringTimeStep`: a scalar representing the value of the constant filtering time step
         - :literal:`initialTime`: a scalar representing the value of the initial time
         - :literal:`initialStateVector`: a vector representing the initial (estimated) state of the system; it is used as first a-priori estimate of the state vector
         - :literal:`initialCovarianceMatrix`: a matrix representing the initial (estimated) covariance of the system; it is used as first a-priori estimate of the covariance matrix
         - :literal:`integratorSettings`: a pointer to integration settings defining the integrator to be used to propagate the state

Unscented Kalman Filter
~~~~~~~~~~~~~~~~~~~~~~~

.. class:: UnscentedKalmanFilter

   .. code-block:: cpp

      boost::shared_ptr< FilterBase< IndependentVariableType, DependentVariableType > > filterObject;
      filterObject = createFilter( filterSettings,
                                   systemFunction,
                                   measurementFunction )

   Here, the inputs are as follows:

      - :literal:`filterSettings`: a pointer to a :class:`UnscentedKalmanFilterSettings` object, described at the end of this section
      - :literal:`systemFunction`: a function returning the state as a function of time and state vector; this can be a differential equation if the :literal:`integratorSettings` is set (i.e., if it is not a :literal:`nullptr`)
      - :literal:`measurementFunction`: a function returning the measurement as a function of time and state

   The filter settings above, can be defined for an extended Kalman filter as:

   .. class:: UnscentedKalmanFilterSettings

      .. code-block:: cpp

         boost::shared_ptr< FilterSettings< IndependentVariableType, DependentVariableType > > filterSettings = 
            boost::make_shared< UnscentedKalmanFilterSettings< IndependentVariableType, DependentVariableType > >( 
                        systemUncertainty,
                        measurementUncertainty,
                        filteringTimeStep,
                        initialTime,
                        initialStateVector,
                        initialCovarianceMatrix,
                        integratorSettings,
                        constantValueReference,
                        customConstantParameters );

      where each element is defined as:

         - :literal:`systemUncertainty`: a matrix defining the uncertainty in modeling of the system
         - :literal:`measurementUncertainty`: a matrix defining the uncertainty in modeling of the measurements
         - :literal:`filteringTimeStep`: a scalar representing the value of the constant filtering time step
         - :literal:`initialTime`: a scalar representing the value of the initial time
         - :literal:`initialStateVector`: a vector representing the initial (estimated) state of the system; it is used as first a-priori estimate of the state vector
         - :literal:`initialCovarianceMatrix`: a matrix representing the initial (estimated) covariance of the system; it is used as first a-priori estimate of the covariance matrix
         - :literal:`integratorSettings`: a pointer to integration settings defining the integrator to be used to propagate the state
         - :literal:`constantValueReference`: reference to be used for the values of the :math:`\alpha` and :math:`\kappa` parameters; this variable has to be part of the :literal:`ConstantParameterReferences` enumeration (custom parameters are supported)
         - :literal:`customConstantParameters`: values of the constant parameters :math:`\alpha` and :math:`\kappa`, input as a pair, in case the :literal:`custom_parameters` enumeration is used in the previous field

      The supported values for :literal:`constantValueReference` are based on UKF literature, and correspond to the values of :math:`\alpha` and :math:`\kappa` shown in the table below (where :math:`N` denotes the length of the state vector). Note that the values of these two parameters have the following meaning: 

         - :math:`\alpha` is used to distribute the sigma points around the a-priori estimate
         - :math:`\kappa` is a secondary scaling parameter, also used to distribute the sigma points around the a-priori estimate, but it has a smaller influence

      ====================================================  ==============  ===============
      :literal:`constantValueReference`                     :math:`\alpha`  :math:`\kappa`
      ====================================================  ==============  ===============
      :literal:`reference_Wan_and_Van_der_Merwe` (default)  0.003           0.0
      :literal:`reference_Lisano_and_Born_and_Axelrad`      1.0             3.0 - :math:`N`
      :literal:`reference_Challa_and_Moore_and_Rogers`      0.001           0.1
      :literal:`customConstantParameters`                   --              --
      ====================================================  ==============  ===============

Using a Kalman Filter Object
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To use the filters after their creation, you will need to update the system with the data corresponding to the new time step. This can be done by using the function :literal:`updateFilter`, which takes one inputs:

   - :literal:`currentMeasurementVector`: a vector denoting the current external measurement which will be used to correct the estimated a-priori state (thus to obtain the a-posteriori estimate)

Note that the filter object comes equipped with two noise generators. These produce random Gaussian noise based on the system and measuremenet uncertainty properties input by you. They can be retrieved with the commands :literal:`produceSystemNoise` and :literal:`produceMeasurementNoise`. You will find cases of how to use the noise and the other features mentioned on this page, in::

   /Tudat/Mathematics/Filters/UnitTests

where a few examples for each filtering technique are shown. One of the test cases shown for the extended (EKF) and unscented (UKF) filters is also shown in an example application, which you can find in :ref:`walkthroughsFiltering`.

The last element that should be discussed is the control system. It may be the case that you also need a control input together with the time and state. This can be added with the creation of a :class:`ControlWrapper` of the type explained in the next section, :ref:`tudatFeaturesFiltersControlSystem`.

These are the only major steps that you will need to take to keep the filter running. At the end of the estimation process, however, you can retrieve, however, retireve the history of the estimated states and covariance matrices by using the functions :literal:`getEstimatedStateHistory` and :literal:`getEstimatedCovarianceHistory`, respectively.

.. _tudatFeaturesFiltersControlSystem:

Adding a Control Input
~~~~~~~~~~~~~~~~~~~~~~

It may be that for your application, your system needs to be controlled. This can be achieved by creating a control system (:class:`ControlSystem`, or :class:`ControlWrapper` in the unit tests), that provides the state function (and possibly the state Jacobian) with the current commanded vector. 

The way this is done is by adding an extra input to the state function, i.e., 

   .. code-block:: cpp

      Eigen::Vector3d stateFunction( const double currentTime, const Eigen::Vector3d& currentState, const Eigen::Vector3d& currentControl )

The last input of the :literal:`stateFunction` is another vector which denotes the commanded control vector. Note how this function, however is incompatible with the definition of :literal:`systemFunction`, i.e., the input to the :class:`ExtendedKalmanFilter` and :class:`UnscentedKalmanFilter` classes. In fact, this should be of type

   .. code-block:: cpp

      boost::function< DependentVector( const IndependentVariableType, const DependentVector& ) >

which one has two inputs (a :literal:`double` and an :literal:`Eigen::Vector3d`, in our case). The extra vector can be added by using the very handy :literal:`boost::bind` command. This function allows us to bind the output of a function as an input to another function. Thus, by using the control class mentioned above, one can replace the :literal:`systemFunction` input with:

   .. code-block:: cpp

      boost::bind( &stateFunction, _1, _2, boost::bind( &ControlSystem::getCurrentStateVector, controlSystem ) )

where the :literal:`controlSystem` object is of type :class:`ControlSystem`. Using this method, if the control system is regularly updated (possibly every time step), the control vector will be automatically retrieved and parsed by the filter.

In :ref:`walkthroughsFiltering` you will find an example where a control system is added to both an EKF and a UKF, which are used to estimate the state of a simple estimation application.

.. tip:: In general, you can apply the principle above to feed any other variable or object to the state and/or measurement functions. 