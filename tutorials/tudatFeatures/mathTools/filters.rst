.. _tudatFeaturesFilters:

Filters
=======
Filters, in astrodynamics, are generally used to estimate the state of a system, based on a model description and measurement data. Tudat provides a few filters, based on the Kalman filter concept, that you can use for these or any other purpose you see fit. The filtering techinques currently offered in Tudat are:

   - **Linear Kalman Filter**: insert basic description (LKF)
   - **Extended Kalman Filter**: insert basic description (EKF)
   - **Unscented Kalman Filter**: insert basic description (UKF)

You can find more information on the mathematical background of these filters in virtually any guidance, navigation and control book. However, the references used to define the filters in Tudat are:

   - **LKF**: E. Mooij, AE4870B - Re-entry Systems, Lecture Notes, Delft University of Technology, 2016
   - **EKF**: Ogata, K., Discrete-Time Control Systems, 2nd ed. Pearson Education Asia, 2002.
   - **UKF**: Wan, E. and Van Der Merwe, R., “The Unscented Kalman Filter for Nonlinear Estimation,” in Adaptive Systems for Signal Processing, Communications, and Control Symposium. Institute of Electrical and Electronics Engineers, 2000, pp. 153–158.

On this page, you will find a description of how to create a Kalman filter object, either based on the linear, extended, or unscented principle. In general, you will notice that the creation of a filter, requires the input of two template arguments. These are called :literal:`IndependentVariableType` and :literal:`DependentVariableType`, and they respectively denote the type (i.e., :literal:`double`, :literal:`int`, etc.) for the independent variable (usually time) and the dependent variables. Note that the dependent variable is internally defined as an :literal:`Eigen::Matrix< DependentVariableType, Eigen::Dynamic, 1 >` object (which corresponds to a vector of unspecified length).

Linear Kalman Filter
~~~~~~~~~~~~~~~~~~~~

.. class:: LinearKalmanFilter

   .. code-block:: cpp

      boost::shared_pointer< LinearKalmanFilter< IndependentVariableType, DependentVariableType > > linearFilter = 
         boost::make_shared< LinearKalmanFilter< IndependentVariableType, DependentVariableType > >(
                     stateTransitionMatrix,
                     controlMatrix,
                     measurementMatrix,
                     systemUncertainty, measurementUncertainty,
                     initialTime, initialEstimatedStateVector, initialEstimatedStateCovarianceMatrix );

   In the code above, the inputs are given by:

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
                                   measurementNoiseJacobianFunction )

   Here, the inputs are as follows:

      - :literal:`filterSettings`: a pointer to a :class:`ExtendedKalmanFilterSettings` object, described at the end of this section.
      - :literal:`systemFunction`:
      - :literal:`measurementFunction`:
      - :literal:`stateJacobianFunction`:
      - :literal:`stateNoiseJacobianFunction`:
      - :literal:`measurementJacobianFunction`:
      - :literal:`measurementNoiseJacobianFunction`:

   The filter settings above, can be defined for an extended Kalman filter as:

   .. class:: ExtendedKalmanFilterSettings

   .. code-block:: cpp

      boost::shared_ptr< ExtendedKalmanFilterSettings< IndependentVariableType, DependentVariableType > > filterSettings = 
         boost::make_shared< ExtendedKalmanFilterSettings< IndependentVariableType, DependentVariableType > >( 
                     systemUncertainty,
                     measurementUncertainty,
                     initialTime,
                     initialStateVector,
                     initialCovarianceMatrix,
                     integratorSettings );

   where each element is defined as:

      - :literal:`systemUncertainty`: 
      - :literal:`measurementUncertainty`: 
      - :literal:`initialTime`: 
      - :literal:`initialStateVector`: 
      - :literal:`initialCovarianceMatrix`: 
      - :literal:`integratorSettings`: 

Unscented Kalman Filter
~~~~~~~~~~~~~~~~~~~~~~~

.. class:: UnscentedKalmanFilter

   .. code-block:: cpp

      boost::shared_ptr< FilterBase< IndependentVariableType, DependentVariableType > > filterObject;
      filterObject = createFilter( filterSettings,
                                   systemFunction,
                                   measurementFunction )

   Here, the inputs are as follows:

      - :literal:`filterSettings`: a pointer to a :class:`UnscentedKalmanFilterSettings` object, described at the end of this section.
      - :literal:`systemFunction`: the system function, which can be passed as a :literal:`boost::function` or as a 
      - :literal:`measurementFunction`:

   The filter settings above, can be defined for an extended Kalman filter as:

   .. class:: UnscentedKalmanFilterSettings

   .. code-block:: cpp

      boost::shared_ptr< UnscentedKalmanFilterSettings< IndependentVariableType, DependentVariableType > > filterSettings = 
         boost::make_shared< UnscentedKalmanFilterSettings< IndependentVariableType, DependentVariableType > >( 
                     systemUncertainty,
                     measurementUncertainty,
                     initialTime,
                     initialStateVector,
                     initialCovarianceMatrix,
                     integratorSettings,
                     constantValueReference,
                     customConstantParameters );

   where each element is defined as:

      - :literal:`systemUncertainty`: 
      - :literal:`measurementUncertainty`: 
      - :literal:`initialTime`: 
      - :literal:`initialStateVector`: 
      - :literal:`initialCovarianceMatrix`: 
      - :literal:`integratorSettings`: 
      - :literal:`constantValueReference`: 
      - :literal:`customConstantParameters`: 

Control System
~~~~~~~~~~~~~~

