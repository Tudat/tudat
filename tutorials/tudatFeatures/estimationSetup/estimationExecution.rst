.. _estimationExecution:

Parameter Estimation 
=========================

.. _estimationObjectCreation:

Creating the estimation object
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The estimation of orbits and physical parameters in Tudat is done through the use of an :class:`OrbitDeterminationManager` object. This object is created once for a given set of estimated parameters, and may be (re-)used to reset and/or estimate parameters from a given set of data.

.. class:: OrbitDeterminationManager

Creating an object of this type automatically propagates both the equations of motion and variational equations using the input to its constructor. Additionally, a set of :class:`ObservationSimulatorBase` objects is created, as well as all required objects to compute the partial derivatives of the observations.

An object is created as follows:

   .. code-block:: cpp

      OrbitDeterminationManager orbitDeterminationManager =
            std::make_shared< ObservationViabilitySettings >( 
                OrbitDeterminationManager( bodyMap, parametersToEstimate, observationSettingsMap, integratorSettings, propagatorSettings );
                
   The input is:

   - :literal:`bodyMap`

      A :literal:`NamedBodyMap` variable which contains all information on the physical environment
      
   - :literal:`parametersToEstimate`

      A :literal:`std::shared_ptr< EstimatableParameterSet< ObservationScalarType > >` variable which contains the full list of parameters that are to be estimated.
      
   - :literal:`observationSettingsMap`

      A list of :class:`ObservationSettings` that may be provided as either :literal:`std::multimap< LinkEnds, std::shared_ptr< ObservationSettings > >` or a :literal:`std::map< ObservableType, std::map< LinkEnds, std::shared_ptr< ObservationSettings > >`. In the former case, the list is sorter by :literal:`LinkEnds` only, in the latter by both :literal:`ObservableType` and :literal:`LinkEnds`
      
   - :literal:`integratorSettings`

      A :literal:`std::shared_ptr< IntegratorSettings< TimeType > >` object which contains all settings for the numerical integrator, used to numerically calcualte the dynamics and variational equations. 

   - :literal:`propagatorSettings`

      A :literal:`std::shared_ptr< PropagatorSettings< ObservationScalarType > >` object which contains all propagation settings. 
      
      .. warning::
      
         The settings in :class:`PropagatorSettings` and :class:`EstimatableParameterSet` must be consistently defined: any initial state that is propagated *must* also be estimated, there is as yet no possibility to only estimate some of the propagates states.
         
After the creation of the :class:`OrbitDeterminationManager` object, a number of objects are created internally that can be used for various purposes:

* An object with base class :class:`VariationalEquationsSolver`, which contains the numerically propagated variational equations and dynamics, can be retrieved using the :literal:`getVariationalEquationsSolver` member function. 
* A list of objects with base class :class:`ObservationSimulatorBase` (one per observable type) to simulate observations, discussed in more detail on the page on :ref:`creatingObservationSimulators`
* A list of objects with base class :class:`ObservationManagerBase` (one per observable type) to simulate observations and the associated partial derivatives. These objects are not directly accesed by users. Their output (partial derivatives of observables) are provided *a posterior* through an object of type :class:`PodOutput`, discussed on the page on :ref:`estimationOutput`.

.. _estimationInput:

Defining estimation input
~~~~~~~~~~~~~~~~~~~~~~~~~

The input to the estimation consists of several parts. Firstly, the input data, weights, *etc.* need to be defined, which is done through the :literal:`PodInput` class. 

.. class:: PodInput

   This class is templated by both :literal:`ObservationScalarType` and :literal:`TimeType`. An object of :class:`PodInput` is created as follows:
   
   .. code-block:: cpp

      std::shared_ptr< PodInput< ObservationScalarType, TimeType > > podInput =
            std::make_shared< PodInput< ObservationScalarType, TimeType >  >( 
                observationsAndTimes, numberOfEstimatedParameters, inverseOfAprioriCovariance );
     
   The input is:
   
   - :literal:`observationsAndTimes` A container of type :literal:`std::map< ObservableType, std::map< LinkEnds, std::pair< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >, std::pair< std::vector< TimeType >, LinkEndType > > > >` (the structure of which is described in more detail on the page :ref:`generatingObservations`). This container has both the observables to be used in the estimation, and the assictaed times and link end types.
   
   - :literal:`numberOfEstimatedParameters` An :literal:`int` denoting the length of the vector of estimated paramaters, discussed in more detail on the page :ref:`parameterSettingCreation`.
   
   - :literal:`inverseOfAprioriCovariance` An :literal:`Eigen::MatrixXd` with the inverse of the *a priori* covariance matrix. This input type may be left empty, in which case no *a priori* covariance is used.


.. note::
   
   Currently, Tudat only supports diagonal weight matrices, implicitly assuming independent observation noise in the inversion.

.. _estimationOutput:

Estimation output
~~~~~~~~~~~~~~~~~


When performing the estimation, the code rescales the values of all parameters :math:`p`, where we denote the scaled parameters as :math:`\tilde{h}`, so that all partials :math:`\partial h/\partial\tilde{p}` w.r.t. lie in the range :math:`[-1,1]`. To provide transparency, it is the covariance and partial derivative matrix of these scaled parameters that is saved to the :literal:`PodOutput` object. However, the following functions allow you to retrieve the information w.r.t. the *unscaled* parameters:

* Inverse covariance, obtained using the :literal:`getUnnormalizedInverseCovarianceMatrix` function.
* Covariance, obtained using the :literal:`getUnnormalizedCovarianceMatrix` function. Note that this only produces valid results if the problem is not ill-posed.
* Formal error vector, obtained using the :literal:`getFormalErrorVector` function. Note that this only produces valid results if the problem is not ill-posed.
* Correlation matrix, obtained using the :literal:`getCorrelationMatrix` function. Note that this only produces valid results if the problem is not ill-posed.



