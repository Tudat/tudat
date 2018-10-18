.. _tudatFeaturesVariationalSimulatorCreation:

Performing Variational Equations Propagation
============================================

Up to this point, we have been concerned with propagating states of bodies only. Tudat is also capable of propagating the variational equations associated with the dynamics to produce the state transition matrix :math:`\Phi(t,t_{0})` and sensitivity matrix :math:`S(t)`, which we define here as:

.. math::
      
      \Phi(t,t_{0}) &= \frac{\partial \mathbf{x}(t)}{\partial\mathbf{x}(t_{0})}\\
      S &= \frac{\partial \mathbf{x}(t)}{\partial \mathbf{p  }}\\

where :math:`\mathbf{x}` is the propagated state, :math:`\mathbf{p}` the vector of a parameter vector (e.g. gravity field parameters, rotation model parameters, etc.), and :math:`t_{0}` denotes the initial time.

.. note:: In some literature, the sensitivity matrix is not defined separately, but the state transition matrix :math:`\Phi(t,t_{0})` is defined as :math:`\frac{\partial[\mathbf{x}(t);\text{ }\mathbf{p}]}{\partial[\mathbf{x}(t_{0};\text{ }\mathbf{p}])}`

The propagation of these equations is done similarly as for the dynamics, and is executed by a (derived class of) :class:`VariationalEquationsSolved`

At the moment, the following :class:`VariationalEquationsSolver` options are available or under development in Tudat:

- Single-arc variational equations solver.
- Multi-arc variational equations solver.
- Hybrid variational equations solver (under development).

These are implemented in derived classes and are discussed below. Note that these objects of these classes propagate bothe the variational equations and dynamics (either concurrently or sequentially). 

.. class:: VariationalEquationsSolver

   Base class from which the classes below are derived.

.. class:: SingleArcVariationalEquationsSolver
   
   This derived class simulates single arc dynamics and variational equations and its constructor is:

   .. code-block:: cpp

      SingleArcVariationalEquationsSolver( bodyMap, 
                                           integratorSettings, 
                                           propagatorSettings, 
                                           parametersToEstimate, 
                                           integrateDynamicalAndVariationalEquationsConcurrently,
                                           variationalOnlyIntegratorSettings, 
                                           clearNumericalSolution,
                                           integrateEquationsOnCreation );

   where:

   - :literal:`bodyMap`

      :class:`NamedBodyMap` the map containing all the objects of type :class:`Body` used in the simulation.

   - :literal:`integratorSettings`

      :class:`IntegratorSettings` contains the settings of the integrator used, as discussed in :ref:`tudatFeaturesIntegratorSettings`. 

   - :literal:`propagatorSettings`

      :class:`PropagatorSettings` contains the settings that defines how the dynamics is propagated, as described in :ref:`tudatFeaturesPropagatorSettings`.

      .. warning:: Propagation of variational equations (for translational motion) is only supported for propagated coordinates of size 6. For a description of the difference between conventional and propagated coordinates, see :ref:`tudatFeaturesPropagatorSettingsCoordinates`. On the same page (in :ref:`tudatFeaturesPropagatorSettingsPropagatedCoordinates`) you can find a list of the available propagated coordinates and their respective sizes.
      
   - :literal:`parametersToEstimate`

      :class:`EstimatableParameterSet` contains the settings that define the parameters :math:`\mathbf{x}_{0}` and :math:`\mathbf{p}` for which the variational equations are to be solved, as described in :ref:`parameterEstimationSettings`.
      
   - :literal:`integrateDynamicalAndVariationalEquationsConcurrently`

      Boolean to denote whether the equations of motion and variational equations are to be propagated concurrently (default: true), or if the variational eqautions are to be solved after the equations of motion.
      
   - :literal:`variationalOnlyIntegratorSettings`

     :class:`IntegratorSettings` contains the settings of the integrator used, as discussed in :ref:`tudatFeaturesIntegratorSettings`, when propagating the variational equations separately (if :literal:`integrateDynamicalAndVariationalEquationsConcurrently` is false). This pointer is empty by default, in which case the :literal:`integratorSettings` are used.
      
   - :literal:`clearNumericalSolution`

      Boolean to denote whether numerical solutions of the propagated equations can be retrieved manually from the object after propagation (if false), or if they are cleared from memory (if true).
      
   - :literal:`integrateEquationsOnCreation`

      Boolean to denote whether the equations of motion and variational equations are to be propagated immediately when the object is created (default true).

.. class:: MultiArcVariationalEquationsSolver
   
   This derived class allows the numerical propagation of variational equations for arc-wise dynamics. It is constructed using:

   .. code-block:: cpp
   
      MultiArcVariationalEquationsSolver( bodyMap, 
                                          integratorSettings, 
                                          propagatorSettings, 
                                          parametersToEstimate, 
                                          arcStartTimes,
                                          integrateDynamicalAndVariationalEquationsConcurrently,
                                          variationalOnlyIntegratorSettings, 
                                          clearNumericalSolution,
                                          integrateEquationsOnCreation ) );

   where:

   - :literal:`arcStartTimes`

      :literal:`std::vector< double >` contains the times at which the separate arcs start.

.. class:: HybridArcVariationalEquationsSolver

   Allows some bodies to be propagated in a single arc, and some in a multi-arc fashion. See :class:`HybridArcDynamicsSimulator` for model details. It is constructed using:

   .. code-block:: cpp
   
      HybridArcVariationalEquationsSolver( bodyMap, 
                                          integratorSettings, 
                                          propagatorSettings, 
                                          parametersToEstimate, 
                                          arcStartTimes,
                                          integrateDynamicalAndVariationalEquationsConcurrently,
                                          clearNumericalSolution,
                                          integrateEquationsOnCreation ) );

      

Retrieving the variational equation history
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Once the :class:`VariationalEquationsSolver` object has been created and the equations of motion have been integrated, the propagation history of the selected bodies is stored within the :class:`VariationalEquationsSolver`. To make use of it manually after propagation, such history can to be retrieved and saved to a file.

If the variational equations propagation history needs to be saved, the following code can be used (assuming an object of type :class:`SingleArcVariationalEquationsSolver` called variationalEquationsSimulator has been created):

.. code-block:: cpp

    std::map< double, Eigen::MatrixXd > stateTransitionResult =
            variationalEquationsSimulator.getNumericalVariationalEquationsSolution( ).at( 0 );
    std::map< double, Eigen::MatrixXd > sensitivityResult =
            variationalEquationsSimulator.getNumericalVariationalEquationsSolution( ).at( 1 );
          
Where the state transition and sensitivity matrix are storeed separately. Saving the results to a file is done in the same manner as for the dynamics. Note however, that the matrix entries of the maps in the above are spread out over a single row in the output file. The concatenation of the matrix entries is done row by row. 

The state propagation history can also be retrieved from teh object, as follows:

.. code-block:: cpp

       std::map< double, Eigen::VectorXd > integrationResult =
            variationalEquationsSimulator.getDynamicsSimulator( )->getEquationsOfMotionNumericalSolution( );




