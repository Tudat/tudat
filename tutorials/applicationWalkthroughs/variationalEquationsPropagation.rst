.. _variationalEquationPropagation:

Variational Equations Propagation
=================================

In the previous tutorials, we have focussed on propagating the (translational and mass) state of one of or more bodies. Tudat is also capable of propagating the variational equations for the dynamics, to produce the state transition matrix :math:`\Phi(t,t_{0})` and sensitivity matrix :math:`S(t)`, which we define here as:

.. math::
      
      \Phi(t,t_{0}) &= \frac{\partial \mathbf{x}(t)}{\partial\mathbf{x}(t_{0})}\\
      S &= \frac{\partial \mathbf{x}(t)}{\partial \mathbf{p  }}\\

where :math:`\mathbf{x}` is the propagated state, :math:`\mathbf{p}` the vector of a parameter vector (e.g. gravity field parameters, rotation model parameters, etc.), and :math:`t_{0}` denotes the initial time.

.. note:: In some literature, the sensitivity matrix is not defined separately, but the state transition matrix :math:`\Phi(t,t_{0})` is defined as :math:`\frac{\partial[\mathbf{x}(t);\text{ }\mathbf{p}]}{\partial[\mathbf{x}(t_{0};\text{ }\mathbf{p}])}`

The code for this tutorial is given on Github, and is also located in your Tudat bundle at::

    tudatBundle/tudatExampleApplications/satellitePropagatorExamples/SatellitePropagatorExamples/variationalEquationsPropagatorExample.cpp

As you can see in the code, setting up the simulations is done in a similar manner as the :ref:`walkthroughsPerturbedEarthOrbitingSatellite`: the code start by providing settings for the environment, accelerations, propagation and integration.

To propagate the variational equations, additional information needs to be provided, specifically: for which paramaters are the variational equations to be set up? In this example, the following code fully defines this:

.. code-block:: cpp

   // Define list of parameters to estimate.
   std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;
   parameterNames.push_back( std::make_shared< InitialTranslationalStateEstimatableParameterSettings< double > >(
                                 "Asterix", asterixInitialState, "Earth" ) );
   parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Asterix", radiation_pressure_coefficient ) );
   parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Asterix", constant_drag_coefficient ) );
   parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                                 2, 0, 2, 2, "Earth", spherical_harmonics_cosine_coefficient_block ) );
   parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                                 2, 1, 2, 2, "Earth", spherical_harmonics_sine_coefficient_block ) );

   // Create parameters
   std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate =
           createParametersToEstimate( parameterNames, bodyMap );

   // Print identifiers and indices of parameters to terminal.
   printEstimatableParameterEntries( parametersToEstimate );

Details on the setup, and the options available, for the parameters is given on the page: :ref:`parameterEstimationSettings`. Summarizing, the lines of code above define that the following parameters should be used:

   - Initial state of body Asterix, denoted as :math:`\mathbf{x}_{0}`
   - Radiation pressure coefficient (constant) of body Asterix, denoted as :math:`C_{r}`
   - Drag coefficient (constant) of body Asterix, denoted as :math:`C_{D}`
   - Both cosine and sine spherical harmonic gravity field coefficients of body Earth, at degree 2: :math:`\bar{C}_{2,0..2}`, denoted as :math:`\bar{S}_{2,1..2}`

Similar to other models, we create objects of (derived classes of) :class:`EstimatableParameterSettings`, and then call the :literal:`createParametersToEstimate` function to create the actual objects that are used internally in Tudat, in this case stored in an :class:`EstimatableParameterSet` object. 

The propagation of the variational equations solver is done by a dedicated object:

.. code-block:: cpp

   // Create simulation object and propagate dynamics.
   SingleArcVariationalEquationsSolver< > variationalEquationsSimulator(
               bodyMap, integratorSettings, propagatorSettings, parametersToEstimate, true,
               std::shared_ptr< numerical_integrators::IntegratorSettings< double > >( ), false, true );
               
So, instead of using a :class:`SingleArcDynamicsSimulator` object (which propagates only the state), we use :class:`SingleArcVariationalEquationsSolver` object (which propagates both the state and the variational equations). Details on the input to the object are discussed on the page :ref:`tudatFeaturesVariationalSimulatorCreation`. Using the above, the equations are immediately propagated upon creation of the object. Retrieving the output is done as:

.. code-block:: cpp

   std::map< double, Eigen::MatrixXd > stateTransitionResult =
           variationalEquationsSimulator.getNumericalVariationalEquationsSolution( ).at( 0 );
   std::map< double, Eigen::MatrixXd > sensitivityResult =
           variationalEquationsSimulator.getNumericalVariationalEquationsSolution( ).at( 1 );
   std::map< double, Eigen::VectorXd > integrationResult =
           variationalEquationsSimulator.getDynamicsSimulator( )->getEquationsOfMotionNumericalSolution( );

.. warning:: Propagation of variational equations (for translational motion) is only supported for propagated coordinates of size 6. For a description of the difference between conventional and propagated coordinates, see :ref:`tudatFeaturesPropagatorSettingsCoordinates`. On the same page (in :ref:`tudatFeaturesPropagatorSettingsPropagatedCoordinates`) you can find a list of the available propagated coordinates and their respective sizes.

Finally, these maps are written to files, similarly to the previous examples and discussed in :ref:`tudatFeaturesInputOutput`. Note however, that the matrix entries of the first two maps in the above are spread out over a single row in the output file. The concatenation of the matrix entries is done row by row.

Below, a plot is given of the entries of the state transition matrix as a function of time. The current state entry is indicated by line style, the initial state entry by color.

.. figure:: images/variationalEquationsExample.png

.. tip:: Open the figure in a new tab for more detail.
