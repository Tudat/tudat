.. _walkthroughsFullPropagationPatchedConicsTrajectory:

Full propagation of a patched conics trajectory
===============================================

This example addresses the difference between the trajectory of a spacecraft obtained from the patched conics solution and the one obtained after propagating the full dynamics problem. The code for this tutorial is available on GitHub, and is also located in the Tudat bundle at: ::

    tudatBundle/tudatExampleApplications/satellitePropagatorExamples/SatellitePropagatorExamples/fullPropagationPatchedConicsTrajectory.cpp

The patched conics trajectory is based on the assumption that the spacecraft to be propagated is only affected by the point-mass gravitational acceleration exerted by the central body of the trajectory. This simplifying assumption is however never true in real-world applications and the patched conics trajectory thus differs from the results provided by the propagation of the full dynamics problem.

This example is based on the Messenger trajectory example already described in the :ref:`interplanetaryTrajectoryDesign`. For further details about the way this trajectory is calculated, the reader is referred to this tutorial.

Ideal case corresponding to the patched conics trajectory assumptions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The first part of this example focuses on the case in which the full dynamics problem actually corresponds to the assumptions of the patched conics approach.

The type and number of legs are defined in the exact same way as it is done in the interplanetary design example, and so are the characteristics of the trajectory (trajectory variables vector, minimum periapasis radii vector, semi-major axes and eccentricities at departure and arrival. Their definition will thus not be further described here.

Then, the name of the body to be propagated and of the central body of the trajectory are defined, as well as the names of all the transfer bodies involved in the trajectory.

.. code-block:: cpp
   
   // Define names of the central body and body to be propagated.
   std::vector< std::string > centralBody; centralBody.push_back( "Sun" );
   std::string bodyToPropagate = "spacecraft";

   // Define names of the bodies involved in the trajectory.
   std::vector< std::string > transferBodyTrajectory;
   transferBodyTrajectory.push_back("Earth");
   transferBodyTrajectory.push_back("Earth");
   transferBodyTrajectory.push_back("Venus");
   transferBodyTrajectory.push_back("Venus");
   transferBodyTrajectory.push_back("Mercury");

The environment is which the trajectory is to be calculated has then to be defined. First, the body map corresponding to the trajectory problem has to be created. The function :literal:`setupBodyMapFromEphemeridesForPatchedConicsTrajectory` directly sets up a body map which fulfills the assumptions of the patched conics approach. In a similar way, the accelerations acting on the spacecraft along the trajectory are defined with the function :literal:`setupAccelerationMapPatchedConicsTrajectory` which automatically returns an :literal:`accelerationMap` object corresponding to the patched conics assumptions.

.. code-block:: cpp

   // Create body map.
   simulation_setup::NamedBodyMap bodyMap = propagators::setupBodyMapFromUserDefinedEphemeridesForPatchedConicsTrajectory(centralBody[0],
            bodyToPropagate, nameBodiesTrajectory, ephemerisVectorTransferBodies, gravitationalParametersTransferBodies);

   // Create acceleration map.
   std::vector< basic_astrodynamics::AccelerationMap > accelerationMap = propagators::setupAccelerationMapPatchedConicsTrajectory(
                nameBodiesTrajectory.size(), centralBody[0], bodyToPropagate, bodyMap);

The integrator settings remain to be chosen before the propagation of the spacecraft trajectory can be performed:

.. code-block:: cpp

   // Define integrator settings.
   double fixedStepSize = 1000.0;
   std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings =
       std::make_shared < numerical_integrators::IntegratorSettings < > > ( numerical_integrators::rungeKutta4, initialTime, fixedStepSize);

Finally, the patched conics solution and the propagation results of the full dynamics problem are provided jointly by the function :literal:`fullPropagationPatchedConicsTrajectory`. This function does not return any variable but fills the maps lambertTargeterResultForEachLeg and fullProblemResultForEachLeg (provided as inputs) with the results of the patched conics method and the full problem propagation, respectively. These two maps are constructed in the same way: the first element is the number of the leg to be considered and the second element is the state history of the spacecraft along this leg, either obtained from the patched conics solution or the full problem propagation depending on the map.

.. code-block:: cpp

  // Calculate the patched conics solution and the propagation results of the associated full dynamics problem.


The difference in cartesian state between the full problem and the patched conics approach at departure and arrival of each leg can be retrieved by calling the function :literal:`getDifferenceFullProblemWrtPatchedConicsTrajectory` which returns a map whose first element is the number of the leg and the second element is a pair of vector containing the state difference for this leg at departure and arrival respectively.

.. code-block:: cpp

   // Compute difference between patched conics trajectory and full problem at departure and at arrival for each leg.
   std::map< int, std::map< double, Eigen::Vector6d > > lambertTargeterResultForEachLeg;
   std::map< int, std::map< double, Eigen::Vector6d > > fullProblemResultForEachLeg;

   std::map< int, std::pair< Eigen::Vector6d, Eigen::Vector6d > > differenceStateArrivalAndDeparturePerLeg =
            propagators::getDifferenceFullProblemWrtPatchedConicsTrajectory(bodyMap, accelerationMap, nameBodiesTrajectory,
                            centralBody[0], bodyToPropagate, legTypeVector, variableVector, minimumPericenterRadii,
                            semiMajorAxes, eccentricities, integratorSettings);


Because the dynamical model used to define the full problem has been chosen to exactly correspond to the assumptions of the patched conics method, the differences in state along the trajectory are expected to be extremely low.


Perturbed case
~~~~~~~~~~~~~~

The first part of this tutorial has focused on the ideal case in which the dynamical model used to define the full problem corresponds to the assumptions made in the patched conics approach. This is why no significant differences were observed between the patched conics solution and the results of the full problem propagation. However, these simplifying assumptions might not be realistic for real-world applications and more complete models can be applied to adress the influence of more complex dynamics on a spacecraft trajectory. 

The characteristics of the trajectory and the initial state of the spacecraft do not have to be modified compared to the previous examples. However, the body settings can be manually defined to account for more complex models.

.. code-block:: cpp

   // Create body map.

Similarly, the set of accelerations acting on the spacecraft now include additional perturbations and is no longer restricted to point-mass gravitational acceleration exerted by the central body.

.. code-block:: cpp

   // Define accelerations acting on the spacecraft.

Again, the function :literal:`fullPropagationPatchedConicsTrajectory` is used to retrieve the state history of the spacecraft for each leg, from both the patched conics solution and the full problem propagation.

.. code-block:: cpp

   // Calculate the patched conics trajectory and propagate the full dynamics problem jointly.


When a more complex dynamical model is applied, the differences in cartesian state between the patched conics trajectory and the full problem results become significantly.

