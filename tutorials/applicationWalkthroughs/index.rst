.. _walkthroughsIndex:

Application Walkthroughs
========================
These pages of the wiki will walk you through several example applications included in the Tudat Bundle. These example make use of several Tudat libraries and show how to use the different elements for your own simulation.

Preparation
~~~~~~~~~~~
The tutorials on these pages focus on the elements of Tudat that you'll need from to write your own simulation. The source files for the examples we use can be found in the ``tudatExampleApplications`` folder within the ``tudatBundle``. If you have downloaded the Tudat Bundle, the following directory should be on your computer::

    ../pathTo/tudatExampleApplications/satellitePropagatorExamples/SatellitePropagatorExamples/

where "pathTo" is the directory where you extracted the Tudat Bundle. When browsing the files in this directory you will see, among others, the ``singleSatellitePropagator.cpp`` and ``singlePerturbedSatellitePropagator.cpp`` files, which we will discuss here among other examples.

Example Tutorials
~~~~~~~~~~~~~~~~~
Each example is discussed in detail on a separate page. Note that in the first example, the entire source file is discussed in a step-by-step manner. Subsequent tutorials focus on the aspects that are different, compared to the basic Kepler orbit propagation.

.. toctree::
   :numbered:
   :maxdepth: 1

   tabulatedAtmosphereExamples
   unperturbedEarthOrbitingSatellite
   perturbedEarthOrbitingSatellite
   propagatorTypesComparison
   unguidedCapsuleEntry
   innerSolarSystemPropagation
   useOfThrustThrustForceAlongVelocityVector
   useOfThrustUserDefinedThrustVector
   variationalEquationsPropagation
   earthOrbiterBasicStateEstimation
   earthOrbiterStateEstimation
   interplanetaryTrajectory
   
Optimization Tutorials
~~~~~~~~~~~~~~~~~~~~~~
In the same tudatExampleApplications library on your computer, you can find the library examples. These examples contain examples on how to use Pagmo 2, which is an external optimization library. The first tutorial, Himmelblau optimization, will show you the basics of this library, the subsequent examples go into more depth of the Pagmo 2 library.

.. toctree::
   :numbered:
   :maxdepth: 1

   himmelblau
   cecOptimization
   earthMarsTrans
   mgaTrans
   propTargeting
