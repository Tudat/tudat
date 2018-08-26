.. _tudatFeaturesSimulatorIndex:

Simulator Set-Up
================

One of the core elements of the Tudat libraries is its simulator framework. The goal of this page is discuss the implementation of such framework as well as the numerous options available. The top-level framework of the simulator is shown below:


.. graphviz::

   digraph
   {
      # General diagram settings
      rankdir = "LR";
      #ranksep = "1.0"
      splines = ortho;    
      compound = true;  


      # general node settings 
      node [shape = box, style = filled, width = 1.25, fixedsize = true, color = lightgrey, fontname = FontAwesome, fontsize = 9];


      # specific node color settings
      DynamicsSimulator [color = lightgreen];
      NamedBodyMap, MultiArc, "RungeKuttaVariable\nStepSizeSettings", IntegratorSettings [color = lightblue];
      "DependentVariable(s)\nNumerical Solution", "Equations of Motion\nNumerical Solution" [color = peachpuff];
      
      # Hyperlinks (Sphinx auto referencing not working here, need to link to exact web adres)
      NamedBodyMap [href = "http://tudat.tudelft.nl/tutorials/tudatFeatures/environmentSetup/index.html#NamedBodyMap", target = "_top"];
      DynamicsSimulator [href = "http://tudat.tudelft.nl/tutorials/tudatFeatures/propagationSetup/simulatorCreation.html#DynamicsSimulator", target = "_top"];
      IntegratorSettings [href = "http://tudat.tudelft.nl/tutorials/tudatFeatures/propagationSetup/integratorSettings.html#IntegratorSettings", target = "_top"];
      "RungeKuttaVariable\nStepSizeSettings" [href = "http://tudat.tudelft.nl/tutorials/tudatFeatures/propagationSetup/integratorSettings.html#RungeKuttaVariableStepSizeSettings", target = "_top"];
      TranslationalState [href = "http://tudat.tudelft.nl/tutorials/tudatFeatures/propagationSetup/propagatorSettings.html#TranslationalStatePropagatorSettings", target = "_top"];
      RotationalState [href = "http://tudat.tudelft.nl/tutorials/tudatFeatures/propagationSetup/propagatorSettings.html#RotationalStatePropagatorSettings", target = "_top"];
      Mass [href = "http://tudat.tudelft.nl/tutorials/tudatFeatures/propagationSetup/propagatorSettings.html#MassPropagatorSettings", target = "_top"];
      Custom [href = "http://tudat.tudelft.nl/tutorials/tudatFeatures/propagationSetup/propagatorSettings.html#CustomPropagatorSettings", target = "_top"];
      MultiType [href = "http://tudat.tudelft.nl/tutorials/tudatFeatures/propagationSetup/propagatorSettings.html#MultiTypePropagatorSettings", target = "_top"];      
      MultiArc [href = "http://tudat.tudelft.nl/tutorials/tudatFeatures/propagationSetup/propagatorSettings.html#MultiArcPropagatorSettings", target = "_top"];      

      subgraph clusterIntegratorSettings
      {
         # cluster settings
         label = "IntegratorSettings";
         fontsize = 9;
         style = dashed;
         href = "http://tudat.tudelft.nl/tutorials/tudatFeatures/propagationSetup/integratorSettings.html";
         target = "_top";

 
         # IntegratorSettings input
         integratorType -> IntegratorSettings;
         simulationStartEpoch -> IntegratorSettings;
         "(initial) stepSize" -> IntegratorSettings;


         # RungeKuttaVariableStepSizeSettings input
         integratorType -> "RungeKuttaVariable\nStepSizeSettings";
         simulationStartEpoch -> "RungeKuttaVariable\nStepSizeSettings";
         "(initial) stepSize" -> "RungeKuttaVariable\nStepSizeSettings";
         "Additional \nsettings" -> "RungeKuttaVariable\nStepSizeSettings";
      }


      subgraph clusterSingleArcPropagatorSettings
      {
         # cluster settings
         label = "Single-arc PropagatorSettings";
         fontsize = 9;
         style = dashed;
         compound = true;
         href = "http://tudat.tudelft.nl/tutorials/tudatFeatures/propagationSetup/propagatorSettings.html";
         target = "_top";

         MultiType;
         TranslationalState -> RotationalState [style = invis];
         Mass -> Custom [style = invis];
      }

     
      # DynamicsSimulator input, use points for "or" blocks
      NamedBodyMap -> DynamicsSimulator;
      or0, or1 [shape = circle, width = 0.2, fillcolor = white, color = black, label = or];
      Custom -> or1 -> DynamicsSimulator [ltail = clusterSingleArcPropagatorSettings];
      MultiArc -> or1 [constraint = false];
      IntegratorSettings -> or0 -> DynamicsSimulator;
      "RungeKuttaVariable\nStepSizeSettings" -> or0;

      # DynamicsSimulator output
      DynamicsSimulator -> "Equations of Motion\nNumerical Solution";
      DynamicsSimulator -> "DependentVariable(s)\nNumerical Solution";
      {rank = same; "Equations of Motion\nNumerical Solution", "DependentVariable(s)\nNumerical Solution"};


      # Extra lines for posisitioning of nodes
      {rank = same; NamedBodyMap, DynamicsSimulator}    
      MultiType -> MultiArc [ltail = clusterSingleArcPropagatorSettings];
      Mass -> Custom [style = invis];
   }

.. graphviz::

   digraph
   {
      # General diagram settings
      rankdir = "LR";
      splines = ortho;    
      compound = true;  

      subgraph clusterLegend
      {
      rank = min;
      style = dashed;


     	# general node settings 
     	node [shape = box, style = filled, width = 1.25, fixedsize = true, color = lightgrey, fontname = FontAwesome, fontsize = 9];


   	"Main block" [fillcolor = lightgreen];
     	"Optional input" [style = dotted, fillcolor = lightgrey, color = black];
     	"Input for \nmain block" [fillcolor = lightblue];
      Output [color = peachpuff];
     	"Optional input"-> "Input for \nmain block" -> "Main block" -> "Output" [style = invis];
      }
   }

The top-element in such framework is the :class:`DynamicsSimulator`, which is in charge of propagating the equations of motion using the environment and acceleration models discussed in :ref:`tudatFeaturesEnvironmentIndex` and :ref:`tudatFeaturesAccelerationIndex`, respectively. The orbit propation is done according to the specified :class:`IntegratorSettings` and the :class:`PropagatorSettings`, which are discussed in detail in :ref:`tudatFeaturesIntegratorSettings` and :ref:`tudatFeaturesPropagatorSettings`.

As shown in the figure above, there are various types of :class:`IntegratorSettings` and :class:`PropagatorSettings` depending on the particularities of the application at hand. The different types are implemented by means of derived classes, as discussed below:

**Integrator Settings**
    - :class:`IntegratorSettings`
        This is the default derived class which is used whenever a fixed-step size integrator is used.
    - :class:`RungeKuttaVariableStepSizeSettings`
        This is a special derived class which takes additional arguments and is used whenever a Runge-Kutta variable step-size integrator is used.

**Propagator Settings**
    - :class:`TranslationalStatePropagatorSettings`
        This derived class defines the settings to propagate the translational dynamics.
    - :class:`RotationalStatePropagatorSettings`
        This derived class defines the settings to propagate the rotational dynamics.
    - :class:`MassPropagationSettings`
        This derived class defines the settings to propagate the mass of a body.
    - :class:`CustomStatePropagatorSettings`
        This derived class allows the propagation of user-defined dynamics.
    - :class:`MultiTypePropagatorSettings`
        This derived class allows to propagate simultaneously the various propagator classes defined above.
    - :class:`MultiArcPropagatorSettings`
        This derived class allows to propagate the dynamics for multiple arcs.

The reader is referred to the following sections to examine in detail how to create the dynamics simulator:

.. toctree::

   integratorSettings
   propagatorSettings
   propagatorSettingsTermination
   propagatorSettingsDependentVariables
   propagatorSettingsCoordinates
   simulatorCreation
   TimeStateTemplates
   variationalSimulatorCreation
