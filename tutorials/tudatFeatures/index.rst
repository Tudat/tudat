.. _tudatFeaturesIndex:

Tudat Libraries
===============

These pages of the wiki provide further details about critical libraries necessary when setting up a simulation. A graphical representation is shown below, linked to the corresponding section of the Tudat libraries. Below the diagram the list of content are presented in the preferred order of reading and they follow how they should be typically be added to the simulation source file.

.. graphviz::

   digraph 
   {
      # general graph settings    
      rankdir = "LR";
      splines = ortho;
      

      # general node settings 
      node [shape = box, style = filled, width = 1.25, fixedsize = true, color = lightgrey, fontname = FontAwesome, fontsize = 8];


      # specific node color settings
      DynamicsSimulator [color = lightgreen];
      IntegratorSettings, NamedBodyMap, PropagatorSettings [color = lightblue];


      # Hyperlinks (Sphinx auto referencing not working here, need to link to exact web adres)
      DynamicsSimulator [href = "http://tudat.tudelft.nl/tutorials/tudatFeatures/propagationSetup/simulatorCreation.html#DynamicsSimulator", target = "_top"]; 
      IntegratorSettings [href = "http://tudat.tudelft.nl/tutorials/tudatFeatures/propagationSetup/integratorSettings.html#IntegratorSettings", target = "_top"];
      NamedBodyMap [href = "http://tudat.tudelft.nl/tutorials/tudatFeatures/environmentSetup/index.html#NamedBodyMap", target = "_top"];
      PropagatorSettings [href = "http://tudat.tudelft.nl/tutorials/tudatFeatures/propagationSetup/propagatorSettings.html#PropagatorSettings", target = "_top"];
      AccelerationModel [href = "http://tudat.tudelft.nl/tutorials/tudatFeatures/accelerationSetup/frameworkAcceleration.html#SelectedAccelerationMap", target = "_top"];
      BodySettings [href = "http://tudat.tudelft.nl/tutorials/tudatFeatures/environmentSetup/index.html#BodySettings", target = "_top"];
      AccelerationSettings [href = "http://tudat.tudelft.nl/tutorials/tudatFeatures/accelerationSetup/frameworkAcceleration.html#AccelerationSettings", target = "_top"];
      TerminationSettings [href = "http://tudat.tudelft.nl/tutorials/tudatFeatures/propagationSetup/propagatorSettingsTermination.html#PropagationTerminationSettings", target = "_top"];
      

      # DynamicsSimulatorInput
      NamedBodyMap -> DynamicsSimulator;
      PropagatorSettings -> DynamicsSimulator;
      IntegratorSettings -> DynamicsSimulator;


      # PropagatorSettings input     
      NamedBodyMap -> PropagatorSettings;
      AccelerationModel -> PropagatorSettings;
      initialBodyStates -> PropagatorSettings;
      bodiesToPropagate -> PropagatorSettings;
      TerminationSettings -> PropagatorSettings;
      centralBodies -> PropagatorSettings;
      dependentVariables -> PropagatorSettings;

      subgraph clusterPropagatorSettings
      {   
         style = dashed;
         dependentVariables, TerminationSettings, initialBodyStates;
         dependentVariables [style = dotted , fillcolor = lightgrey, color = black];
      }


      # BodySettings input
      "User-defined \nenvironment settings" -> BodySettings;
      getDefaultSettings -> BodySettings;
      bodyNames -> getDefaultSettings;


      # NamedBodyMap input
      BodySettings -> NamedBodyMap;

      subgraph clusterNamedBodyMap
      {   
         style = dashed;
         "User-defined \nenvironment settings" [style = dotted , fillcolor = lightgrey, color = black];
         { rank = same; "User-defined \nenvironment settings", BodySettings } 
         { rank = same; bodyNames, getDefaultSettings } 
      }
     

      # AccelerationSettings input
      accelerationType -> AccelerationSettings;
      "Additional \ninformation" -> AccelerationSettings;

      subgraph clusterAccelerationSettings
      {
         style = dashed;
         "Additional \ninformation",accelerationType;
         "Additional \ninformation" [style = dotted, color = black];
      } 


      # AccelerationModel input
      AccelerationSettings -> AccelerationModel;
      bodiesToPropagate -> AccelerationModel;
      centralBodies -> AccelerationModel;
      NamedBodyMap -> AccelerationModel;

      # IntegratorSettings input
      integratorType -> IntegratorSettings;
      simulationStartEpoch -> IntegratorSettings;
      "(initial) stepSize" -> IntegratorSettings;
     
      subgraph clusterIntegratorSettings
      {
         style = dashed;
         integratorType, simulationStartEpoch, "(initial) stepSize";
      } 
               
      # Define rank of some blocks for nice alignment 
      { rank = same; PropagatorSettings, IntegratorSettings, NamedBodyMap, AccelerationModel }
   }

.. toctree::
   :numbered:
   :maxdepth: 2

   astroTools/index
   mathTools/index
   environmentSetup/index
   accelerationSetup/index
   propagationSetup/index
   estimationSetup/index
   otherLibraries/index
